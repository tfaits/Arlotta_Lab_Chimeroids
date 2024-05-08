## Monocle runs on Chimeroid datasets
library(monocle3)
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggridges)
library(plyr)
library(patchwork)
library(magrittr)
library(circlize)

TypeCols <- c("#8dd3c7", "#fff7bc", "#bebada", "#fcbba1", "#80b1d3", "#fdb462",
              "#fa9fb5", "#04cf18", "#a6d854", "#fccde5", "#d9d9d9")
names(TypeCols) <- c("aRG","oRG","IP","PN","CFuPN","CPN",
                     "IN","ChPL","CR","Cortical hem","Unknown")
DonorCols <- c("darkred","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702")
names(DonorCols) <- c("11a","CW","GM","H1","PGP1","Mito210")
ProtCols <- c()
ProtCols["Single"] <- "#8cdb86"
ProtCols["Multi"] <- "#ce367c"
ProtCols["SingleDonor"] <- ProtCols["Single"]
ProtCols["MultiDonor"] <- ProtCols["Multi"]
ProtCols["SingleDonorNSC"] <- ProtCols["Single"]
ProtCols["MultiDonorNSC"] <- ProtCols["Multi"]
ProtCols["NPC"] <- "#cd8432"
ProtCols["Atlas"] <- "#196099"
ProtCols["Velasco"] <- "#196099"
col_fun = colorRamp2(c(-1, 0, 1), c("#b95858", "white", "#648ab1"))
`%ni%` = Negate(`%in%`)
                                 

## This is a test
CheckPackage <- function(package, repository, ...) {
  if (!requireNamespace(package = basename(path = package), quietly = TRUE)) {
    if (interactive()) {
      message("Package ", package, " is not yet installed")
      message("Install now?")
      choice <- menu(choices = c('yes', 'no'))
      if (choice == 1) {
        repository <- match.arg(
          arg = tolower(x = repository),
          choices = c('github', 'bioconductor', 'cran')
        )
        switch(
          EXPR = repository,
          'github' = remotes::install_github(repo = package, ...),
          'bioconductor' = BiocManager::install(pkgs = package, ...),
          'cran' = install.packages(pkgs = package, ...),
          stop("Unknown repository ", repository, call. = FALSE)
        )
        return(invisible(x = NULL))
      }
    }
    stop("Unable to find package ", package, ", please install", call. = FALSE)
  }
}
AssociatedDimReducs <- function(
    object,
    assay = DefaultAssay(object = object),
    global = TRUE
) {
  return(Filter(
    f = function(x) {
      check <- DefaultAssay(object = object[[x]]) == assay
      if (global) {
        check <- c(check, IsGlobal(object = object[[x]]))
      }
      return(any(check))
    },
    x = Reductions(object = object)
  ))
}

as.cell_data_set.Seurat <- function(
    x,
    assay = DefaultAssay(object = x),
    reductions = AssociatedDimReducs(object = x, assay = assay),
    default.reduction = DefaultDimReduc(object = x, assay = assay),
    graph = paste0(assay, '_snn'),
    group.by = NULL,
    ...
) {
  cds <- as(
    object = as.SingleCellExperiment(x = x, assay = assay),
    Class = 'cell_data_set'
  )
  if (is.null(x = SummarizedExperiment::assays(x = cds)$counts)) {
    SummarizedExperiment::assays(x = cds)$counts <- SummarizedExperiment::assays(x = cds)[[1]]
  }
  if (!"Size_Factor" %in% colnames(x = SummarizedExperiment::colData(x = cds))) {
    size.factor <- paste0('nCount_', assay)
    if (size.factor %in% colnames(x = x[[]])) {
      SummarizedExperiment::colData(x = cds)$Size_Factor <- x[[size.factor, drop = TRUE]]
    }
  }
  SingleCellExperiment::reducedDims(x = cds)[SingleCellExperiment::reducedDimNames(x = cds)] <- NULL
  reductions <- intersect(
    x = reductions,
    y = AssociatedDimReducs(object = x, assay = assay)
  )
  for (reduc in reductions) {
    SingleCellExperiment::reducedDims(x = cds)[[toupper(x = reduc)]] <- Embeddings(object = x[[reduc]])
    loadings <- Loadings(object = x[[reduc]])
    if (!IsMatrixEmpty(x = loadings)) {
      slot(object = cds, name = 'reduce_dim_aux')[['gene_loadings']] <- loadings
    }
    stdev <- Stdev(object = x[[reduc]])
    if (length(x = stdev)) {
      slot(object = cds, name = 'reduce_dim_aux')[['prop_var_expl']] <- stdev
    }
  }
  if (!is.null(x = group.by)) {
    Idents(object = x) <- group.by
  }
  clusters.list <- list()
  if (length(x = clusters.list)) {
    slot(object = cds, name = 'clusters')[[toupper(x = default.reduction)]] <- clusters.list
  }
  return(cds)
}


RunMonocle <- function(seur_obj, root) {
  cds = as.cell_data_set.Seurat(seur_obj)
  cds = cluster_cells(cds, reduction_method="UMAP")
  cds = learn_graph(cds, use_partition=T)
  print(paste0("using the cluster ", root, " as the root_cells"))
  #start_cells = WhichCells(seur_obj, expression=CellType==root)
  start_cells = WhichCells(seur_obj, idents=root)
  cds = order_cells(cds, reduction_method="UMAP", root_cells=start_cells)
  return(cds)
}

seur <- readRDS("NSC_SD-MD_Control-VPA-EtOH.rds")
Idents(seur) <- seur$MergedAnnotation
MD_Con <- subset(seur, subset=Protocol=="MultiDonorNSC" & Treatment=="Control")
SD_Con <- subset(seur, subset=Protocol=="SingleDonorNSC" & Treatment=="Control")
MD_VPA <- subset(seur, subset=Protocol=="MultiDonorNSC" & Treatment=="VPA")
SD_VPA <- subset(seur, subset=Protocol=="SingleDonorNSC" & Treatment=="VPA")
rm(seur)

mdc <- RunMonocle(MD_Con, "aRG")
MD_Con <- AddMetaData(object=MD_Con, metadata=mdc@principal_graph_aux@listData$UMAP$pseudotime, col.name="pseudotime")
MD_Con$pseudotime[MD_Con$pseudotime > 25] <- 25
FeaturePlot(MD_Con, features="pseudotime", reduction='umap')
plot_cells(mdc, color_cells_by = "MergedAnnotation", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
ggplot(MD_Con@meta.data, aes(x=pseudotime, color=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density(alpha=0.6, bw=1) +
  scale_fill_manual(values=TypeCols[unique(MD_Con$MergedAnnotation)]) +
  scale_color_manual(values=TypeCols[unique(MD_Con$MergedAnnotation)]) +
  ggtitle("Multidonor Controls") +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey85", linewidth=0.3),
        plot.title = element_text(face="bold", hjust=0.5))
MD_Con <- MD_Con@meta.data

sdc <- RunMonocle(SD_Con, "aRG")
SD_Con <- AddMetaData(object=SD_Con, metadata=sdc@principal_graph_aux@listData$UMAP$pseudotime, col.name="pseudotime")
SD_Con$pseudotime[SD_Con$pseudotime == Inf] <- 25
FeaturePlot(SD_Con, features="pseudotime", reduction='umap')
plot_cells(sdc, color_cells_by = "MergedAnnotation", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
ggplot(SD_Con@meta.data, aes(x=pseudotime, color=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density(alpha=0.6, bw=1) +
  scale_fill_manual(values=TypeCols[unique(SD_Con$MergedAnnotation)]) +
  scale_color_manual(values=TypeCols[unique(SD_Con$MergedAnnotation)]) +
  ggtitle("Single-Donor Controls") +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey85", linewidth=0.3),
        plot.title = element_text(face="bold", hjust=0.5))
SD_Con <- SD_Con@meta.data

tmp <- rbind(MD_Con, SD_Con)
ggplot(tmp, aes(x=pseudotime, color=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density(alpha=0.6, bw=1) +
  scale_fill_manual(values=TypeCols[unique(SD_Con$MergedAnnotation)]) +
  scale_color_manual(values=TypeCols[unique(SD_Con$MergedAnnotation)]) +
  ggtitle("Single-Donor Controls") +
  facet_grid(Protocol~1) +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey85", linewidth=0.3),
        plot.title = element_text(face="bold", hjust=0.5))

ggplot(tmp, aes(x=pseudotime, color=Genotype, fill=Genotype)) +
  geom_density(alpha=0.6, bw=1) +
  scale_fill_manual(values=DonorCols[unique(SD_Con$Genotype)]) +
  scale_color_manual(values=DonorCols[unique(SD_Con$Genotype)]) +
  ggtitle("Single-Donor Controls") +
  facet_grid(Protocol~1) +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey85", linewidth=0.3),
        plot.title = element_text(face="bold", hjust=0.5))

atlas <- readRDS("Atlas_3mo.rds")
atlas <- subset(atlas, subset=dataset %in% c("4","5","6","7"))
Idents(atlas) <- atlas$FinalName
adc <- RunMonocle(atlas, "aRG")
atlas <- AddMetaData(object=atlas, metadata=adc@principal_graph_aux@listData$UMAP$pseudotime, col.name="pseudotime")
FeaturePlot(atlas, features="pseudotime", reduction='umap')
plot_cells(adc, color_cells_by = "FinalName", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)

atlas <- atlas@meta.data
atlas$MergedAnnotation <- atlas$FinalName
atlas$Protocol <- "Velasco"

SD_Con$ScaledPseudotime <- SD_Con$pseudotime/max(SD_Con$pseudotime)
MD_Con$ScaledPseudotime <- MD_Con$pseudotime/max(MD_Con$pseudotime)
atlas$ScaledPseudotime <- atlas$pseudotime/max(atlas$pseudotime)

npc <- readRDS("M3M4_NPC_NSC.rds")
npc <- subset(npc, subset=orig.ident %in% c("August2022_Mix1_D23_1", "August2022_Mix1_D23_2"))
Idents(npc) <- npc$MergedAnnotation
NPC <- RunMonocle(npc, "aRG")
npc <- AddMetaData(object=npc, metadata=NPC@principal_graph_aux@listData$UMAP$pseudotime, col.name="pseudotime")
FeaturePlot(npc, features="pseudotime", reduction='umap')
plot_cells(NPC, color_cells_by = "MergedAnnotation", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)
npc <- npc@meta.data
npc$Protocol <- "NPC"

npc$ScaledPseudotime <- npc$pseudotime/max(npc$pseudotime)
atlas$Genotype <- sapply(atlas$org, function(x){strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1]})

tmp <- rbind(atlas[,c("Protocol","MergedAnnotation","Genotype","ScaledPseudotime")],
             MD_Con[,c("Protocol","MergedAnnotation","Genotype","ScaledPseudotime")],
             npc[,c("Protocol","MergedAnnotation","Genotype","ScaledPseudotime")],
             SD_Con[,c("Protocol","MergedAnnotation","Genotype","ScaledPseudotime")])
tmp$MergedAnnotation[tmp$MergedAnnotation=="Immature IN"] <- "IN"
tmp <- tmp[tmp$MergedAnnotation %in% c("aRG","oRG","IP","PN","CFuPN","CPN","IN"),]
tmp$MergedAnnotation <- factor(tmp$MergedAnnotation, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN"))
tmp$Protocol <- factor(tmp$Protocol, levels=c("Velasco","SingleDonorNSC","MultiDonorNSC","NPC"))
ggplot(tmp, aes(x=ScaledPseudotime, color=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density(alpha=0.6, bw=0.04) +
  scale_fill_manual(values=TypeCols[levels(tmp$MergedAnnotation)]) +
  scale_color_manual(values=TypeCols[levels(tmp$MergedAnnotation)]) +
  ggtitle("Pseudotime By Cell Type") +
  facet_grid(Protocol~1) +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey85", linewidth=0.3),
        plot.title = element_text(face="bold", hjust=0.5))

ggplot(tmp, aes(x=ScaledPseudotime, y=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density_ridges(size=0.3) + 
  scale_fill_manual(values=TypeCols[unique(tmp$MergedAnnotation)], name="Cell Type") + 
  theme_classic() +
  labs(x = "Scaled Pseudotime", y="Density of Cells per Cell Type") +
  facet_grid(~Protocol)

ggplot(MD_Con[MD_Con$MergedAnnotation %in% c("aRG","oRG","IP","CFuPN","CPN","PN","IN","ChPL","CR"),], aes(x=ScaledPseudotime, y=MergedAnnotation, fill=MergedAnnotation)) +
  geom_density_ridges(size=0.3) + 
  scale_fill_manual(values=TypeCols[unique(MD_Con$MergedAnnotation)], name="Cell Type") + 
  theme_classic() +
  labs(x = "Scaled Pseudotime", y="Density of Cells per Cell Type") +
  facet_grid(~Genotype)

tmp$Genotype <- factor(tmp$Genotype, levels=rev(c("CW","H1","Mito210","PGP1","11a","GM")))

ggplot(tmp, aes(x=ScaledPseudotime, y=Genotype, fill=Genotype)) +
  geom_density_ridges(size=0.3) + 
  scale_fill_manual(values=DonorCols[levels(tmp$Genotype)], name="Donor") + 
  theme_classic() +
  labs(x = "Scaled Pseudotime", y="Density of Cells per Donor") +
  facet_grid(~Protocol)



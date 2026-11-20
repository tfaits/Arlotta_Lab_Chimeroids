## Pseudobulk, but pulling in the Atlas (PGP1 and Mito210) and NPC orgs?
library(Seurat)
library(ggplot2)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(ComplexHeatmap)
library(svglite)
library(circlize)
library(openxlsx)
library(reshape2)

setwd("/Users/tfaits/Documents/Chimera/OctoberSequencing2023")

TypeCols <- c("#8dd3c7", "#fff7bc", "#bebada", "#fcbba1", "#80b1d3", "#fdb462",
              "#fa9fb5", "#04cf18", "#a6d854", "#fccde5", "#d9d9d9")
names(TypeCols) <- c("aRG","oRG","IP","PN","CFuPN","CPN",
                      "IN","ChPL","CR","Cortical hem","Unknown")
DonorCols <- c("darkred","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702")
names(DonorCols) <- c("11a","CW","GM","H1","PGP1","Mito210")
DonorCols["Mito"] <- DonorCols["Mito210"]
DonorCols["PGP"] <- DonorCols["PGP1"]
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

# # Load in the Atlas Mito210 and PGP1 data:
# atlas <- readRDS("/Users/tfaits/Desktop/Documents/Atlas_RNA_data/ForLabelTransfer/3mo_harm_111120.rds")
# atlas <- subset(atlas, subset=dataset %in% c("4","5","6","7"))
# atlas$Pseudobulk_CellType <- paste(atlas$org, atlas$FinalName, sep="-")
# ## For the cell-type split, we'll want to remove any with fewer than 20 cells.
# allGroups <- table(atlas$Pseudobulk_CellType)
# goodGroups <- names(allGroups)[allGroups >= 20]
# print(paste(length(allGroups), length(goodGroups))) #Looks like 77 out of 79 samples make it?
# pbulk2 <- matrix(ncol=length(goodGroups), nrow=nrow(atlas@assays$RNA@counts))
# colnames(pbulk2) <- goodGroups
# rownames(pbulk2) <- rownames(atlas@assays$RNA@counts)
# idx <- 0
# for(i in goodGroups){
#   idx <- idx+1
#   pbulk2[,i] <- rowSums(atlas@assays$RNA@counts[,atlas$Pseudobulk_CellType==i])
#   print(idx)
# }
# saveRDS(pbulk2, "PseudoBulk_CellType_AtlasOnly.RDS")
# metabulk2 <- data.frame(Sample=colnames(pbulk2),
#                        Protocol="Velasco",
#                        Genotype=sapply(colnames(pbulk2), function(x){strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1]}),
#                        CellType=sapply(colnames(pbulk2), function(x){strsplit(x,"-")[[1]][2]}),
#                        Treatment="Control",
#                        Mix="Atlas",
#                        Batch=sapply(colnames(pbulk2), function(x){paste(strsplit(x,"_")[[1]][1:(length(strsplit(x,"_")[[1]])-1)], collapse="_")}),
#                        Replicate=sapply(colnames(pbulk2), function(x){strsplit(strsplit(x,"_")[[1]][length(strsplit(x,"_")[[1]])],"-")[[1]][1]}),
#                        Organoid=sapply(colnames(pbulk2), function(x){strsplit(x,"-")[[1]][1]}),
#                        Individual=sapply(colnames(pbulk2), function(x){strsplit(x,"-")[[1]][1]}),
#                        ReadCount=colSums(pbulk2),
#                        CellCount=sapply(colnames(pbulk2), function(x){sum(atlas$Pseudobulk_CellType==x)}))
# saveRDS(metabulk2,"PseudoBulk_CellType_Metadata_AtlasOnly.RDS")


# # Load in the Atlas Mito210 and PGP1 data:
# npc <- readRDS("../Seurats/Finalized_EtOH_VPA_Objects/Mix3-Mix4_D15-D23_Integrated.RDS")
# npcMeta <- readRDS("../Seurats/Finalized_EtOH_VPA_Objects/D23_D15_M3M4_Trimmed_Metadata.RDS")
# npcMeta <- npcMeta[npcMeta$Sample %in% c("Mix3 D23 1", "Mix3 D23 2"),]
# npc <- subset(npc, cells=which(colnames(npc) %in% rownames(npcMeta)))
# npc@meta.data <- npcMeta
# head(npcMeta)
# #Batch, Mix, Treatment, Replicate, Genotype
# npc$Individual <- paste("NPC-A", "Mix3", "Control", sapply(npc$orig.ident, function(x){strsplit(x,"_")[[1]][4]}),npc$Genotype, sep="_")
# npc$Pseudobulk_CellType <- paste(npc$Individual, npc$CellType, sep="_")
# ## For the cell-type split, we'll want to remove any with fewer than 20 cells.
# allGroups <- table(npc$Pseudobulk_CellType)
# goodGroups <- names(allGroups)[allGroups >= 20]
# print(paste(length(allGroups), length(goodGroups))) #Looks like 40 out of 54 samples make it
# pbulk3 <- matrix(ncol=length(goodGroups), nrow=nrow(npc@assays$RNA@counts))
# colnames(pbulk3) <- goodGroups
# rownames(pbulk3) <- rownames(npc@assays$RNA@counts)
# idx <- 0
# for(i in goodGroups){
#   idx <- idx+1
#   pbulk3[,i] <- rowSums(npc@assays$RNA@counts[,npc$Pseudobulk_CellType==i])
#   print(idx)
# }
# saveRDS(pbulk3, "PseudoBulk_CellType_NPCOnly.RDS")
# metabulk3 <- data.frame(Sample=colnames(pbulk3),
#                         Protocol="NPC",
#                         Genotype=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][5]}),
#                         CellType=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][6]}),
#                         Treatment=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][3]}),
#                         Mix=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][2]}),
#                         Batch=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][1]}),
#                         Replicate=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][4]}),
#                         Organoid=sapply(colnames(pbulk3), function(x){paste(strsplit(x,"_")[[1]][c(1:4)], collapse="_")}),
#                         Individual=sapply(colnames(pbulk3), function(x){paste(strsplit(x,"_")[[1]][c(1:5)], collapse="_")}),
#                         ReadCount=colSums(pbulk3),
#                         CellCount=sapply(colnames(pbulk3), function(x){sum(seur$Pseudobulk_CellType==x)}))
# saveRDS(metabulk3,"PseudoBulk_CellType_Metadata_NPCOnly.RDS")

# npc <- readRDS("../Seurats/Finalized_EtOH_VPA_Objects/Mix3-Mix4_D15-D23_Integrated.RDS")
# npcMeta <- readRDS("../Seurats/Finalized_EtOH_VPA_Objects/D23_D15_M3M4_Trimmed_Metadata.RDS")
# npcMeta <- npcMeta[npcMeta$Sample %in% c("Mix3 D23 1", "Mix3 D23 2"),]
# npc <- subset(npc, cells=which(colnames(npc) %in% rownames(npcMeta)))
# npc@meta.data <- npcMeta
# head(npcMeta)
# #Batch, Mix, Treatment, Replicate, Genotype
# npc$Individual <- paste("NPC-A", "Mix3", "Control", sapply(npc$orig.ident, function(x){strsplit(x,"_")[[1]][4]}),npc$Genotype, sep="_")
# ## For the cell-type split, we'll want to remove any with fewer than 20 cells.
# allGroups <- table(npc$Individual)
# goodGroups <- names(allGroups)[allGroups >= 20]
# print(paste(length(allGroups), length(goodGroups))) #Looks like 40 out of 54 samples make it
# pbulk3 <- matrix(ncol=length(goodGroups), nrow=nrow(npc@assays$RNA@counts))
# colnames(pbulk3) <- goodGroups
# rownames(pbulk3) <- rownames(npc@assays$RNA@counts)
# idx <- 0
# for(i in goodGroups){
#   idx <- idx+1
#   pbulk3[,i] <- rowSums(npc@assays$RNA@counts[,npc$Individual==i])
#   print(idx)
# }
# saveRDS(pbulk3, "GlobalPseudoBulk_CellType_NPCOnly.RDS")
# metabulk3 <- data.frame(Sample=colnames(pbulk3),
#                         Protocol="NPC",
#                         Genotype=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][5]}),
#                         Treatment=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][3]}),
#                         Mix=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][2]}),
#                         Batch=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][1]}),
#                         Replicate=sapply(colnames(pbulk3), function(x){strsplit(x,"_")[[1]][4]}),
#                         Organoid=sapply(colnames(pbulk3), function(x){paste(strsplit(x,"_")[[1]][c(1:4)], collapse="_")}),
#                         Individual=sapply(colnames(pbulk3), function(x){paste(strsplit(x,"_")[[1]][c(1:5)], collapse="_")}),
#                         ReadCount=colSums(pbulk3),
#                         CellCount=sapply(colnames(pbulk3), function(x){sum(npc$Individual==x)}))
# saveRDS(metabulk3,"GlobalPseudoBulk_CellType_Metadata_NPCOnly.RDS")

## Load in the Pseudobulk matrix and metadata for the NSC Chimeroids (Single-donor and Multi-donor)
# pbulk <- readRDS("PseudoBulk_CellType_UsingAnnotationFromNovember17b.RDS")
# pbulk2 <- readRDS("PseudoBulk_CellType_AtlasOnly.RDS")
# pbulk3 <- readRDS("PseudoBulk_CellType_NPCOnly.RDS")
# metabulk <- readRDS("PseudoBulk_CellType_Metadata_UsingAnnotationFromNovember17b.RDS")
# metabulk$Genotype[metabulk$Genotype=="CIRM17"] <- "CW"
# metabulk$Mix[metabulk$Mix=="CIRM17"] <- "CW"
# metabulk2 <- readRDS("PseudoBulk_CellType_Metadata_AtlasOnly.RDS")
# metabulk3 <- readRDS("PseudoBulk_CellType_Metadata_NPCOnly.RDS")
# 
# sharedGenes <- rownames(pbulk)[rownames(pbulk) %in% rownames(pbulk2) & rownames(pbulk) %in% rownames(pbulk3)]
# print(paste(nrow(pbulk), nrow(pbulk2), nrow(pbulk3), length(sharedGenes))) # Only 19,601 gene shared between all three datasets?
# 
# pbulk <- cbind(pbulk[sharedGenes,], pbulk2[sharedGenes,], pbulk3[sharedGenes,])
# metabulk <- rbind(metabulk, metabulk2, metabulk3)
# saveRDS(pbulk, "PseudBulk_FourProtocols.RDS")
# saveRDS(metabulk, "PseudoBulk_Metadata_FourProtocols.RDS")
pbulk <- readRDS("PseudBulk_FourProtocols.RDS")
metabulk <- readRDS("PseudoBulk_Metadata_FourProtocols.RDS")
metabulk$Batch[metabulk$Batch=="SingleDonorA"] <- "BatchB"
metabulk$Batch[metabulk$Batch=="SingleDonorB"] <- "BatchC"

## Functions!
### Subset the pseudobulk dataset based on user-selected traits.
subsetPseudobulk <- function(Treatments=c("Control","VPA"),
                             CellTypes=c("aRG"),
                             Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                             Protocols=c("MultiDonor"),
                             Mixes=c("Mix1","Mix2","CW","H1","Mito210","PGP1","11a","Mix3","Atlas")){
  tmpMeta <- metabulk[metabulk$Treatment %in% Treatments &
                        metabulk$CellType %in% CellTypes &
                        metabulk$Genotype %in% Genotypes &
                        metabulk$Protocol %in% Protocols &
                        metabulk$Mix %in% Mixes,]
  tmpPseudo <- pbulk[,metabulk$Treatment %in% Treatments &
                       metabulk$CellType %in% CellTypes &
                       metabulk$Genotype %in% Genotypes &
                       metabulk$Protocol %in% Protocols &
                       metabulk$Mix %in% Mixes]
  tmpPseudo <- tmpPseudo[rowSums(tmpPseudo) > 10,]
  return(list(PB=tmpPseudo, META=tmpMeta))
}

## Global:
gPB <- readRDS("PseudoBulk_Global_UsingAnnotationFromNovember17b.RDS")
gMeta <- readRDS("PseudoBulk_Global_Metadata_UsingAnnotationFromNovember17b.RDS")
gMeta$Genotype[gMeta$Genotype=="CIRM17"] <- "CW"
gMeta$Mix[gMeta$Mix=="CIRM17"] <- "CW"
pbulk3 <- readRDS("GlobalPseudoBulk_CellType_NPCOnly.RDS")
metabulk3 <- readRDS("GlobalPseudoBulk_CellType_Metadata_NPCOnly.RDS")
sGs <- rownames(gPB)[rownames(gPB) %in% rownames(pbulk3)]
gPB <- cbind(gPB[sGs,], pbulk3[sGs,])
gMeta <- rbind(gMeta, metabulk3)
saveRDS(gPB, "PseudoBulk_Global_UsingAnnotationFromNovember17c.RDS")
saveRDS(gMeta, "PseudoBulk_Global_Metadata_UsingAnnotationFromNovember17c.RDS")
gMeta$Batch[gMeta$Batch=="SingleDonorA"] <- "BatchB"
gMeta$Batch[gMeta$Batch=="SingleDonorB"] <- "BatchC"

subsetGlobalPseudobulk <- function(Treatments=c("Control","VPA"),
                                   Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                   Protocols=c("MultiDonor"),
                                   Mixes=c("Mix1","Mix2","CW","H1","Mito210","PGP1","11a")){
  tmpMeta <- gMeta[gMeta$Treatment %in% Treatments &
                     gMeta$Genotype %in% Genotypes &
                     gMeta$Protocol %in% Protocols &
                     gMeta$Mix %in% Mixes,]
  tmpPseudo <- gPB[,gMeta$Treatment %in% Treatments &
                     gMeta$Genotype %in% Genotypes &
                     gMeta$Protocol %in% Protocols &
                     gMeta$Mix %in% Mixes]
  tmpPseudo <- tmpPseudo[rowSums(tmpPseudo) > 10,]
  return(list(PB=tmpPseudo, META=tmpMeta))
}

### Run DESeq2 on a properly subsetted pseudobulk object. Can be customized for any design matrix you want.
runDESeq <- function(tmpPseudo, tmpMeta, # The properly filtered/subsetted pseudobulk matrix
                     myEffect="Protocol",
                     addCovs=""){
  print("Run DESeq2")
  design <- as.formula(paste("~",myEffect, sep=""))
  if(nchar(addCovs)>0){
    design <- as.formula(paste("~",addCovs," + ", myEffect, sep=""))
  }
  print(design)
  dds <- DESeqDataSetFromMatrix(countData=tmpPseudo, colData=tmpMeta, design=design)
  dds <- DESeq(dds)
  resNam <- resultsNames(dds)
  res <- lfcShrink(dds, coef=length(resNam))
  return(as.data.frame(res))
}

### Get the DESeq results table and also run a full GO-term overrepresentation analysis:
GetDE_Results <- function(Treatments=c("Control","VPA"),
                          CellTypes=c("aRG"),
                          Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                          Protocols=c("MultiDonor"),
                          Mixes=c("Mix1","Mix2","CW","H1","Mito210","PGP1","11a","Mix3","Atlas"),
                          myEffect="Protocol",
                          addCovs="",
                          IsGlobal=FALSE){
  if(IsGlobal){
    myPB <- subsetGlobalPseudobulk(Treatments=Treatments,
                                   Genotypes=Genotypes,
                                   Protocols=Protocols,
                                   Mixes=Mixes)
  }else{
    myPB <- subsetPseudobulk(Treatments=Treatments,
                             CellTypes=CellTypes,
                             Genotypes=Genotypes,
                             Protocols=Protocols,
                             Mixes=Mixes)
  }
  print(table(myPB$META$Protocol, myPB$META$Genotype))
  myRes <- runDESeq(myPB$PB,
                    myPB$META,
                    myEffect=myEffect,
                    addCovs=addCovs)
  SigRes <- myRes[!is.na(myRes$padj) & myRes$padj < 0.1,]
  UpGenes <- rownames(SigRes)[SigRes$log2FoldChange > 0.5]
  DownGenes <- rownames(SigRes)[SigRes$log2FoldChange < -0.5]
  upGO <- enrichGO(gene=UpGenes,
                   OrgDb = org.Hs.eg.db,
                   keyType="SYMBOL",
                   ont="BP",
                   universe=rownames(myPB$PB))
  downGO <- enrichGO(gene=DownGenes,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     ont="BP",
                     universe=rownames(myPB$PB))
  tryCatch( { upSimp <- simplify(upGO, cutoff=0.4) },
            warning = function(w) {upSimp <<- "None"},
            error = function(e) {upSimp <<- "None"})
  tryCatch( { downSimp <- simplify(downGO, cutoff=0.4) },
            warning = function(w) {downSimp <<- "None"},
            error = function(e) {downSimp <<- "None"})
  return(list(Results=as.data.frame(myRes), UpGo=upGO, DownGo=downGO, SimpUp=upSimp, SimpDown=downSimp))
}

### Draw a heatmap of a subsetted pseudobulk dataset; can select genes as desired.
Make_Heatmap <- function(myGenes=c("SATB2","FEZF2","SOX2","EOMES","DLX5"),
                         Treatments=c("Control","VPA"),
                         CellTypes=c("aRG"),
                         Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                         Protocols=c("MultiDonor"),
                         Mixes=c("Mix1","Mix2","CW","H1","Mito210","PGP1","11a","Mix3","Atlas"),
                         myTitle="Test Title",
                         IsGlobal=FALSE,
                         textSize=6,
                         scaleData=TRUE){
  if(IsGlobal){
    myPB <- subsetGlobalPseudobulk(Treatments=Treatments,
                                   Genotypes=Genotypes,
                                   Protocols=Protocols,
                                   Mixes=Mixes)
    topAnnotation <- HeatmapAnnotation(Treatment=myPB$META$Treatment,
                                       Genotype=myPB$META$Genotype,
                                       Batch=myPB$META$Batch,
                                       Protocol=myPB$META$Protocol,
                                       col = list(Treatment=c("VPA"="#F2CC54", "Control"="#EAEAEA", "EtOH"="#B497D4"),
                                                  Genotype=c("CW"="#0a9396", "GM"="#94d2bd", "H1"="#E9d8a6", "PGP1"="#ee9b00", "Mito210"="#ca6702","11a"="darkred"),
                                                  Batch=c("BatchA"="blue","BatchB"="green","BatchC"="red","SingleDonorA"="darkblue","SingleDonorB"="darkgreen", "NPC-A"="darkred","DorsalKadoshima.Mito210_3m"="#ba5700", "DorsalKadoshima.Mito210_3m_PTEN"="#904000","DorsalKadoshima.PGP1_b1_3mon"="#ffab10","DorsalKadoshima.PGP1_b2_3mon"="#dd8b00"),
                                                  Protocol=c("SingleDonor"="#9cbd86", "MultiDonor"="#ce367c","Velasco"="darkblue",NPC="purple")))
  }else{
    myPB <- subsetPseudobulk(Treatments=Treatments,
                             CellTypes=CellTypes,
                             Genotypes=Genotypes,
                             Protocols=Protocols,
                             Mixes=Mixes)
    topAnnotation <- HeatmapAnnotation(Treatment=myPB$META$Treatment,
                                       Genotype=myPB$META$Genotype,
                                       Batch=myPB$META$Batch,
                                       Protocol=myPB$META$Protocol,
                                       CellTypes=myPB$META$CellType,
                                       col = list(Treatment=c("VPA"="#F2CC54", "Control"="#EAEAEA", "EtOH"="#B497D4"),
                                                  Genotype=c("CW"="#0a9396", "GM"="#94d2bd", "H1"="#E9d8a6", "PGP1"="#ee9b00", "Mito210"="#ca6702","11a"="darkred"),
                                                  Batch=c("BatchA"="blue","BatchB"="green","BatchC"="red","SingleDonorA"="darkblue","SingleDonorB"="darkgreen", "NPC-A"="darkred","DorsalKadoshima.Mito210_3m"="#ba5700", "DorsalKadoshima.Mito210_3m_PTEN"="#904000","DorsalKadoshima.PGP1_b1_3mon"="#ffab10","DorsalKadoshima.PGP1_b2_3mon"="#dd8b00"),
                                                  Protocol=c("SingleDonor"="#9cbd86", "MultiDonor"="#ce367c","Velasco"="darkblue",NPC="purple"),
                                                  CellTypes=TypeCols))
  }
  myPB$PB <- log(t(t(myPB$PB)/colSums(myPB$PB)*1000000)+1)#Turn into logCPM
  if(scaleData==TRUE){
    myExp <- t(scale(t(myPB$PB[myGenes,])))
  }else{
    myExp <- myPB$PB[myGenes,]
  }
  myPlot <- Heatmap(myExp, top_annotation = topAnnotation,
                    show_column_names = FALSE, column_title = myTitle, row_names_gp = gpar(fontsize=textSize))
  return(myPlot)
}

### Get a list of the top upregulated and downregulated genes from a results table
getGenes <- function(inList, numGenes=100){
  return(c(rownames(inList$Results[!is.na(inList$Results$padj) & inList$Results$log2FoldChange > 0,])[order(inList$Results$padj[!is.na(inList$Results$padj) & inList$Results$log2FoldChange > 0])[1:min(floor(numGenes/2),sum(!is.na(inList$Results$padj) & inList$Results$log2FoldChange > 0))]],
           rownames(inList$Results[!is.na(inList$Results$padj) & inList$Results$log2FoldChange < 0,])[order(inList$Results$padj[!is.na(inList$Results$padj) & inList$Results$log2FoldChange < 0])[1:min(ceiling(numGenes/2),sum(!is.na(inList$Results$padj) & inList$Results$log2FoldChange > 0))]]))
}

getSigGenes <- function(inList){
  temp <- inList$Results[!is.na(inList$Results$padj) & inList$Results$padj < 0.1 & abs(inList$Results$log2FoldChange) > 0.5,]
  return(temp[order(temp$padj),])
}

## Find the correlation between logFCs between two DEG results tables
Correlate_Results <- function(inRes1, inRes2){
  inRes1 <- inRes1$Results
  inRes2 <- inRes2$Results
  inRes1 <- inRes1[!is.na(inRes1$pvalue),]
  inRes2 <- inRes2[!is.na(inRes2$pvalue),]
  SGs <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[SGs,]
  inRes2 <- inRes2[SGs,]
  return(cor(inRes1$log2FoldChange,
             inRes2$log2FoldChange))
}

Fisher_Results <- function(inRes1, inRes2){
  inRes1 <- inRes1$Results
  inRes2 <- inRes2$Results
  inRes1 <- inRes1[!is.na(inRes1$pvalue),]
  inRes2 <- inRes2[!is.na(inRes2$pvalue),]
  SGs <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[SGs,]
  inRes2 <- inRes2[SGs,]
  sharedSig <- sum(inRes1$padj < 0.1 & !is.na(inRes1$padj) &
                     inRes2$padj < 0.1 & !is.na(inRes2$padj) &
                     ((inRes1$log2FoldChange > 0.5 & inRes2$log2FoldChange > 0.5) |
                        (inRes1$log2FoldChange < -0.5 & inRes2$log2FoldChange < -0.5)))
  sig1 <- sum(inRes1$padj < 0.1 & !is.na(inRes1$padj) & abs(inRes1$log2FoldChange) > 0.5)
  sig2 <- sum(inRes2$padj < 0.1 & !is.na(inRes2$padj) & abs(inRes2$log2FoldChange) > 0.5)
  allG <- length(SGs)
  tmpMat <- matrix(ncol=2, nrow=2, data=c(sharedSig, sig1-sharedSig, sig2-sharedSig, allG+sharedSig-sig1-sig2))
  return(fisher.test(tmpMat)$p.value)
}

## Let's run some tests on some controls!
## Specifically, we first want to see whether the differences between protocols in PGP1 are the same differences...
## ...as those in Mito210 (between protocols).

myThing <- list()
for(myProtocol in c("NPC","SingleDonor")){
  for(myDonor in c("PGP1","Mito210")){
    print(paste("Global",myProtocol,myDonor, sep="-"))
    myThing[[paste("Global",myProtocol,myDonor, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                            IsGlobal=TRUE,
                                                                            Genotypes=c(myDonor),
                                                                            Protocols=c(myProtocol,"MultiDonor"),
                                                                            myEffect="Protocol")
  }
}

myPGP <- myThing$`Global-SingleDonor-PGP1`$Results
myPGP$Sig <- "Not Sig"
myPGP$Sig[myPGP$log2FoldChange > 0.5 & myPGP$padj < 0.05] <- "Up in SD"
myPGP$Sig[myPGP$log2FoldChange < -0.5 & myPGP$padj < 0.05] <- "Down in SD"

ggplot(myPGP, aes(x=log2FoldChange, y=-log10(pvalue), color=Sig)) +
  geom_point() +
  scale_color_manual(values=c("darkblue","grey80", "darkred")) +
  theme(panel.background = element_blank(),
        plot.title=element_text(face="bold", hjust=0.5)) +
  geom_abline(slope=0, intercept=3.3581, linetype="dashed") +
  geom_vline(xintercept=-0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  ggtitle("Global DEGs Between\nSingle-Donor and Multi-Donor NSC Chimeroids\nin PGP1 Cells")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_PGP1_SD_vs_MD.svg", height=4, width=6)


for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  for(myProtocol in c("Velasco","NPC","SingleDonor")){
    for(myDonor in c("PGP1","Mito210")){
      print(paste(CellType,myProtocol,myDonor, sep="-"))
      if(sum(metabulk$Protocol==myProtocol & metabulk$Genotype==myDonor & metabulk$Treatment=="Control" & metabulk$CellType==CellType) > 1){
        myThing[[paste(CellType,myProtocol,myDonor, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                                CellTypes=c(CellType),
                                                                                Genotypes=c(myDonor),
                                                                                Protocols=c(myProtocol,"MultiDonor"),
                                                                                myEffect="Protocol")
      }else{
        print("Not enough replicates.")
      }
    }
  }
}

myThing$Global_SingleDonor_All <- GetDE_Results(Treatments="Control",
                                                IsGlobal=TRUE,
                                                Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                Protocols=c("SingleDonor","MultiDonor"),
                                                myEffect="Protocol",
                                                addCovs = "Genotype + Batch")
for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  myThing[[paste0(CellType,"_SingleDonor_All")]] <- GetDE_Results(Treatments="Control",
                                                  IsGlobal=FALSE,
                                                  CellTypes = c(CellType),
                                                  Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                  Protocols=c("SingleDonor","MultiDonor"),
                                                  myEffect="Protocol",
                                                  addCovs = "Genotype + Batch")
}

SigColors <- c("darkblue","grey80","darkred")
names(SigColors) <- c("Down in SD", "Not Sig", "Up in SD")
makeVolcano <- function(inTab, myTitle="No Title"){
  inTab$Sig <- "Not Sig"
  inTab$Sig[inTab$log2FoldChange > 0.5 & inTab$padj < 0.05] <- "Up in SD"
  inTab$Sig[inTab$log2FoldChange < -0.5 & inTab$padj < 0.05] <- "Down in SD"
  myLabs <- c(paste0("\n% Down: ",round(sum(inTab$Sig=="Down in SD")/nrow(inTab)*100,digits=1)),
              paste0("% Up: ",round(sum(inTab$Sig=="Up in SD")/nrow(inTab)*100,digits=1)),
              paste0("\n%Not Sig: ",round(sum(inTab$Sig == "Not Sig")/nrow(inTab)*100,digits=1)))
  ggplot(inTab, aes(x=log2FoldChange, y=-log10(pvalue), color=Sig)) +
    geom_point(size=1, alpha=0.7) +
    scale_color_manual(values=SigColors[unique(inTab$Sig)]) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    coord_cartesian(clip = 'off') +
    geom_abline(slope=0, intercept=min(-log10(inTab$pvalue[inTab$padj < 0.05 & !is.na(inTab$padj)])), linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_vline(xintercept=0.5, linetype="dashed") +
    annotate("text",  x=c(-Inf,Inf,0), y = Inf, label = myLabs, hjust=c(-1,1,0), vjust=1) +
    ggtitle(myTitle)
}
makeVolcano2 <- function(inTab, myTitle="No Title"){
  inTab$Sig <- "Not Sig"
  inTab$Sig[inTab$log2FoldChange > 0.5 & inTab$padj < 0.05] <- "Up in SD"
  inTab$Sig[inTab$log2FoldChange < -0.5 & inTab$padj < 0.05] <- "Down in SD"
  ggplot(inTab, aes(x=log2FoldChange, y=-log10(padj), color=Sig)) +
    geom_point(size=1, alpha=0.7) +
    coord_cartesian(clip = 'off') +
    scale_color_manual(values=SigColors[unique(inTab$Sig)]) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    geom_abline(slope=0, intercept=-log10(0.05), linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_vline(xintercept=0.5, linetype="dashed") +
    ggtitle(myTitle)
}
makeVolcano(myThing$Global_SingleDonor_All$Results, "Global DEGs; MD vs SD NSCs\nAll Donors")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_Global_SD_vs_MD.svg", height=4, width=6)
makeVolcano(myThing$CFuPN_SingleDonor_All$Results, "CFuPN DEGs MD vs SD NSCs\nAll Donors")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_AllDonors_SD_vs_MD.svg", height=4, width=6)
makeVolcano(myThing$CPN_SingleDonor_All$Results, "CPN DEGs MD vs SD NSCs\nAll Donors")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CPN_AllDonors_SD_vs_MD.svg", height=4, width=6)
makeVolcano(myThing$IP_SingleDonor_All$Results, "IP DEGs MD vs SD NSCs\nAll Donors")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_IP_AllDonors_SD_vs_MD.svg", height=4, width=6)


## Let's make a specific effort with the single-donor vs multi-donor samples to include Batch as a covariate...
## ...since we actually have matched batches there.
SingleMulti <- list()
myProtocol <- "SingleDonor"
for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  for(myDonor in c("CW","H1","Mito210","PGP1","11a")){
    print(paste(CellType,myProtocol,myDonor, sep="_"))
    if(sum(metabulk$Protocol==myProtocol & metabulk$Genotype==myDonor & metabulk$Treatment=="Control" & metabulk$CellType==CellType) > 1){
      SingleMulti[[paste(CellType,myProtocol,myDonor, sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                                                  CellTypes=c(CellType),
                                                                                  Genotypes=c(myDonor),
                                                                                  Protocols=c(myProtocol,"MultiDonor"),
                                                                                  myEffect="Protocol",
                                                                                  addCovs = "Batch")
    }else{
      print("Not enough replicates.")
    }
  }
}
for(myDonor in c("CW","H1","Mito210","PGP1","11a")){
  SingleMulti[[paste("Global",myProtocol,myDonor, sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                                            IsGlobal = TRUE,
                                                                            Genotypes=c(myDonor),
                                                                            Protocols=c(myProtocol,"MultiDonor"),
                                                                            myEffect="Protocol",
                                                                            addCovs = "Batch")
}
for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  SingleMulti[[paste(CellType,myProtocol,"Global", sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                                              CellTypes=c(CellType),
                                                                              Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                                              Protocols=c(myProtocol,"MultiDonor"),
                                                                              myEffect="Protocol",
                                                                              addCovs = "Genotype + Batch")
}
SingleMulti$Global_SingleDonor_Global <- GetDE_Results(Treatments=c("Control"),
                                                       IsGlobal=TRUE,
                                                       Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                       Protocols=c(myProtocol,"MultiDonor"),
                                                       myEffect="Protocol",
                                                       addCovs = "Genotype + Batch")

## Global SD vs MD, one donor at a time:
makeVolcano(SingleMulti$Global_SingleDonor_CW$Results, "Global DEGs; MD vs SD NSCs\nCW-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_CW_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$Global_SingleDonor_H1$Results, "Global DEGs; MD vs SD NSCs\nH1-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_H1_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$Global_SingleDonor_Mito210$Results, "Global DEGs; MD vs SD NSCs\nMito210-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_Mito210_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$Global_SingleDonor_PGP1$Results, "Global DEGs; MD vs SD NSCs\nPGP1-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_PGP1_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$Global_SingleDonor_11a$Results, "Global DEGs; MD vs SD NSCs\n11a-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_Global_11a_SD_vs_MD_BatchCovariate.svg", height=4, width=6)

## CFuPN SD vs MD, one donor at a time:
makeVolcano(SingleMulti$CFuPN_SingleDonor_CW$Results, "CFuPN DEGs; MD vs SD NSCs\nCW-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_CW_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$CFuPN_SingleDonor_H1$Results, "CFuPN DEGs; MD vs SD NSCs\nH1-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_H1_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$CFuPN_SingleDonor_Mito210$Results, "CFuPN DEGs; MD vs SD NSCs\nMito210-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_Mito210_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$CFuPN_SingleDonor_PGP1$Results, "CFuPN DEGs; MD vs SD NSCs\nPGP1-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_PGP1_SD_vs_MD_BatchCovariate.svg", height=4, width=6)
makeVolcano(SingleMulti$CFuPN_SingleDonor_11a$Results, "CFuPN DEGs; MD vs SD NSCs\n11a-only")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_CFuPN_11a_SD_vs_MD_BatchCovariate.svg", height=4, width=6)

NumSigDiff <- matrix(ncol=6, nrow=7, data=NA)
colnames(NumSigDiff) <- c("Global","CW","H1","Mito210","PGP1","11a")
rownames(NumSigDiff) <- c("Global","aRG","oRG","IP","PN","CFuPN","CPN")
for(j in colnames(NumSigDiff)){
  for(i in rownames(NumSigDiff)){
    NumSigDiff[i,j] <- nrow(getSigGenes(SingleMulti[[paste(i,"SingleDonor",j, sep="_")]]))
  }
}

NumSigDiff2 <- matrix(ncol=2, nrow=6, data=NA)
colnames(NumSigDiff2) <- c("Mito210","PGP1")
rownames(NumSigDiff2) <- c("aRG","oRG","IP","PN","CFuPN","CPN")
for(j in colnames(NumSigDiff2)){
  for(i in rownames(NumSigDiff2)){
    NumSigDiff2[i,j] <- nrow(getSigGenes(myThing[[paste(i,"SingleDonor",j, sep="-")]]))
  }
}

library(openxlsx)

CW_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  CW_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_CW")]])
  tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_CW")]]$SimpUp@result
  tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_CW")]]$SimpDown@result
  tmp1$Direction <- rep("Up", nrow(tmp1))
  tmp2$Direction <- rep("Down", nrow(tmp2))
  CW_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(CW_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_CW.xlsx", rowNames=TRUE)
H1_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  H1_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_H1")]])
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_H1")]]$UpGo)){
    tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_H1")]]$SimpUp@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_H1")]]$DownGo)){
    tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_H1")]]$SimpDown@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  H1_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(H1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_H1.xlsx", rowNames=TRUE)
PGP1_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  PGP1_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_PGP1")]])
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_PGP1")]]$UpGo)){
    tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_PGP1")]]$SimpUp@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_PGP1")]]$DownGo)){
    tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_PGP1")]]$SimpDown@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  PGP1_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(PGP1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_PGP1.xlsx", rowNames=TRUE)
Mito210_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  Mito210_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_Mito210")]])
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_Mito210")]]$UpGo)){
    tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_Mito210")]]$SimpUp@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_Mito210")]]$DownGo)){
    tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_Mito210")]]$SimpDown@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  Mito210_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(Mito210_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_Mito210.xlsx", rowNames=TRUE)
a11a_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  a11a_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_11a")]])
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_11a")]]$UpGo)){
    tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_11a")]]$SimpUp@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_11a")]]$DownGo)){
    tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_11a")]]$SimpDown@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  a11a_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(a11_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_11a.xlsx", rowNames=TRUE)
All_SigGenes <- list()
for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
  All_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_Global")]])
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_Global")]]$UpGo)){
    tmp1 <- SingleMulti[[paste0(i,"_SingleDonor_Global")]]$SimpUp@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(SingleMulti[[paste0(i,"_SingleDonor_Global")]]$DownGo)){
    tmp2 <- SingleMulti[[paste0(i,"_SingleDonor_Global")]]$SimpDown@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  All_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(All_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_and_GO_AllDonors.xlsx", rowNames=TRUE)



GetDE_Results2 <- function(Donor1, Donor2){
  mySamps <- ((gMeta$Treatment == "Control" & gMeta$Protocol=="MultiDonor" & gMeta$Genotype==Donor1) |
    (gMeta$Treatment == "Control" & gMeta$Protocol=="SingleDonor" & gMeta$Genotype==Donor2))
  myPB <- list(PB = gPB[,mySamps], META = gMeta[mySamps,])
  print(paste(Donor1, Donor2))
  print(table(myPB$META$Protocol, myPB$META$Genotype))
  myRes <- runDESeq(myPB$PB,
                    myPB$META,
                    myEffect="Protocol",
                    addCovs="Batch")
  SigRes <- myRes[!is.na(myRes$padj) & myRes$padj < 0.1,]
  UpGenes <- rownames(SigRes)[SigRes$log2FoldChange > 0.5]
  DownGenes <- rownames(SigRes)[SigRes$log2FoldChange < -0.5]
  upGO <- enrichGO(gene=UpGenes,
                   OrgDb = org.Hs.eg.db,
                   keyType="SYMBOL",
                   ont="BP",
                   universe=rownames(myPB$PB))
  downGO <- enrichGO(gene=DownGenes,
                     OrgDb = org.Hs.eg.db,
                     keyType="SYMBOL",
                     ont="BP",
                     universe=rownames(myPB$PB))
  tryCatch( { upSimp <- simplify(upGO, cutoff=0.4) },
            warning = function(w) {upSimp <<- "None"},
            error = function(e) {upSimp <<- "None"})
  tryCatch( { downSimp <- simplify(downGO, cutoff=0.4) },
            warning = function(w) {downSimp <<- "None"},
            error = function(e) {downSimp <<- "None"})
  return(list(Results=as.data.frame(myRes), UpGo=upGO, DownGo=downGO, SimpUp=upSimp, SimpDown=downSimp))
}

ComparingDonorsBetweenProtocols <- list()
for(i in myDonors){
  for(j in myDonors){
    ComparingDonorsBetweenProtocols[[paste(i,j,"_")]] <- GetDE_Results2(i, j)
  }
}
mySigMat <- matrix(ncol=5, nrow=5)
colnames(mySigMat) <- rownames(mySigMat) <- myDonors
for(i in myDonors){
  for(j in myDonors){
    mySigMat[i,j] <- nrow(getSigGenes(ComparingDonorsBetweenProtocols[[paste(i,j,"_")]]))
  }
}
meanLFMat <- matrix(ncol=5, nrow=5)
colnames(meanLFMat) <- rownames(meanLFMat) <- myDonors
for(i in myDonors){
  for(j in myDonors){
    tmp1 <- abs(ComparingDonorsBetweenProtocols[[paste(i,j,"_")]]$Results$log2FoldChange)
    tmp1 <- tmp1[!is.na(tmp1)]
    meanLFMat[i,j] <- mean(tmp1)
  }
}
myDF <- data.frame(MD="None", SD="None", NumSig=0, meanFC=0)
for(i in myDonors){
  for(j in myDonors){
    tmp1 <- abs(ComparingDonorsBetweenProtocols[[paste(i,j,"_")]]$Results$log2FoldChange)
    tmp1 <- tmp1[!is.na(tmp1)]
    myDF <- rbind(myDF, c(i, j, nrow(getSigGenes(ComparingDonorsBetweenProtocols[[paste(i,j,"_")]])), mean(tmp1)))
  }
}
myDF <- myDF[-1,]
myDF$NumSig <- as.numeric(myDF$NumSig)
myDF$meanFC <- as.numeric(myDF$meanFC)
myDF$MD <- factor(myDF$MD, levels=myDonors)
myDF$SD <- factor(myDF$SD, levels=rev(myDonors))
ggplot(myDF, aes(x=MD, y=SD, color=-meanFC, size=1/NumSig)) +
  geom_point() +
  theme(panel.background = element_blank()) +
  scale_size_continuous(range=c(0,5)) +
  scale_color_gradient(low="blue", high="red")
ggplot(myDF, aes(x=MD, y=SD, color=meanFC)) +
  geom_point(size=7) +
  theme(panel.background = element_blank()) +
  scale_color_gradient(low="blue", high="red")
ggplot(myDF, aes(x=MD, y=SD, color=NumSig)) +
  geom_point(size=7) +
  theme(panel.background = element_blank()) +
  scale_size_continuous(range=c(0,5)) +
  scale_color_gradient(low="blue", high="red")
Heatmap(mySigMat, column_order = 1:5, row_order=1:5)

topAnno <- HeatmapAnnotation(Single_Donor=myDonors,
                             col=list(Single_Donor=c("CW"="#0A9396",
                                              "H1"="#E9D8A6",
                                              "Mito210"="#CA6702",
                                              "PGP1"="#EE9B00",
                                              "11a"="darkred")))
sideAnno <- rowAnnotation(Multi_Donor=myDonors,
                          col=list(Multi_Donor=c("CW"="#0A9396",
                                           "H1"="#E9D8A6",
                                           "Mito210"="#CA6702",
                                           "PGP1"="#EE9B00",
                                           "11a"="darkred")))
p1 <- Heatmap(meanLFMat, column_order = 1:5, row_order=1:5,
        top_annotation = topAnno, left_annotation = sideAnno,
        col=circlize::colorRamp2(c(min(meanLFMat), median(meanLFMat), max(meanLFMat)), c("red", "white", "blue")),
        name="Average\nLog(FoldChange)")
svglite("/Users/tfaits/Documents/Chimera/Plots/SD_MD_Donor-Donor_Comparison_meanLFC_Heatmap.svg", height=4.5, width=7)
grid.newpage()
draw(p1)
dev.off()



NPC_vs_MD <- list()
for(i in c("Global","oRG","IP","PN","CFuPN","CPN")){
  j <- "AllDonors"
  if(i == "Global"){
    NPC_vs_MD[[paste(i,j,sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                IsGlobal=TRUE,
                                                Genotypes=c("H1","Mito210","PGP1"),
                                                Protocols=c("MultiDonor","NPC"),
                                                myEffect = "Protocol",
                                                addCovs = "Genotype")
  }else{
    NPC_vs_MD[[paste(i,j,sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                CellTypes = i,
                                                Genotypes=c("H1","Mito210","PGP1"),
                                                Protocols=c("MultiDonor","NPC"),
                                                myEffect = "Protocol",
                                                addCovs = "Genotype")
  }
  for(j in c("H1","Mito210","PGP1")){
    if(i == "Global"){
      NPC_vs_MD[[paste(i,j,sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                       IsGlobal=TRUE,
                                                       Genotypes=c(j),
                                                       Protocols=c("MultiDonor","NPC"),
                                                       myEffect = "Protocol")
    }else{
      NPC_vs_MD[[paste(i,j,sep="_")]] <- GetDE_Results(Treatments=c("Control"),
                                                      CellTypes = i,
                                                      Genotypes=c(j),
                                                      Protocols=c("MultiDonor","NPC"),
                                                      myEffect = "Protocol")
    }
  }
}

H1_SigGenes <- list()
for(i in c("Global","oRG","IP","PN","CFuPN","CPN")){
  H1_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(NPC_vs_MD[[paste0(i,"_H1")]])
  if(!is.null(NPC_vs_MD[[paste0(i,"_H1")]]$UpGo)){
    tmp1 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_H1")]]$UpGo, cutoff=0.4)@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(NPC_vs_MD[[paste0(i,"_H1")]]$DownGo)){
    tmp2 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_H1")]]$DownGo, cutoff=0.4)@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  H1_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(H1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/NPC_vs_NSC_DEGs_and_GO_H1.xlsx", rowNames=TRUE)
Mito210_SigGenes <- list()
for(i in c("Global","oRG","IP","PN","CFuPN","CPN")){
  Mito210_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(NPC_vs_MD[[paste0(i,"_Mito210")]])
  if(!is.null(NPC_vs_MD[[paste0(i,"_Mito210")]]$UpGo)){
    tmp1 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_Mito210")]]$UpGo, cutoff=0.4)@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(NPC_vs_MD[[paste0(i,"_Mito210")]]$DownGo)){
    tmp2 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_Mito210")]]$DownGo, cutoff=0.4)@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  Mito210_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(Mito210_SigGenes, "/Users/tfaits/Documents/Chimera/Results/NPC_vs_NSC_DEGs_and_GO_Mito210.xlsx", rowNames=TRUE)
PGP1_SigGenes <- list()
for(i in c("Global","oRG","IP","PN","CFuPN","CPN")){
  PGP1_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(NPC_vs_MD[[paste0(i,"_PGP1")]])
  if(!is.null(NPC_vs_MD[[paste0(i,"_PGP1")]]$UpGo)){
    tmp1 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_PGP1")]]$UpGo, cutoff=0.4)@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(NPC_vs_MD[[paste0(i,"_PGP1")]]$DownGo)){
    tmp2 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_PGP1")]]$DownGo, cutoff=0.4)@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  PGP1_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(PGP1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/NPC_vs_NSC_DEGs_and_GO_PGP1.xlsx", rowNames=TRUE)
All_SigGenes <- list()
for(i in c("Global","oRG","IP","PN","CFuPN","CPN")){
  All_SigGenes[[paste0(i,"_genes")]] <- getSigGenes(NPC_vs_MD[[paste0(i,"_AllDonors")]])
  if(!is.null(NPC_vs_MD[[paste0(i,"_AllDonors")]]$UpGo)){
    tmp1 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_AllDonors")]]$UpGo, cutoff=0.4)@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
  }else{
    tmp1 <- matrix(ncol=10)
    colnames(tmp1) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  if(!is.null(NPC_vs_MD[[paste0(i,"_AllDonors")]]$DownGo)){
    tmp2 <- clusterProfiler::simplify(NPC_vs_MD[[paste0(i,"_AllDonors")]]$DownGo, cutoff=0.4)@result
    tmp2$Direction <- rep("Down", nrow(tmp2))
  }else{
    tmp2 <- matrix(ncol=10)
    colnames(tmp2) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Direction")
  }
  All_SigGenes[[paste0(i,"_GO")]] <- rbind(tmp1,tmp2)
}
write.xlsx(All_SigGenes, "/Users/tfaits/Documents/Chimera/Results/NPC_vs_NSC_DEGs_and_GO_AllDonors.xlsx", rowNames=TRUE)



p1 <- makeVolcano(NPC_vs_MD$Global_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - Global Differences") + NoLegend()
p2 <- makeVolcano(NPC_vs_MD$oRG_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - oRGs") + NoLegend()
p3 <- makeVolcano(NPC_vs_MD$IP_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - IPs") + NoLegend()
p4 <- makeVolcano(NPC_vs_MD$PN_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - PNs") + NoLegend()
p5 <- makeVolcano(NPC_vs_MD$CFuPN_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - CFuPNs") + NoLegend()
p6 <- makeVolcano(NPC_vs_MD$CPN_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - CPNs") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_NPCvsNSC_AllDonors_EachCellType_Feb5.svg", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_NPCvsNSC_AllDonors_EachCellType_Feb5.jpeg", res=500, units="in", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()

p1 <- makeVolcano(NPC_vs_MD$Global_AllDonors$Results, "NSC vs NPC Chimeroids\nAll Donors - Global Differences") + NoLegend()
p2 <- makeVolcano(NPC_vs_MD$Global_H1$Results, "NSC vs NPC Chimeroids\nH1 - Global Differences") + NoLegend()
p3 <- makeVolcano(NPC_vs_MD$Global_Mito210$Results, "NSC vs NPC Chimeroids\nMito210 - Global Differences") + NoLegend()
p4 <- makeVolcano(NPC_vs_MD$Global_PGP1$Results, "NSC vs NPC Chimeroids\nPGP1 - Global Differences") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_NPCvsNSC_EachDonor_Global_Feb5.svg", height=8, width=8)
grid.newpage()
grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_NPCvsNSC_EachDonor_Global_Feb5.jpeg", res=500, units="in", height=8, width=8)
grid.newpage()
grid.arrange(p1,p2,p3,p4, ncol=2)
dev.off()




# H1_SigGenes <- list()
# for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
#   H1_SigGenes[[i]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_H1")]])
# }
# write.xlsx(H1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_H1_Global.xlsx", rowNames=TRUE)
# PGP1_SigGenes <- list()
# for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
#   PGP1_SigGenes[[i]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_PGP1")]])
# }
# write.xlsx(PGP1_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_PGP1_Global.xlsx", rowNames=TRUE)
# Mito_SigGenes <- list()
# for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
#   Mito_SigGenes[[i]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_Mito210")]])
# }
# write.xlsx(Mito_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_Mito210_Global.xlsx", rowNames=TRUE)
# a11_SigGenes <- list()
# for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
#   a11_SigGenes[[i]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_11a")]])
# }
# write.xlsx(a11_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_11a_Global.xlsx", rowNames=TRUE)
# Global_SigGenes <- list()
# for(i in c("Global","aRG","oRG","IP","PN","CFuPN","CPN")){
#   Global_SigGenes[[i]] <- getSigGenes(SingleMulti[[paste0(i,"_SingleDonor_Global")]])
# }
# write.xlsx(Global_SigGenes, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_AllDonors_Global.xlsx", rowNames=TRUE)

saveRDS(SingleMulti, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_MultiDonor_DEGs_FullObject.RDS")


myMat <- matrix(ncol=length(myThing), nrow=length(myThing), data=NA)
colnames(myMat) <- names(myThing)
rownames(myMat) <- names(myThing)
for(j in colnames(myMat)){
  for(i in rownames(myMat)){
    myMat[i,j] <- Correlate_Results(myThing[[i]],myThing[[j]])
  }
}
tmp <- names(myThing)
myMet <- data.frame(CellType = sapply(tmp, function(x){strsplit(x,"-")[[1]][1]}),
                    Protocol = sapply(tmp, function(x){strsplit(x,"-")[[1]][2]}),
                    Donor = sapply(tmp, function(x){strsplit(x,"-")[[1]][3]}))
topAnnotation <- HeatmapAnnotation(CellType=myMet$CellType,
                                   Protocol=myMet$Protocol,
                                   Donor=myMet$Donor,
                                   col=list(CellType=TypeCols[unique(myMet$CellType)],
                                         Protocol=ProtCols[unique(myMet$Protocol)],
                                         Donor=DonorCols[unique(myMet$Donor)]))
svglite("Plots/CorrHeatmap_BetweenProtocols.svg",height=8,width=8)
grid.newpage()
p0 <- Heatmap(myMat, top_annotation = topAnnotation, col=col_fun,
        show_row_names = FALSE, show_column_names = FALSE,
        column_title = "Correlation of LogFoldChange Gene Expression\nBetween NSC Multidonor Chimeroids and Other Protocols\nIn Mito210 and PGP1 Lines, in Major Cell Types")
draw(p1)
dev.off()

tA1 <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="Velasco"],
                         Donor=myMet$Donor[myMet$Protocol=="Velasco"],
                         col=list(CellType=TypeCols[unique(myMet$CellType)],
                                  Donor=DonorCols[unique(myMet$Donor)]))
p1 <- Heatmap(myMat[myMet$Protocol=="Velasco",myMet$Protocol=="Velasco"], top_annotation = tA1, col=col_fun,
        show_row_names = FALSE, show_column_names = FALSE, row_order=1:12, column_order=1:12,
        column_title = "Correlation of LogFoldChange Gene Expression\nBetween NSC Multidonor Chimeroids and Velasco\nIn Mito210 and PGP1 Lines, in Major Cell Types")

tA1b <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="Velasco" & myMet$Donor=="PGP1"],
                         Donor=myMet$Donor[myMet$Protocol=="Velasco" & myMet$Donor=="PGP1"],
                         col=list(CellType=TypeCols[unique(myMet$CellType)],
                                  Donor=DonorCols[unique(myMet$Donor)]))
sA1b <- rowAnnotation(CellType=myMet$CellType[myMet$Protocol=="Velasco" & myMet$Donor=="Mito210"],
                      Donor=myMet$Donor[myMet$Protocol=="Velasco" & myMet$Donor=="Mito210"],
                      col=list(CellType=TypeCols[unique(myMet$CellType)],
                               Donor=DonorCols[unique(myMet$Donor)]))
p1b <- Heatmap(myMat[myMet$Protocol=="Velasco" & myMet$Donor=="Mito210",myMet$Protocol=="Velasco" & myMet$Donor=="PGP1"], top_annotation = tA1b, left_annotation = sA1b, col=col_fun,
              show_row_names = FALSE, show_column_names = FALSE,
              column_title = "Correlation Between PGP1 and Mito210\nlogFC Between Velasco and Multi-Donor NSC Chimeroids")
svglite("Plots/CorrHeatmap_Velasco_vs_MultiDonor.svg",height=5.5,width=6)
grid.newpage()
draw(p1b)
dev.off()

tA2b <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="NPC" & myMet$Donor=="PGP1"],
                          Donor=myMet$Donor[myMet$Protocol=="NPC" & myMet$Donor=="PGP1"],
                          col=list(CellType=TypeCols[unique(myMet$CellType)],
                                   Donor=DonorCols[unique(myMet$Donor)]))
sA2b <- rowAnnotation(CellType=myMet$CellType[myMet$Protocol=="NPC" & myMet$Donor=="Mito210"],
                      Donor=myMet$Donor[myMet$Protocol=="NPC" & myMet$Donor=="Mito210"],
                      col=list(CellType=TypeCols[unique(myMet$CellType)],
                               Donor=DonorCols[unique(myMet$Donor)]))
p2b <- Heatmap(myMat[myMet$Protocol=="NPC" & myMet$Donor=="Mito210",myMet$Protocol=="NPC" & myMet$Donor=="PGP1"], top_annotation = tA2b, left_annotation = sA2b, col=col_fun,
               show_row_names = FALSE, show_column_names = FALSE,
               column_title = "Correlation Between PGP1 and Mito210\nlogFC Between NPC and Multi-Donor NSC Chimeroids")
svglite("Plots/CorrHeatmap_NPC_vs_MultiDonor.svg",height=5,width=5.5)
grid.newpage()
draw(p2b)
dev.off()

tA3b <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="SingleDonor" & myMet$Donor=="PGP1"],
                          Donor=myMet$Donor[myMet$Protocol=="SingleDonor" & myMet$Donor=="PGP1"],
                          col=list(CellType=TypeCols[unique(myMet$CellType)],
                                   Donor=DonorCols[unique(myMet$Donor)]))
sA3b <- rowAnnotation(CellType=myMet$CellType[myMet$Protocol=="SingleDonor" & myMet$Donor=="Mito210"],
                      Donor=myMet$Donor[myMet$Protocol=="SingleDonor" & myMet$Donor=="Mito210"],
                      col=list(CellType=TypeCols[unique(myMet$CellType)],
                               Donor=DonorCols[unique(myMet$Donor)]))
p3b <- Heatmap(myMat[myMet$Protocol=="SingleDonor" & myMet$Donor=="Mito210",myMet$Protocol=="SingleDonor" & myMet$Donor=="PGP1"], top_annotation = tA3b, left_annotation = sA3b, col=col_fun,
               show_row_names = FALSE, show_column_names = FALSE,
               column_title = "Correlation Between PGP1 and Mito210\nlogFC Between SingleDonor and Multi-Donor NSC Chimeroids")
svglite("Plots/CorrHeatmap_SingleDonor_vs_MultiDonor.svg",height=5.5,width=6)
grid.newpage()
draw(p3b)
dev.off()

sG1 <- rownames(myThing$`CPN-Velasco-Mito210`$Results)
sG1 <- sG1[sG1 %in% rownames(myThing$`CPN-Velasco-PGP1`$Results) &
             sG1 %in% rownames(myThing$`CPN-SingleDonor-PGP1`$Results) &
             sG1 %in% rownames(myThing$`CPN-SingleDonor-Mito210`$Results) &
             sG1 %in% rownames(myThing$`CPN-NPC-PGP1`$Results) &
             sG1 %in% rownames(myThing$`CPN-NPC-Mito210`$Results)]

tempLFC <- data.frame(Gene=sG1,
                      Velasco_Mito=myThing$`CPN-Velasco-Mito210`$Results[sG1,"log2FoldChange"],
                      Velasco_PGP1=myThing$`CPN-Velasco-PGP1`$Results[sG1,"log2FoldChange"],
                      NPC_Mito=myThing$`CPN-NPC-Mito210`$Results[sG1,"log2FoldChange"],
                      NPC_PGP1=myThing$`CPN-NPC-PGP1`$Results[sG1,"log2FoldChange"],
                      SingleDonor_Mito=myThing$`CPN-SingleDonor-Mito210`$Results[sG1,"log2FoldChange"],
                      SingleDonor_PGP1=myThing$`CPN-SingleDonor-PGP1`$Results[sG1,"log2FoldChange"],
                      Velasco_Mito_Sig=myThing$`CPN-Velasco-Mito210`$Results[sG1,"padj"],
                      Velasco_PGP1_Sig=myThing$`CPN-Velasco-PGP1`$Results[sG1,"padj"],
                      NPC_Mito_Sig=myThing$`CPN-NPC-Mito210`$Results[sG1,"padj"],
                      NPC_PGP1_Sig=myThing$`CPN-NPC-PGP1`$Results[sG1,"padj"],
                      SingleDonor_Mito_Sig=myThing$`CPN-SingleDonor-Mito210`$Results[sG1,"padj"],
                      SingleDonor_PGP1_Sig=myThing$`CPN-SingleDonor-PGP1`$Results[sG1,"padj"])


p1 <- ggplot(tempLFC, aes(x=Velasco_Mito, y=Velasco_PGP1)) +
  geom_point(color="#196099") +
  ggtitle(paste0("MultiDonor NSC Chimeroid CPNs\nvs Velasco CPNs (r^2=",round(cor(tempLFC$Velasco_Mito, tempLFC$Velasco_PGP1),digits=3),")")) +
  xlab("logFC in Mito210") + ylab("logFC in PGP1") +
  xlim(c(-16,11)) + ylim(c(-16,11)) +
  theme_classic() +
  theme(plot.title=element_text(face="bold",hjust=0.5)) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", linewidth=0.5)
p1
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Velasco_vs_MD_Mito_vs_PGP_Scatterplot.svg", height=4, width=4)

p2 <- ggplot(tempLFC, aes(x=NPC_Mito, y=NPC_PGP1)) +
  geom_point(color="#cd8432") +
  ggtitle(paste0("MultiDonor NSC Chimeroid CPNs\nvs NPC Chimeroid CPNs (r^2=",round(cor(tempLFC$NPC_Mito, tempLFC$NPC_PGP1),digits=3),")")) +
  xlab("logFC in Mito210") + ylab("logFC in PGP1") +
  xlim(c(-16,11)) + ylim(c(-16,11)) +
  theme_classic() +
  theme(plot.title=element_text(face="bold",hjust=0.5)) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", linewidth=0.5)
p2
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/NPC_vs_MD_Mito_vs_PGP_Scatterplot.svg", height=4, width=4)

p3 <- ggplot(tempLFC, aes(x=SingleDonor_Mito, y=SingleDonor_PGP1)) +
  geom_point(color="#8cdb86") +
  ggtitle(paste0("MultiDonor NSC Chimeroid CPNs\nvs SingleDonor NSC CPNs (r^2=",round(cor(tempLFC$SingleDonor_Mito, tempLFC$SingleDonor_PGP1),digits=3),")")) +
  xlab("logFC in Mito210") + ylab("logFC in PGP1") +
  xlim(c(-16,11)) + ylim(c(-16,11)) +
  theme_classic() +
  theme(plot.title=element_text(face="bold",hjust=0.5)) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", linewidth=0.5)
p3
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Velasco_vs_MD_Mito_vs_PGP_Scatterplot.svg", height=4, width=4)

cowplot::plot_grid(p1,p2,p3, ncol=3)
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/EachProtocol_vs_MD_Mito_vs_PGP_Scatterplot.svg", height=4, width=12)

tA2 <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="NPC"],
                         Donor=myMet$Donor[myMet$Protocol=="NPC"],
                         col=list(CellType=TypeCols[unique(myMet$CellType)],
                                  Donor=DonorCols[unique(myMet$Donor)]))
p2 <- Heatmap(myMat[myMet$Protocol=="NPC",myMet$Protocol=="NPC"], top_annotation = tA1, col=col_fun,
              show_row_names = FALSE, show_column_names = FALSE,
              column_title = "Correlation of LogFoldChange Gene Expression\nBetween NSC Multidonor Chimeroids and NPC\nIn Mito210 and PGP1 Lines, in Major Cell Types")
tA3 <- HeatmapAnnotation(CellType=myMet$CellType[myMet$Protocol=="SingleDonor"],
                         Donor=myMet$Donor[myMet$Protocol=="SingleDonor"],
                         col=list(CellType=TypeCols[unique(myMet$CellType)],
                                  Donor=DonorCols[unique(myMet$Donor)]))
p3 <- Heatmap(myMat[myMet$Protocol=="SingleDonor",myMet$Protocol=="SingleDonor"], top_annotation = tA1, col=col_fun,
              show_row_names = FALSE, show_column_names = FALSE,
              column_title = "Correlation of LogFoldChange Gene Expression\nBetween NSC Multidonor Chimeroids and SingleDonor\nIn Mito210 and PGP1 Lines, in Major Cell Types")

p1+p2+p3


## Interesting... The DEGs between MD and SD are a lot less consistent between cell types and between PGP1 and Mito210...
## ...compared to the DEGs between MD and Velasco or MD and NPCs.
## Could this be because the MDs and SDs are actually so similar that we're not really finding DEGs?
## ...Check!
ggplot(myThing$`CPN-SingleDonor-PGP1`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()
ggplot(myThing$`CPN-Velasco-PGP1`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()
ggplot(myThing$`CPN-NPC-PGP1`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()

ggplot(myThing$`CPN-SingleDonor-Mito210`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()
ggplot(myThing$`CPN-Velasco-Mito210`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()
ggplot(myThing$`CPN-NPC-Mito210`$Results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point()





##
myThing2 <- list()
for(CellType in c("oRG","IP","PN","CPN","CFuPN")){
  for(myProtocol in c("Velasco","NPC","SingleDonor","MultiDonor")){
      print(paste(CellType,myProtocol, sep="-"))
      myThing2[[paste(CellType,myProtocol,myDonor, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                              CellTypes=c(CellType),
                                                                              Genotypes=c("PGP1","Mito210"),
                                                                              Protocols=c(myProtocol),
                                                                              myEffect="Genotype")
  }
}
CellType <- "aRG"
for(myProtocol in c("Velasco","SingleDonor","MultiDonor")){
  print(paste(CellType,myProtocol, sep="-"))
  myThing2[[paste(CellType,myProtocol,myDonor, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                           CellTypes=c(CellType),
                                                                           Genotypes=c("PGP1","Mito210"),
                                                                           Protocols=c(myProtocol),
                                                                           myEffect="Genotype")
}

myMat2 <- matrix(ncol=length(myThing2), nrow=length(myThing2), data=NA)
colnames(myMat2) <- names(myThing2)
rownames(myMat2) <- names(myThing2)
for(j in colnames(myMat2)){
  for(i in rownames(myMat2)){
    myMat2[i,j] <- Correlate_Results(myThing2[[i]],myThing2[[j]])
  }
}
tmp <- names(myThing2)
myMet2 <- data.frame(CellType = sapply(tmp, function(x){strsplit(x,"-")[[1]][1]}),
                    Protocol = sapply(tmp, function(x){strsplit(x,"-")[[1]][2]}))
colnames(myMat2) <- sapply(colnames(myMat2), function(x){paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
rownames(myMat2) <- sapply(rownames(myMat2), function(x){paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
topAnnotation2 <- HeatmapAnnotation(CellType=myMet2$CellType,
                                   Protocol=myMet2$Protocol,
                                   col=list(CellType=TypeCols[unique(myMet2$CellType)],
                                            Protocol=ProtCols[unique(myMet2$Protocol)]))
svglite("Plots/CorrHeatmap_BetweenDonors_PGP1Mito210.svg",height=5.5,width=6)
grid.newpage()
p1 <- Heatmap(myMat2, top_annotation = topAnnotation2, col=col_fun,
        show_row_names = FALSE, show_column_names = FALSE,
        column_title = "Correlation of LogFoldChange Gene Expression\nBetween Mito210 and PGP1 Lines\nAcross Protocols in Major Cell Types")
draw(p1)
dev.off()

myTmp <- myMat2[myMet2$CellType=="aRG",myMet2$CellType=="aRG"]
myTmp <- rbind(myTmp[1,],c(0,0,0),myTmp[2:3,])
myTmp <- cbind(myTmp[,1],c(0,0,0,0), myTmp[,2:3])
colnames(myTmp)[1:2] <- c("aRG-Velasco","aRG-NPC")
rownames(myTmp)[1:2] <- c("aRG-Velasco","aRG-NPC")
myTmp <- cbind(myTmp,
               myMat2[myMet2$CellType=="oRG",myMet2$CellType=="oRG"],
               myMat2[myMet2$CellType=="IP",myMet2$CellType=="IP"],
               myMat2[myMet2$CellType=="PN",myMet2$CellType=="PN"],
               myMat2[myMet2$CellType=="CPN",myMet2$CellType=="CPN"],
               myMat2[myMet2$CellType=="CFuPN",myMet2$CellType=="CFuPN"])
topAnnotation2 <- HeatmapAnnotation(Protocol=sapply(colnames(myTmp), function(x){strsplit(x,"-")[[1]][2]}),
                                    CellType=sapply(colnames(myTmp), function(x){strsplit(x,"-")[[1]][1]}),
                                    col=list(Protocol=ProtCols[unique(myMet2$Protocol)],
                                             CellType=TypeCols[unique(myMet2$CellType)]))
p1 <- Heatmap(myTmp, top_annotation = topAnnotation2, col=col_fun,
              show_row_names = FALSE, show_column_names = FALSE,
              column_order = c(1,5,9,13,17,21,2,6,10,14,18,22,3,7,11,15,19,23,4,8,12,16,20,24),
              column_title = "Correlation of LogFoldChange Gene Expression\nBetween Mito210 and PGP1 Lines\nAcross Protocols in Major Cell Types")




myThing3 <- list()
for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  for(myProtocol in c("Velasco","NPC","SingleDonor")){
    print(paste(CellType,myProtocol, sep="-"))
    try(myThing3[[paste(CellType,myProtocol, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                         CellTypes=c(CellType),
                                                                         Genotypes=c("PGP1","Mito210"),
                                                                         Protocols=c(myProtocol,"MultiDonor"),
                                                                         myEffect="Protocol",
                                                                         addCovs = "Genotype"))
  }
}
dim(getSigGenes(myThing3$`aRG-SingleDonor`))

Velasco_Vs_mdNSC <- list()
sdNSC_Vs_mdNSC <- list()
NPC_Vs_mdNSC <- list()
UpGeneCounts <- matrix(ncol=6, nrow=3, data=0)
DownGeneCounts <- matrix(ncol=6, nrow=3, data=0)
colnames(UpGeneCounts) <- colnames(DownGeneCounts) <- c("aRG","oRG","IP","PN","CFuPN","CPN")
rownames(UpGeneCounts) <- rownames(DownGeneCounts) <- c("Velasco","NPC","SingleDonor")
Velasco_GO <- list()
NPC_GO <- list()
SD_GO <- list()
simVelasco_GO <- list()
simNPC_GO <- list()
simSD_GO <- list()
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN")){
  Velasco_Vs_mdNSC[[i]] <- getSigGenes(myThing3[[paste(i, "Velasco", sep="-")]])
  UpGeneCounts[1,i] <- sum(Velasco_Vs_mdNSC[[i]]$log2FoldChange > 0)
  DownGeneCounts[1,i] <- sum(Velasco_Vs_mdNSC[[i]]$log2FoldChange < 0)
  sdNSC_Vs_mdNSC[[i]] <- getSigGenes(myThing3[[paste(i, "SingleDonor", sep="-")]])
  UpGeneCounts[3,i] <- sum(sdNSC_Vs_mdNSC[[i]]$log2FoldChange > 0)
  DownGeneCounts[3,i] <- sum(sdNSC_Vs_mdNSC[[i]]$log2FoldChange < 0)
  if(i != "aRG"){
    NPC_Vs_mdNSC[[i]] <- getSigGenes(myThing3[[paste(i, "NPC", sep="-")]])
    UpGeneCounts[2,i] <- sum(NPC_Vs_mdNSC[[i]]$log2FoldChange > 0)
    DownGeneCounts[2,i] <- sum(NPC_Vs_mdNSC[[i]]$log2FoldChange < 0)
  }
}
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN")){
  Velasco_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "Velasco", sep="-")]]$UpGo@result
  Velasco_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "Velasco", sep="-")]]$DownGo@result
  simVelasco_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "Velasco", sep="-")]]$SimpUp@result
  simVelasco_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "Velasco", sep="-")]]$SimpDown@result
  SD_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "SingleDonor", sep="-")]]$UpGo@result
  SD_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "SingleDonor", sep="-")]]$DownGo@result
  simSD_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "SingleDonor", sep="-")]]$SimpUp@result
  simSD_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "SingleDonor", sep="-")]]$SimpDown@result
  if(i != "aRG"){
    NPC_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "NPC", sep="-")]]$UpGo@result
    NPC_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "NPC", sep="-")]]$DownGo@result
    simNPC_GO[[paste(i,"Up", sep="_")]] <- myThing3[[paste(i, "NPC", sep="-")]]$SimpUp@result
    simNPC_GO[[paste(i,"Down", sep="_")]] <- myThing3[[paste(i, "NPC", sep="-")]]$SimpDown@result
  }
}

write.xlsx(Velasco_Vs_mdNSC, "/Users/tfaits/Documents/Chimera/Results/DEG_Table_Velasco_vs_Multidonor.xlsx",rowNames=TRUE)
write.xlsx(sdNSC_Vs_mdNSC, "/Users/tfaits/Documents/Chimera/Results/DEG_Table_Singledonor_vs_Multidonor.xlsx",rowNames=TRUE)
write.xlsx(NPC_Vs_mdNSC, "/Users/tfaits/Documents/Chimera/Results/DEG_Table_NPC_vs_Multidonor.xlsx",rowNames=TRUE)
write.csv(UpGeneCounts,"/Users/tfaits/Documents/Chimera/Results/DEG_Table_GeneCounts_Up.csv", quote=FALSE, col.names = TRUE, row.names=TRUE)
write.csv(DownGeneCounts,"/Users/tfaits/Documents/Chimera/Results/DEG_Table_GeneCounts_Down.csv", quote=FALSE, col.names = TRUE, row.names=TRUE)

write.xlsx(Velasco_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_Velasco_vs_Multidonor.xlsx", rowNames=TRUE)
write.xlsx(NPC_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_NPC_vs_Multidonor.xlsx", rowNames=TRUE)
write.xlsx(SD_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_Singledonor_vs_Multidonor.xlsx", rowNames=TRUE)
write.xlsx(simVelasco_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_Simplified_Velasco_vs_Multidonor.xlsx", rowNames=TRUE)
write.xlsx(simNPC_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_Simplified_NPC_vs_Multidonor.xlsx", rowNames=TRUE)
write.xlsx(simSD_GO, "/Users/tfaits/Documents/Chimera/Results/GOterms_Simplified_Singledonor_vs_Multidonor.xlsx", rowNames=TRUE)

## Find the number of gene detected in each cell type in each protocol
## Looking only in the Mito210/PGP1 Control cells:
detectedGenes <- matrix(ncol=6, nrow=4)
colnames(detectedGenes) <- c("aRG","oRG","IP","PN","CFuPN","CPN")
rownames(detectedGenes) <- c("Velasco","NPC","SingleDonor","MultiDonor")
for(proto in c("Velasco","NPC","SingleDonor","MultiDonor")){
  for(celltype in c("aRG","oRG","IP","PN","CFuPN","CPN")){
    tmp <- pbulk[,metabulk$Genotype %in% c("Mito210","PGP1") &
                   metabulk$Protocol == proto &
                   metabulk$CellType == celltype &
                   metabulk$Treatment == "Control"]
    tmp <- tmp[apply(tmp,1,sum) > 0,] # To be "detected", needs more than a 10-count
    detectedGenes[proto, celltype] <- nrow(tmp)
  }
}
detectedGenes

aRGVel <- josd

## Something to look at: % nonzero expression of genes in single-cell
## ...and/or just looking for overlap in expressed genes overall


## Are these genes the same across lines?
myThing4 <- list()
for(CellType in c("aRG","oRG","IP","PN","CPN","CFuPN")){
  for(myProtocol in c("Velasco","NPC","SingleDonor")){
    for(myDonor in c("Mito210","PGP1")){
      print(paste(CellType, myProtocol, myDonor, sep="-"))
      try(myThing4[[paste(CellType,myProtocol,myDonor, sep="-")]] <- GetDE_Results(Treatments=c("Control"),
                                                                           CellTypes=c(CellType),
                                                                           Genotypes=c(myDonor),
                                                                           Protocols=c(myProtocol,"MultiDonor"),
                                                                           myEffect="Protocol"))
    }
  }
}

Velasco_Matrix <- matrix(ncol=12, nrow=12,data=NA)
SD_Matrix <- matrix(ncol=12, nrow=12,data=NA)
colnames(Velasco_Matrix) <- rownames(Velasco_Matrix) <- c(paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"Velasco","PGP1",sep="-"),paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"Velasco","Mito210",sep="-"))
colnames(SD_Matrix) <- rownames(SD_Matrix) <- c(paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"SingleDonor","PGP1",sep="-"),paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"SingleDonor","Mito210",sep="-"))
for(i in 1:11){
  for(j in (i+1):12){
    Velasco_Matrix[i,j] <- Fisher_Results(myThing4[[colnames(Velasco_Matrix)[i]]], myThing4[[colnames(Velasco_Matrix)[j]]])
    Velasco_Matrix[j,i] <- Velasco_Matrix[i,j]
    SD_Matrix[i,j] <- Fisher_Results(myThing4[[colnames(SD_Matrix)[i]]], myThing4[[colnames(SD_Matrix)[j]]])
    SD_Matrix[j,i] <- SD_Matrix[i,j]
  }
}



Velasco_Matrix <- matrix(ncol=12, nrow=12,data=NA)
SD_Matrix <- matrix(ncol=12, nrow=12,data=NA)
colnames(Velasco_Matrix) <- rownames(Velasco_Matrix) <- c(paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"Velasco","PGP1",sep="-"),paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"Velasco","Mito210",sep="-"))
colnames(SD_Matrix) <- rownames(SD_Matrix) <- c(paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"SingleDonor","PGP1",sep="-"),paste(c("aRG","oRG","IP","PN","CPN","CFuPN"),"SingleDonor","Mito210",sep="-"))
for(i in 1:11){
  for(j in (i+1):12){
    tmpGenesV1 <- rownames(getSigGenes(myThing3[[paste(strsplit(colnames(Velasco_Matrix)[i],"-")[[1]][1:2],collapse="-")]]))
    tmpGenesV2 <- rownames(getSigGenes(myThing3[[paste(strsplit(colnames(Velasco_Matrix)[j],"-")[[1]][1:2],collapse="-")]]))
    tmpGenesS1 <- rownames(getSigGenes(myThing3[[paste(strsplit(colnames(SD_Matrix)[i],"-")[[1]][1:2],collapse="-")]]))
    tmpGenesS2 <- rownames(getSigGenes(myThing3[[paste(strsplit(colnames(SD_Matrix)[j],"-")[[1]][1:2],collapse="-")]]))
    tmpGenesV <- unique(tmpGenesV1, tmpGenesV2)
    tmpGenesS <- unique(tmpGenesS1, tmpGenesS2)
    tmpGenesV <- tmpGenesV[tmpGenesV %in% rownames(myThing4[[colnames(Velasco_Matrix)[i]]]$Results) &
                             tmpGenesV %in% rownames(myThing4[[colnames(Velasco_Matrix)[j]]]$Results)]
    tmpGenesS <- tmpGenesS[tmpGenesS %in% rownames(myThing4[[colnames(SD_Matrix)[i]]]$Results) &
                             tmpGenesS %in% rownames(myThing4[[colnames(SD_Matrix)[j]]]$Results)]
    Velasco_Matrix[i,j] <- cor(myThing4[[colnames(Velasco_Matrix)[i]]]$Results[tmpGenesV,"log2FoldChange"], myThing4[[colnames(Velasco_Matrix)[j]]]$Results[tmpGenesV,"log2FoldChange"])
    Velasco_Matrix[j,i] <- Velasco_Matrix[i,j]
    SD_Matrix[i,j] <- cor(myThing4[[colnames(SD_Matrix)[i]]]$Results[tmpGenesS,"log2FoldChange"], myThing4[[colnames(SD_Matrix)[j]]]$Results[tmpGenesS,"log2FoldChange"])
    SD_Matrix[j,i] <- SD_Matrix[i,j]
  }
}

NPC_Matrix <- matrix(ncol=10, nrow=10,data=NA)
colnames(NPC_Matrix) <- rownames(NPC_Matrix) <- c(paste(c("oRG","IP","PN","CPN","CFuPN"),"NPC","PGP1",sep="-"),paste(c("oRG","IP","PN","CPN","CFuPN"),"NPC","Mito210",sep="-"))


p1 <- VlnPlot(atlas, features=c("CADM2"), group.by="FinalName")
p2 <- VlnPlot(seur, features=c("CADM2"), group.by="MergedAnnotation")
p3 <- VlnPlot(atlas, features=c("OPCML"), group.by="FinalName")
p4 <- VlnPlot(seur, features=c("OPCML"), group.by="MergedAnnotation")
p5 <- VlnPlot(atlas, features=c("EIF3L"), group.by="FinalName")
p6 <- VlnPlot(seur, features=c("EIF3L"), group.by="MergedAnnotation")
cowplot::plot_grid(p1,p2,p3,p4,p5,p6, ncol=2)

p1 <- FeaturePlot(atlas, features=c("CADM2"))
p2 <- FeaturePlot(seur, features=c("CADM2"))
p3 <- FeaturePlot(atlas, features=c("OPCML"))
p4 <- FeaturePlot(seur, features=c("OPCML"))
p5 <- FeaturePlot(atlas, features=c("EIF3L"))
p6 <- FeaturePlot(seur, features=c("EIF3L"))
cowplot::plot_grid(p1,p2,p3,p4,p5,p6, ncol=2)

Make_Heatmap(myGenes=rownames(Velasco_Vs_mdNSC$CPN)[1:100],
             Treatments="Control",
             CellTypes = c("CPN"),
             Genotypes=c("PGP1","Mito210"),
             Protocols=c("Velasco","MultiDonor"),
             myTitle = "DEGs Between Velasco and MultiDonor NSCs in CPNs", scaleData = FALSE)


aRG_MvAtlas_PGP1 <- GetDE_Results(Treatments=c("Control"),
                         CellTypes=c("aRG"),
                         Genotypes = c("PGP1","Mito210"),
                         Protocols=c("Velasco","MultiDonor"),
                         myEffect = "Protocol",
                         addCovs = "Genotype")
# aRG_MvNPC_PGP1 <- GetDE_Results(Treatments=c("Control"),
#                                   CellTypes=c("aRG"),
#                                   Genotypes = c("PGP1"),
#                                   Protocols=c("NPC","MultiDonor"),
#                                   myEffect = "Protocol")
aRG_MvSD_PGP1 <- GetDE_Results(Treatments=c("Control"),
                                  CellTypes=c("aRG"),
                                  Genotypes = c("PGP1","Mito210"),
                                  Protocols=c("SingleDonor","MultiDonor"),
                                  myEffect = "Protocol",
                               addCovs = "Genotype")






## oRG
oRG_MvAtlas_PGP1 <- GetDE_Results(Treatments=c("Control"),
                                  CellTypes=c("oRG"),
                                  Genotypes = c("PGP1"),
                                  Protocols=c("Velasco","MultiDonor"),
                                  myEffect = "Protocol")
oRG_MvAtlas_Mito <- GetDE_Results(Treatments=c("Control"),
                                  CellTypes=c("oRG"),
                                  Genotypes = c("Mito210"),
                                  Protocols=c("Velasco","MultiDonor"),
                                  myEffect = "Protocol")
oRG_MvNPC_PGP1 <- GetDE_Results(Treatments=c("Control"),
                                  CellTypes=c("oRG"),
                                  Genotypes = c("PGP1"),
                                  Protocols=c("NPC","MultiDonor"),
                                  myEffect = "Protocol")
oRG_MvSD_PGP1 <- GetDE_Results(Treatments=c("Control"),
                               CellTypes=c("oRG"),
                               Genotypes = c("PGP1"),
                               Protocols=c("SingleDonor","MultiDonor"),
                               myEffect = "Protocol")


## Begin loooking at VPA:
CPN_Mito <- GetDE_Results(Treatments=c("Control","VPA"),
                          CellTypes = c("CPN"),
                          Genotypes="Mito210",
                          Protocols="MultiDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")
CPN_CW <- GetDE_Results(Treatments=c("Control","VPA"),
                        CellTypes = c("CPN"),
                        Genotypes="CW",
                        Protocols="MultiDonor",
                        myEffect="Treatment",
                        addCovs = "Batch")
CPN_PGP1 <- GetDE_Results(Treatments=c("Control","VPA"),
                          CellTypes = c("CPN"),
                          Genotypes="PGP1",
                          Protocols="MultiDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")
CPN_H1 <- GetDE_Results(Treatments=c("Control","VPA"),
                          CellTypes = c("CPN"),
                          Genotypes="H1",
                          Protocols="MultiDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")
CPN_11a <- GetDE_Results(Treatments=c("Control","VPA"),
                          CellTypes = c("CPN"),
                          Genotypes="11a",
                          Protocols="MultiDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")

CPN_MitoSD <- GetDE_Results(Treatments=c("Control","VPA"),
                            CellTypes = c("CPN"),
                          Genotypes="Mito210",
                          Protocols="SingleDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")
CPN_CWSD <- GetDE_Results(Treatments=c("Control","VPA"),
                          CellTypes = c("CPN"),
                          Genotypes="CW",
                          Protocols="SingleDonor",
                          myEffect="Treatment",
                          addCovs = "Batch")
CPN_PGP1SD <- GetDE_Results(Treatments=c("Control","VPA"),
                            CellTypes = c("CPN"),
                            Genotypes="PGP1",
                            Protocols="SingleDonor",
                            myEffect="Treatment",
                            addCovs = "Batch")
CPN_H1SD <- GetDE_Results(Treatments=c("Control","VPA"),
                            CellTypes = c("CPN"),
                            Genotypes="H1",
                            Protocols="SingleDonor",
                            myEffect="Treatment",
                            addCovs = "Batch")
CPN_11aSD <- GetDE_Results(Treatments=c("Control","VPA"),
                            CellTypes = c("CPN"),
                            Genotypes="11a",
                            Protocols="SingleDonor",
                            myEffect="Treatment",
                            addCovs = "Batch")

sharedGenes <- rownames(CPN_CW$Results)
print(length(sharedGenes)) #17737
sharedGenes <- sharedGenes[sharedGenes %in% rownames(CPN_PGP1$Results) &
                             sharedGenes %in% rownames(CPN_PGP1SD$Results) &
                             sharedGenes %in% rownames(CPN_CWSD$Results) &
                             sharedGenes %in% rownames(CPN_MitoSD$Results) &
                             sharedGenes %in% rownames(CPN_H1SD$Results) &
                             sharedGenes %in% rownames(CPN_11aSD$Results) &
                             sharedGenes %in% rownames(CPN_Mito$Results) &
                             sharedGenes %in% rownames(CPN_H1$Results) &
                             sharedGenes %in% rownames(CPN_11a$Results)]
print(length(sharedGenes)) #14949

myFC <- data.frame(CW_FC=CPN_CW$Results[sharedGenes, "log2FoldChange"],
                   CW_padj=CPN_CW$Results[sharedGenes, "padj"],
                   CWsd_FC=CPN_CWSD$Results[sharedGenes, "log2FoldChange"],
                   CWsd_padj=CPN_CWSD$Results[sharedGenes, "padj"],
                   PGP_FC=CPN_PGP1$Results[sharedGenes, "log2FoldChange"],
                   PGP_padj=CPN_PGP1$Results[sharedGenes, "padj"],
                   PGPsd_FC=CPN_PGP1SD$Results[sharedGenes, "log2FoldChange"],
                   PGPsd_padj=CPN_PGP1SD$Results[sharedGenes, "padj"],
                   Mito_FC=CPN_Mito$Results[sharedGenes, "log2FoldChange"],
                   Mito_padj=CPN_Mito$Results[sharedGenes, "padj"],
                   Mitosd_FC=CPN_MitoSD$Results[sharedGenes, "log2FoldChange"],
                   Mitosd_padj=CPN_MitoSD$Results[sharedGenes, "padj"],
                   H1_FC=CPN_H1$Results[sharedGenes, "log2FoldChange"],
                   H1_padj=CPN_H1$Results[sharedGenes, "padj"],
                   H1sd_FC=CPN_H1SD$Results[sharedGenes, "log2FoldChange"],
                   H1sd_padj=CPN_H1SD$Results[sharedGenes, "padj"],
                   a11_FC=CPN_11a$Results[sharedGenes, "log2FoldChange"],
                   a11_padj=CPN_11a$Results[sharedGenes, "padj"],
                   a11sd_FC=CPN_11aSD$Results[sharedGenes, "log2FoldChange"],
                   a11sd_padj=CPN_11aSD$Results[sharedGenes, "padj"])
myFC$CW_sig <- ifelse(myFC$CW_padj < 0.1, 1, 0)
myFC$CWsd_sig <- ifelse(myFC$CWsd_padj < 0.1, 2, 0)
myFC$CombSig <- myFC$CW_sig + myFC$CWsd_sig

ggplot(myFC, aes(x=CW_FC, y=CWsd_FC, size=-log10(CW_padj), color=-log10(CWsd_padj))) +
  geom_point()
ggplot(myFC, aes(size=CW_FC, color=CWsd_FC, x=-log10(CW_padj)*ifelse(CW_FC>0,1,-1), y=-log10(CWsd_padj)*ifelse(CWsd_FC>0,1,-1),
                 pch=factor(CombSig))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0)



ggplot(myFC, aes(size=CW_FC, color=PGP_FC, x=-log10(CW_padj)*ifelse(CW_FC>0,1,-1), y=-log10(PGP_padj)*ifelse(PGP_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0)

ggplot(myFC, aes(color=(PGP_FC-CW_FC)/(mean(c(PGP_FC,CW_FC))), x=-log10(CW_padj)*ifelse(CW_FC>0,1,-1), y=-log10(PGP_padj)*ifelse(PGP_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0)




p1 <- ggplot(myFC, aes(color=(abs(PGP_FC)-abs(CW_FC)), x=-log10(CW_padj)*ifelse(CW_FC>0,1,-1), y=-log10(PGP_padj)*ifelse(PGP_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0, name="Relative difference in FC") +
  ylim(c(-15,15)) + xlim(c(-20,20)) +
  theme_classic() +
  xlab("Significance (and direction) in CW") + ylab("Significance (and direction) in PGP1") +
  ggtitle("PGP1 vs CW (MultiDonors)")

p2 <- ggplot(myFC, aes(color=(abs(CWsd_FC)-abs(CW_FC)), x=-log10(CW_padj)*ifelse(CW_FC>0,1,-1), y=-log10(CWsd_padj)*ifelse(CWsd_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0, name="Relative difference in FC") +
  ylim(c(-15,15)) + xlim(c(-20,20)) +
  theme_classic() +
  xlab("Significance (and direction) in MultiDonor CW") + ylab("Significance (and direction) in SingleDonor CW") +
  ggtitle("CW MultiDonor vs CW SingleDonor")

p3 <- ggplot(myFC, aes(color=(abs(PGPsd_FC)-abs(PGP_FC)), x=-log10(PGP_padj)*ifelse(PGP_FC>0,1,-1), y=-log10(PGPsd_padj)*ifelse(PGPsd_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0, name="Relative difference in FC") +
  ylim(c(-15,15)) + xlim(c(-20,20)) +
  theme_classic() +
  xlab("Significance (and direction) in MultiDonor PGP1") + ylab("Significance (and direction) in SingleDonor PGP1") +
  ggtitle("PGP1 MultiDonor vs PGP1 SingleDonor")

p4 <- ggplot(myFC, aes(color=(abs(PGPsd_FC)-abs(CWsd_FC)), x=-log10(CWsd_padj)*ifelse(CWsd_FC>0,1,-1), y=-log10(PGPsd_padj)*ifelse(PGPsd_FC>0,1,-1))) +
  geom_point() +
  scale_color_gradient2(low="darkblue",high="darkred",mid="white",midpoint=0, name="Relative difference in FC") +
  ylim(c(-15,15)) + xlim(c(-20,20)) +
  theme_classic() +
  xlab("Significance (and direction) in CW") + ylab("Significance (and direction) in PGP1") +
  ggtitle("PGP1 vs CW (SingleDonors)")

cowplot::plot_grid(p1,p2,p4,p3, ncol=2)

## Is the difference in VPA-induced-logFC between PGP1 and CW similar between protocols? ...No.
ggplot(myFC, aes(x=PGP_FC-CW_FC, y=PGPsd_FC-CWsd_FC)) +
  geom_point()
## Is the difference in VPA-induced-sigDEGs between PGP1 and CW similar between protocols? ...No.
ggplot(myFC, aes(x=log10(PGP_padj)-log10(CW_FC), y=log10(PGPsd_padj)-log10(CWsd_padj))) +
  geom_point()

sigUpCW <- CPN_CW$UpGo@gene
sigUpCWsd <- CPN_CWSD$UpGo@gene
sigUpPGP <- CPN_PGP1$UpGo@gene
sigUpPGPsd <- CPN_PGP1SD$UpGo@gene
sigUpMito <- CPN_Mito$UpGo@gene
sigUpMitosd <- CPN_MitoSD$UpGo@gene
tmp <- list(sigUpCW, sigUpCWsd, sigUpPGP, sigUpPGPsd, sigUpMito, sigUpMitosd)
names(tmp) <- c("sigUpCW", "sigUpCWsd", "sigUpPGP", "sigUpPGPsd", "sigUpMito","sigUpMitosd")

sigUpCW <- CPN_CW$UpGo@gene
sigUpCWsd <- CPN_CWSD$UpGo@gene
sigUpPGP <- CPN_PGP1$UpGo@gene
sigUpPGPsd <- CPN_PGP1SD$UpGo@gene
sigUpMito <- CPN_Mito$UpGo@gene
sigUpMitosd <- CPN_MitoSD$UpGo@gene
sigUpH1 <- CPN_H1$UpGo@gene
sigUpH1sd <- CPN_H1SD$UpGo@gene
sigUp11a <- CPN_11a$UpGo@gene
sigUp11asd <- CPN_11aSD$UpGo@gene
tmp <- list(sigUpCW, sigUpCWsd, sigUpPGP, sigUpPGPsd, sigUpMito, sigUpMitosd, sigUpH1, sigUpH1sd, sigUp11a, sigUp11asd)
names(tmp) <- c("sigUpCW", "sigUpCWsd", "sigUpPGP", "sigUpPGPsd", "sigUpMito","sigUpMitosd", "sigUpH1", "sigUpH1sd", "sigUp11a", "sigUp11asd")

sigDownCW <- CPN_CW$DownGo@gene
sigDownCWsd <- CPN_CWSD$DownGo@gene
sigDownPGP <- CPN_PGP1$DownGo@gene
sigDownPGPsd <- CPN_PGP1SD$DownGo@gene
sigDownMito <- CPN_Mito$DownGo@gene
sigDownMitosd <- CPN_MitoSD$DownGo@gene
sigDownH1 <- CPN_H1$DownGo@gene
sigDownH1sd <- CPN_H1SD$DownGo@gene
sigDown11a <- CPN_11a$DownGo@gene
sigDown11asd <- CPN_11aSD$DownGo@gene
tmpD <- list(sigDownCW, sigDownCWsd, sigDownPGP, sigDownPGPsd, sigDownMito, sigDownMitosd, sigDownH1, sigDownH1sd, sigDown11a, sigDown11asd)
names(tmpD) <- c("sigDownCW", "sigDownCWsd", "sigDownPGP", "sigDownPGPsd", "sigDownMito","sigDownMitosd", "sigDownH1", "sigDownH1sd", "sigDown11a", "sigDown11asd")


myOverlap <- matrix(ncol=10,nrow=10)
colnames(myOverlap) <- names(tmp)
rownames(myOverlap) <- names(tmp)
myOverlap3 <- matrix(ncol=10,nrow=10)
colnames(myOverlap3) <- names(tmp)
rownames(myOverlap3) <- names(tmp)
for(i in names(tmp)){
  for(j in names(tmp)){
    myOverlap[i,j] <- sum(tmp[[i]] %in% tmp[[j]])
    myOverlap3[i,j] <- -log10(fisher.test(matrix(ncol=2, data=c(sum(sharedGenes %ni% c(tmp[[i]],tmp[[j]])),
                                                         sum(tmp[[i]]%ni%tmp[[j]]),
                                                         sum(tmp[[j]]%ni%tmp[[i]]),
                                                         sum(tmp[[i]] %in% tmp[[j]]))))$p.value)
  }
}

mO <- data.frame(SDprot = c("Fake"), MDprot = c("Fake"), FisherSig = c(1), Overlap = 1, ProportionalOverlap = 1)
for(i in c("CW","PGP","Mito","H1","11a")){
  for(j in c("CW","PGP","Mito","H1","11a")){
    mO <- rbind(mO, c(i, j, myOverlap3[paste0("sigUp",i,"sd"),paste0("sigUp",j)],
                      myOverlap[paste0("sigUp",i,"sd"),paste0("sigUp",j)],
                      myOverlap[paste0("sigUp",i,"sd"),paste0("sigUp",j)]/min(length(tmp[paste0("sigUp",j)][[1]]),length(tmp[paste0("sigUp",i,"sd")][[1]]))))
  }
}
mO <- mO[-1,]
mO$FisherSig <- as.numeric(mO$FisherSig)
mO$Overlap <- as.numeric(mO$Overlap)
mO$ProportionalOverlap <- as.numeric(mO$ProportionalOverlap)
mO$Within_Donor <- mO$SDprot == mO$MDprot
mO$SDprot <- factor(mO$SDprot, levels=c("CW","H1","Mito","PGP","11a"))
ggplot(mO, aes(x=Within_Donor, y=FisherSig, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mO$SDprot)])

ggplot(mO, aes(x=Within_Donor, y=ProportionalOverlap, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mO$SDprot)])



## The downs:
myOverlap4 <- matrix(ncol=10,nrow=10)
colnames(myOverlap4) <- names(tmpD)
rownames(myOverlap4) <- names(tmpD)
myOverlap5 <- matrix(ncol=10,nrow=10)
colnames(myOverlap5) <- names(tmpD)
rownames(myOverlap5) <- names(tmpD)
for(i in names(tmpD)){
  for(j in names(tmpD)){
    myOverlap4[i,j] <- sum(tmpD[[i]] %in% tmpD[[j]])
    myOverlap5[i,j] <- -log10(fisher.test(matrix(ncol=2, data=c(sum(sharedGenes %ni% c(tmpD[[i]],tmpD[[j]])),
                                                                sum(tmpD[[i]]%ni%tmpD[[j]]),
                                                                sum(tmpD[[j]]%ni%tmpD[[i]]),
                                                                sum(tmpD[[i]] %in% tmpD[[j]]))))$p.value)
  }
}

mOd <- data.frame(SDprot = c("Fake"), MDprot = c("Fake"), FisherSig = c(1), Overlap = 1, ProportionalOverlap = 1)
for(i in c("CW","PGP","Mito","H1","11a")){
  for(j in c("CW","PGP","Mito","H1","11a")){
    mOd <- rbind(mOd, c(i, j, myOverlap5[paste0("sigDown",i,"sd"),paste0("sigDown",j)],
                      myOverlap4[paste0("sigDown",i,"sd"),paste0("sigDown",j)],
                      myOverlap4[paste0("sigDown",i,"sd"),paste0("sigDown",j)]/min(length(tmpD[paste0("sigDown",j)][[1]]),length(tmpD[paste0("sigDown",i,"sd")][[1]]))))
  }
}
mOd <- mOd[-1,]
mOd$FisherSig <- as.numeric(mOd$FisherSig)
mOd$Overlap <- as.numeric(mOd$Overlap)
mOd$ProportionalOverlap <- as.numeric(mOd$ProportionalOverlap)
mOd$Within_Donor <- mOd$SDprot == mOd$MDprot
mOd$SDprot <- factor(mOd$SDprot, levels=c("CW","H1","Mito","PGP","11a"))
mOd$MDprot <- factor(mOd$MDprot, levels=c("CW","H1","Mito","PGP","11a"))
ggplot(mOd, aes(x=Within_Donor, y=FisherSig, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mO$SDprot)])
ggplot(mOd, aes(x=Within_Donor, y=FisherSig, fill=MDprot)) +
  geom_boxplot() +
  facet_grid(~MDprot) +
  scale_fill_manual(values=DonorCols[levels(mOd$MDprot)])

p1 <- ggplot(mOd, aes(x=Within_Donor, y=ProportionalOverlap, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mOd$SDprot)])
p2 <- ggplot(mOd, aes(x=Within_Donor, y=ProportionalOverlap, fill=MDprot)) +
  geom_boxplot() +
  facet_grid(~MDprot) +
  scale_fill_manual(values=DonorCols[levels(mOd$MDprot)])
p1 + p2


## Both?
tmp3 <- list()
myThings <- c("CW","CWsd","PGP","PGPsd","Mito","Mitosd","H1","H1sd","11a","11asd")
for(i in myThings){
  tmp3[[i]] <- c(tmp[paste0("sigUp",i)][[1]], tmpD[paste0("sigDown",i)][[1]])
}
## The downs:
myOverlap6 <- matrix(ncol=10,nrow=10)
colnames(myOverlap6) <- names(tmp3)
rownames(myOverlap6) <- names(tmp3)
myOverlap7 <- matrix(ncol=10,nrow=10)
colnames(myOverlap7) <- names(tmp3)
rownames(myOverlap7) <- names(tmp3)
for(i in names(tmp3)){
  for(j in names(tmp3)){
    myOverlap6[i,j] <- sum(tmp3[[i]] %in% tmp3[[j]])
    myOverlap7[i,j] <- -log10(fisher.test(matrix(ncol=2, data=c(sum(sharedGenes %ni% c(tmp3[[i]],tmp3[[j]])),
                                                                sum(tmp3[[i]]%ni%tmp3[[j]]),
                                                                sum(tmp3[[j]]%ni%tmp3[[i]]),
                                                                sum(tmp3[[i]] %in% tmp3[[j]]))))$p.value)
  }
}

mOb <- data.frame(SDprot = c("Fake"), MDprot = c("Fake"), FisherSig = c(1), Overlap = 1, ProportionalOverlap = 1)
for(i in c("CW","PGP","Mito","H1","11a")){
  for(j in c("CW","PGP","Mito","H1","11a")){
    mOb <- rbind(mOb, c(i, j, myOverlap7[paste0(i,"sd"),paste0(j)],
                        myOverlap6[paste0(i,"sd"),paste0(j)],
                        myOverlap6[paste0(i,"sd"),paste0(j)]/min(length(tmp3[[paste0(j)]]),length(tmp3[[paste0(i,"sd")]]))))
  }
}
mOb <- mOb[-1,]
mOb$FisherSig <- as.numeric(mOb$FisherSig)
mOb$Overlap <- as.numeric(mOb$Overlap)
mOb$ProportionalOverlap <- as.numeric(mOb$ProportionalOverlap)
mOb$Within_Donor <- mOb$SDprot == mOb$MDprot
mOb$SDprot <- factor(mOb$SDprot, levels=c("CW","H1","Mito","PGP","11a"))
mOb$MDprot <- factor(mOb$MDprot, levels=c("CW","H1","Mito","PGP","11a"))
ggplot(mOb, aes(x=Within_Donor, y=FisherSig, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mO$SDprot)])
# ggplot(mOb, aes(x=Within_Donor, y=FisherSig, fill=MDprot)) +
#   geom_boxplot() +
#   facet_grid(~MDprot) +
#   scale_fill_manual(values=DonorCols[levels(mOb$MDprot)])

p1 <- ggplot(mOb, aes(x=Within_Donor, y=ProportionalOverlap, fill=SDprot)) +
  geom_boxplot() +
  facet_grid(~SDprot) +
  scale_fill_manual(values=DonorCols[levels(mOb$SDprot)]) +
  ggtitle("Proportional Overlap of Genes Significantly Dysregulated by VPA\nBetween Each Single-Donor and Each Donor in Multi-Donor Chimeroids") +
  theme(plot.title=element_text(face="bold", hjust=0.5))
p2 <- ggplot(mOb, aes(x=Within_Donor, y=ProportionalOverlap, fill=MDprot)) +
  geom_boxplot() +
  facet_grid(~MDprot) +
  scale_fill_manual(values=DonorCols[levels(mOb$MDprot)])
p1 + p2
p1

ggplot(mOb, aes(x=SDprot, y=ProportionalOverlap, color=MDprot)) +
  geom_point(size=5) +
  scale_color_manual(values=DonorCols[levels(mOb$MDprot)]) +
  ggtitle("Proportional Overlap of Genes Significantly Dysregulated by VPA\nBetween Each Single-Donor and Each Donor in Multi-Donor Chimeroids") +
  theme(plot.title=element_text(face="bold", hjust=0.5))

## Genes that are *completely* unique to CW, and exist in BOTH SD and MD: GRM3, LPCAT2
CWgenes <- sigUpCW[sigUpCW %in% sigUpCWsd & sigUpCW %ni% c(sigUpPGPsd, sigUpMitosd, sigUpH1sd, sigUp11asd,
                                                           sigUpPGP, sigUpMito, sigUpH1, sigUp11a)]
## Genes that are *completely* unique to PGP1, and exist in BOTH SD and MD: None
PGPgenes <- sigUpPGP[sigUpPGP %in% sigUpPGPsd & sigUpPGP %ni% c(sigUpCWsd, sigUpMitosd, sigUpH1sd, sigUp11asd,
                                                           sigUpCW, sigUpMito, sigUpH1, sigUp11a)]
## Genes that are *completely* unique to Mito, and exist in BOTH SD and MD: ZNF208, ZNF676
Mitogenes <- sigUpMito[sigUpMito %in% sigUpMitosd & sigUpMito %ni% c(sigUpPGPsd, sigUpCWsd, sigUpH1sd, sigUp11asd,
                                                           sigUpPGP, sigUpCW, sigUpH1, sigUp11a)]
## Genes that are *completely* unique to Mito, and exist in BOTH SD and MD: ZNF208, ZNF676
H1genes <- sigUpMito[sigUpH1 %in% sigUpH1sd & sigUpH1 %ni% c(sigUpPGPsd, sigUpCWsd, sigUpMitosd, sigUp11asd,
                                                                     sigUpPGP, sigUpCW, sigUpMito, sigUp11a)]
## Genes that are *completely* unique to Mito, and exist in BOTH SD and MD: ZNF208, ZNF676
a11genes <- sigUp11a[sigUp11a %in% sigUp11asd & sigUp11a %ni% c(sigUpPGPsd, sigUpCWsd, sigUpH1sd, sigUpMitosd,
                                                                     sigUpPGP, sigUpCW, sigUpH1, sigUpMito)]

sigUpCW[sigUpCW %in% sigUpPGPsd & sigUpCW %ni% c(sigUpCWsd, sigUpMitosd, sigUpH1sd, sigUp11asd,
                                                sigUpPGP, sigUpMito, sigUpH1, sigUp11a)]


ggplot(mO, aes(x=Within_Donor, y=FisherSig, fill=MDprot)) +
  geom_boxplot() +
  facet_grid(~MDprot)
ggplot(mO, aes(x=Within_Donor, y=FisherSig, color=SDprot)) +
  geom_jitter()


ggplot(myFC, aes(x=CW_FC, y=PGP_FC)) +
  geom_point()

ggplot(myFC, aes(x=PGP_FC, y=PGPsd_FC)) +
  geom_point()

ggplot(myFC, aes(x=CWsd_FC, y=CW_FC)) +
  geom_point()




PGP <- Correlate_Results(CPN_PGP1, CPN_PGP1SD)
Mito <- Correlate_Results(CPN_Mito, CPN_MitoSD)
PGP_Mito_MD <- Correlate_Results(CPN_PGP1, CPN_Mito)
PGP_Mito_SD <- Correlate_Results(CPN_PGP1SD, CPN_MitoSD)

Make_Heatmap(myGenes=getGenes(global_MvS),
             Treatments=c("Control"),
             IsGlobal = TRUE,
             Protocols = c("SingleDonor","MultiDonor"),
             myTitle = "Top DEGs Between Single-Donor and Multi-Donor\nGlobally in Control-Treated Chimeroids")

# Make_Heatmap(myGenes=getGenes(global_MvS),
#              Treatments=c("Control"),
#              IsGlobal = FALSE,
#              CellTypes = c("aRG","oRG","IP","CPN","CFuPN","CR","ChPL","IN"),
#              Protocols = c("SingleDonor","MultiDonor"),
#              myTitle = "Top DEGs Between Single-Donor and Multi-Donor\nGlobally in All Cell Types of Control-Treated Chimeroids")

Make_Heatmap(myGenes=getGenes(aRG_MvS),
             Treatments=c("Control"),
             CellTypes=c("aRG"),
             Protocols = c("SingleDonor","MultiDonor"),
             myTitle = "Top DEGs Between Single-Donor and Multi-Donor\nin Control-Treated aRGs")

Make_Heatmap(myGenes=getGenes(oRG_MvS),
             Treatments=c("Control"),
             CellTypes=c("oRG"),
             Protocols = c("SingleDonor","MultiDonor"),
             myTitle = "Top DEGs Between Single-Donor and Multi-Donor\nin Control-Treated oRGs")

Make_Heatmap(myGenes=getGenes(IP_MvS),
             Treatments=c("Control"),
             CellTypes=c("IP"),
             Protocols = c("SingleDonor","MultiDonor"),
             myTitle = "Top DEGs Between Single-Donor and Multi-Donor\nin Control-Treated IPs")


## Functionize everything?
VPA_Res <- list()
VPA_Res[["Global"]] <- list()
for(donor in c("CW","H1","Mito210","PGP1","11a")){
  VPA_Res[["Global"]][[donor]] <- list()
  for(protocol in c("SingleDonor","MultiDonor")){
    VPA_Res[["Global"]][[donor]][[protocol]] <- GetDE_Results(Treatments=c("Control","VPA"),
                                                              CellTypes=celltype,
                                                              IsGlobal = TRUE,
                                                              Genotypes=donor,
                                                              Protocols=protocol,
                                                              myEffect="Treatment",
                                                              addCovs="Batch")
  }
}
for(celltype in c("aRG","oRG","IP","CFuPN","CPN")){
  VPA_Res[[celltype]] <- list()
  for(donor in c("CW","H1","Mito210","PGP1","11a")){
    VPA_Res[[celltype]][[donor]] <- list()
    for(protocol in c("SingleDonor","MultiDonor")){
      VPA_Res[[celltype]][[donor]][[protocol]] <- GetDE_Results(Treatments=c("Control","VPA"),
                                                                CellTypes=celltype,
                                                                Genotypes=donor,
                                                                Protocols=protocol,
                                                                myEffect="Treatment",
                                                                addCovs="Batch")
    }
  }
}
saveRDS(VPA_Res, "/Users/tfaits/Documents/Chimera/Results/VPA_DEGs_SD_MD_FullObject.RDS")

EtOH_Res <- list()
EtOH_Res[["Global"]] <- list()
for(donor in c("CW","H1","Mito210","PGP1","11a")){
  EtOH_Res[["Global"]][[donor]] <- list()
  for(protocol in c("SingleDonor","MultiDonor")){
    EtOH_Res[["Global"]][[donor]][[protocol]] <- GetDE_Results(Treatments=c("Control","EtOH"),
                                                              CellTypes=celltype,
                                                              IsGlobal = TRUE,
                                                              Genotypes=donor,
                                                              Protocols=protocol,
                                                              myEffect="Treatment",
                                                              addCovs="Batch")
  }
}
for(celltype in c("aRG","oRG","IP","CFuPN","CPN")){
  print(paste("Starting",celltype))
  EtOH_Res[[celltype]] <- list()
  for(donor in c("CW","H1","Mito210","PGP1","11a")){
    print(paste("Starting",celltype,donor))
    EtOH_Res[[celltype]][[donor]] <- list()
    for(protocol in c("SingleDonor","MultiDonor")){
      print(paste("Starting",celltype,donor,protocol))
      yesBatch <- "Batch"
      if(donor=="11a" & protocol=="MultiDonor"){
        yesBatch <- ""
      }
      EtOH_Res[[celltype]][[donor]][[protocol]] <- GetDE_Results(Treatments=c("Control","EtOH"),
                                                                CellTypes=celltype,
                                                                Genotypes=donor,
                                                                Protocols=protocol,
                                                                myEffect="Treatment",
                                                                addCovs=yesBatch)
    }
  }
}
for(celltype in c("aRG","oRG","IP","CFuPN","CPN")){
  print(paste("Starting",celltype))
  EtOH_Res[[celltype]]$Global <- list()
  for(protocol in c("SingleDonor","MultiDonor")){
    print(paste("Starting",celltype,"Global",protocol))
    EtOH_Res[[celltype]]$Global[[protocol]] <- GetDE_Results(Treatments=c("Control","EtOH"),
                                                             CellTypes=celltype,
                                                             Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                             Protocols=protocol,
                                                             myEffect="Treatment",
                                                             addCovs="Batch")
  }
}
EtOH_Res$Global$Global <- list()
EtOH_Res$Global$Global$MultiDonor <- GetDE_Results(Treatments=c("Control","EtOH"),
                                                   IsGlobal=TRUE,
                                                   Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                   Protocols="MultiDonor",
                                                   myEffect="Treatment",
                                                   addCovs="Batch + Genotype")
EtOH_Res$Global$Global$SingleDonor <- GetDE_Results(Treatments=c("Control","EtOH"),
                                                   IsGlobal=TRUE,
                                                   Genotypes=c("CW","H1","Mito210","PGP1","11a"),
                                                   Protocols="SingleDonor",
                                                   myEffect="Treatment",
                                                   addCovs="Batch + Genotype")
saveRDS(EtOH_Res, "/Users/tfaits/Documents/Chimera/Results/EtOH_DEGs_SD_MD_FullObject.RDS")




getOverlap <- function(Donor1="CW", Protocol1="MultiDonor",
                       Donor2="CW", Protocol2="SingleDonor", CellType="Global",
                       doPrint=FALSE){
  upG1 <- rownames(VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results)[VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$padj<0.1 &
                                                                       ! is.na(VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$padj) &
                                                                       VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$log2FoldChange > 0.5]
  upG2 <- rownames(VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results)[VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$padj<0.1 &
                                                                         ! is.na(VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$padj) &
                                                                         VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$log2FoldChange > 0.5]
  downG1 <- rownames(VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results)[VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$padj<0.1 &
                                                                           ! is.na(VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$padj) &
                                                                           VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results$log2FoldChange < -0.5]
  downG2 <- rownames(VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results)[VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$padj<0.1 &
                                                                           ! is.na(VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$padj) &
                                                                           VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results$log2FoldChange < -0.5]
  if(doPrint){
    print(paste("Genes up in",Donor1,Protocol1,CellType,":",length(upG1)))
    print(paste("Genes down in",Donor1,Protocol1,CellType,":",length(downG1)))
    print(paste("Genes up in",Donor2,Protocol2,CellType,":",length(upG2)))
    print(paste("Genes down in",Donor2,Protocol2,CellType,":",length(downG2)))
    print(paste("Genes in common:", sum(c(upG1 %in% upG2, downG1 %in% downG2))))
  }
  return(c("Overlap"=sum(c(upG1 %in% upG2, downG1 %in% downG2)), "G1"=length(upG1)+length(downG1), "G2"=length(upG2)+length(downG2)))
}
getFisher <- function(Donor1="CW", Protocol1="MultiDonor",
                      Donor2="CW", Protocol2="SingleDonor", CellType="Global",
                      doPrint=FALSE){
  return(Fisher_Results(VPA_Res[[CellType]][[Donor1]][[Protocol1]],
                        VPA_Res[[CellType]][[Donor2]][[Protocol2]]))
}
getCorrelation <- function(Donor1="CW", Protocol1="MultiDonor",
                           Donor2="CW", Protocol2="SingleDonor", CellType="Global",
                           Metric="log2FoldChange"){
  res1 <- VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results
  res2 <- VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results
  res1 <- res1[!is.na(res1$padj),]
  res2 <- res2[!is.na(res2$padj),]
  shGen <- rownames(res1)[rownames(res1) %in% rownames(res2)]
  res1 <- res1[shGen,]
  res2 <- res2[shGen,]
  return(cor(res1[[Metric]], res2[[Metric]]))
}

GlobalOverlap <- matrix(ncol=5, nrow=5)
colnames(GlobalOverlap) <- rownames(GlobalOverlap) <- c("CW","H1","Mito210","PGP1","11a")
GlobalCorFC <- matrix(ncol=5, nrow=5)
colnames(GlobalCorFC) <- rownames(GlobalCorFC) <- c("CW","H1","Mito210","PGP1","11a")
GlobalCorPval <- matrix(ncol=5, nrow=5)
colnames(GlobalCorPval) <- rownames(GlobalCorPval) <- c("CW","H1","Mito210","PGP1","11a")
GlobalFisher <- matrix(ncol=5, nrow=5)
colnames(GlobalFisher) <- rownames(GlobalFisher) <- c("CW","H1","Mito210","PGP1","11a")

for(i in rownames(GlobalOverlap)){
  for(j in colnames(GlobalOverlap)){
    tmpOverlap <- getOverlap(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType="Global")
    GlobalOverlap[i,j] <- tmpOverlap[1]/min(tmpOverlap[2:3])
    GlobalCorFC[i,j] <- getCorrelation(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType="Global",Metric="log2FoldChange")
    GlobalCorPval[i,j] <- getCorrelation(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType="Global",Metric="padj")
    GlobalFisher[i,j] <- getFisher(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType="Global")
  }
}



gf <- melt(GlobalFisher)
colnames(gf) <- c("MultiDonor","SingleDonor","pval")
gf$logP <- -log10(gf$pval)
gf$WithinDonor <- gf$MultiDonor==gf$SingleDonor
ggplot(gf, aes(x=WithinDonor, y=logP, fill=SingleDonor)) +
  geom_boxplot() +
  scale_fill_manual(values=DonorCols[levels(gf$SingleDonor)]) +
  facet_grid(~SingleDonor) +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Fisher's Exact Test for overlapping significant genes between SD and MDs\nGrouped by the SD donor, highlighting the homogenous comparison\n(Global)")



gcfc <- melt(GlobalCorFC)
colnames(gcfc) <- c("MultiDonor","SingleDonor","R2")
gcfc$WithinDonor <- gcfc$MultiDonor==gcfc$SingleDonor
ggplot(gcfc, aes(x=WithinDonor, y=R2, fill=SingleDonor)) +
  geom_boxplot() +
  scale_fill_manual(values=DonorCols[levels(gcfc$SingleDonor)]) +
  facet_grid(~SingleDonor) +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Correlations of VPA-induced log2FC between SD and MDs\nGrouped by the SD donor, highlighting the homogenous comparison\n(Global)")

gcpv <- melt(GlobalCorPval)
colnames(gcpv) <- c("MultiDonor","SingleDonor","R2")
gcpv$WithinDonor <- gcpv$MultiDonor==gcpv$SingleDonor
ggplot(gcpv, aes(x=WithinDonor, y=R2, fill=SingleDonor)) +
  geom_boxplot() +
  scale_fill_manual(values=DonorCols[levels(gcpv$SingleDonor)]) +
  facet_grid(~SingleDonor) +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Correlations of significance values of VPA-induced DEGs between SD and MDs\nGrouped by the SD donor, highlighting the homogenous comparison\n(Global)")

make_Comparison_Plot <- function(inMat, myTitle="NoTitle?"){
  meltMat <- melt(inMat)
  colnames(meltMat) <- c("MultiDonor","SingleDonor","value")
  meltMat$WithinDonor <- meltMat$MultiDonor==meltMat$SingleDonor
  myPlot <- ggplot(meltMat, aes(x=WithinDonor, y=value, fill=SingleDonor)) +
    geom_boxplot() +
    scale_fill_manual(values=DonorCols[levels(meltMat$SingleDonor)]) +
    facet_grid(~SingleDonor) +
    theme(plot.title=element_text(face="bold", hjust=0.5),
          panel.background = element_blank()) +
    ggtitle(myTitle) +
    geom_point(aes(color=MultiDonor), size=3) +
    scale_color_manual(values=DonorCols[levels(meltMat$MultiDonor)])
  return(myPlot)
}


myOverlaps <- list()
myCorFCs <- list()
myCorPvals <- list()
myFishers <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  myOverlaps[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myOverlaps[[celltype]]) <- rownames(myOverlaps[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  myCorFCs[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myCorFCs[[celltype]]) <- rownames(myCorFCs[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  myCorPvals[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myCorPvals[[celltype]]) <- rownames(myCorPvals[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  myFishers[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myFishers[[celltype]]) <- rownames(myFishers[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myOverlaps[[celltype]])){
    for(j in colnames(myOverlaps[[celltype]])){
      tmpOverlap <- getOverlap(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType=celltype)
      myOverlaps[[celltype]][i,j] <- tmpOverlap[1]/min(tmpOverlap[2:3])
      myCorFCs[[celltype]][i,j] <- getCorrelation(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType=celltype,Metric="log2FoldChange")
      myCorPvals[[celltype]][i,j] <- getCorrelation(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType=celltype,Metric="padj")
      myFishers[[celltype]][i,j] <- getFisher(Donor1=i, Donor2=j, Protocol1="MultiDonor", Protocol2="SingleDonor", CellType=celltype)
    }
  }
}

make_Comparison_Plot(myOverlaps$Global, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nGlobal Response to VPA")
make_Comparison_Plot(myOverlaps$CFuPN, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn CFuPN Response to VPA")
make_Comparison_Plot(myOverlaps$aRG, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn aRG Response to VPA")
make_Comparison_Plot(myOverlaps$oRG, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn oRG Response to VPA")
make_Comparison_Plot(myOverlaps$IP, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn IP Response to VPA")
make_Comparison_Plot(myOverlaps$CPN, myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn CPN Response to VPA")

make_Comparison_Plot(-log10(myFishers$Global), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs: Global Response to VPA")
make_Comparison_Plot(-log10(myFishers$CFuPN), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In CFuPN Response to VPA")
make_Comparison_Plot(-log10(myFishers$aRG), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In aRG Response to VPA")
make_Comparison_Plot(-log10(myFishers$oRG), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In oRG Response to VPA")
make_Comparison_Plot(-log10(myFishers$IP), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In IP Response to VPA")
make_Comparison_Plot(-log10(myFishers$CPN), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In CPN Response to VPA")

make_Comparison_Plot(myCorFCs$Global, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nGlobal Response to VPA")
make_Comparison_Plot(myCorFCs$CFuPN, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn CFuPN Response to VPA")
make_Comparison_Plot(myCorFCs$aRG, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn aRG Response to VPA")
make_Comparison_Plot(myCorFCs$oRG, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn oRG Response to VPA")
make_Comparison_Plot(myCorFCs$IP, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn IP Response to VPA")
make_Comparison_Plot(myCorFCs$CPN, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn CPN Response to VPA")

# make_Comparison_Plot(myCorPvals$CFuPN, myTitle="Pval Correlation Between Single-Donor and Multi-Donor DEGs\nIn CFuPN Response to VPA")
# make_Comparison_Plot(myCorPvals$aRG, myTitle="Pval Correlation Between Single-Donor and Multi-Donor DEGs\nIn aRG Response to VPA")
# make_Comparison_Plot(myCorPvals$oRG, myTitle="Pval Correlation Between Single-Donor and Multi-Donor DEGs\nIn oRG Response to VPA")
# make_Comparison_Plot(myCorPvals$IP, myTitle="Pval Correlation Between Single-Donor and Multi-Donor DEGs\nIn IP Response to VPA")
# make_Comparison_Plot(myCorPvals$CPN, myTitle="Pval Correlation Between Single-Donor and Multi-Donor DEGs\nIn CPN Response to VPA")

getUniqueGenes <- function(celltype="CFuPN", MultiDonor="CW", SingleDonor="PGP1",
                           minFC=0.5, maxPadj=0.1){
  myMultis <- c("CW","H1","Mito210","PGP1","11a")
  myMultis <- myMultis[myMultis != MultiDonor]
  mySingles <- c("CW","H1","Mito210","PGP1","11a")
  mySingles <- mySingles[mySingles != SingleDonor]
  mG <- VPA_Res[[celltype]][[MultiDonor]]$MultiDonor$Results
  mG <- rownames(mG)[!is.na(mG$padj) & mG$padj < maxPadj & abs(mG$log2FoldChange) > minFC]
  sG <- VPA_Res[[celltype]][[SingleDonor]]$SingleDonor$Results
  sG <- rownames(sG)[!is.na(sG$padj) & sG$padj < maxPadj & abs(sG$log2FoldChange) > minFC]
  otherGenes <- c()
  for(i in myMultis){
    tG <- VPA_Res[[celltype]][[i]]$MultiDonor$Results
    tG <- rownames(tG)[tG$padj < maxPadj & !is.na(tG$padj) & abs(tG$log2FoldChange) > minFC]
    otherGenes <- unique(c(otherGenes,tG))
  }
  for(i in mySingles){
    tG <- VPA_Res[[celltype]][[i]]$SingleDonor$Results
    tG <- rownames(tG)[tG$padj < maxPadj & !is.na(tG$padj) & abs(tG$log2FoldChange) > minFC]
    otherGenes <- unique(c(otherGenes,tG))
  }
  commonGenes <- mG[mG %in% sG]
  uniqueGenes <- commonGenes[commonGenes %ni% otherGenes]
  return(uniqueGenes)
}


## Try things flipped; split the facets by the multidonors?
make_Comparison_Plot(t(myOverlaps$CFuPN), myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn CFuPN Response to VPA")
make_Comparison_Plot(t(myOverlaps$aRG), myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn aRG Response to VPA")
make_Comparison_Plot(t(myOverlaps$oRG), myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn oRG Response to VPA")
make_Comparison_Plot(t(myOverlaps$IP), myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn IP Response to VPA")
make_Comparison_Plot(t(myOverlaps$CPN), myTitle="Overlaps Between Single-Donor and Multi-Donor DEGs\nIn CPN Response to VPA")

make_Comparison_Plot(t(-log10(myFishers$Global)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In Global Response to VPA")
make_Comparison_Plot(t(-log10(myFishers$CFuPN)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In CFuPN Response to VPA")
make_Comparison_Plot(t(-log10(myFishers$aRG)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In aRG Response to VPA")
make_Comparison_Plot(t(-log10(myFishers$oRG)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In oRG Response to VPA")
make_Comparison_Plot(t(-log10(myFishers$IP)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In IP Response to VPA")
make_Comparison_Plot(t(-log10(myFishers$CPN)), myTitle="Significance of Overlaps Between Single-Donor\nand Multi-Donor DEGs In CPN Response to VPA")

make_Comparison_Plot(myCorFCs$CFuPN, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn CFuPN Response to VPA")
make_Comparison_Plot(myCorFCs$aRG, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn aRG Response to VPA")
make_Comparison_Plot(myCorFCs$oRG, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn oRG Response to VPA")
make_Comparison_Plot(myCorFCs$IP, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn IP Response to VPA")
make_Comparison_Plot(myCorFCs$CPN, myTitle="FoldChange Correlation Between Single-Donor and Multi-Donor DEGs\nIn CPN Response to VPA")

## Try some KS tests. Rather than overlap of the significant genes, look at rankings
runKS_Test <- function(queryRes, refRes, numGenes=200, rankMetric="log2FoldChange", Decreasing=TRUE){
  ##Order the reference results by p-value, then take the top <numGenes> genes:
  refRes <- refRes[order(refRes$pvalue),]
  myGeneList <- rownames(refRes)[1:200]
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  queryRanks <- which(rownames(queryRes) %in% myGeneList)
  return(ks.test(queryRanks, 1:nrow(queryRes)))
}

plotKS <- function(queryRes, refRes, numGenes=200, rankMetric="log2FoldChange", Decreasing=TRUE){
  refRes <- refRes[order(refRes$pvalue),]
  myGeneList <- rownames(refRes)[1:200]
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  queryRanks <- which(rownames(queryRes) %in% myGeneList)
  myP <- hist(queryRanks, breaks=nrow(queryRes))
  return(myP)
}

myKS <- list()
myKSt <- list()
myUniqueGenes <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  myKS[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myKS[[celltype]]) <- rownames(myKS[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  myKSt[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myKSt[[celltype]]) <- rownames(myKSt[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  myUniqueGenes[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myUniqueGenes[[celltype]]) <- rownames(myUniqueGenes[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myKS[[celltype]])){
    for(j in colnames(myKS[[celltype]])){
      myKS[[celltype]][i,j] <- runKS_Test(refRes = VPA_Res[[celltype]][[i]]$MultiDonor$Results,
                                          queryRes = VPA_Res[[celltype]][[j]]$SingleDonor$Results,
                                          rankMetric = "pvalue")$statistic
      myKSt[[celltype]][j,i] <- runKS_Test(refRes = VPA_Res[[celltype]][[j]]$MultiDonor$Results,
                                          queryRes = VPA_Res[[celltype]][[i]]$SingleDonor$Results,
                                          rankMetric = "pvalue")$statistic
      myUniqueGenes[[celltype]][i,j] <- length(getUniqueGenes(celltype=celltype, MultiDonor=i, SingleDonor=j))
    }
  }
}

getUniqueGenes2 <- function(celltype="CFuPN", MultiDonor="CW", SingleDonor="PGP1", numGenes=200){
  myMultis <- c("CW","H1","Mito210","PGP1","11a")
  myMultis <- myMultis[myMultis != MultiDonor]
  mySingles <- c("CW","H1","Mito210","PGP1","11a")
  mySingles <- mySingles[mySingles != SingleDonor]
  mG <- VPA_Res[[celltype]][[MultiDonor]]$MultiDonor$Results
  mG <- mG[order(mG$pvalue),]
  mG <- rownames(mG)[1:numGenes]
  sG <- VPA_Res[[celltype]][[SingleDonor]]$SingleDonor$Results
  sG <- sG[order(sG$pvalue),]
  sG <- rownames(sG)[1:numGenes]
  otherGenes <- c()
  for(i in myMultis){
    tG <- VPA_Res[[celltype]][[i]]$MultiDonor$Results
    tG <- tG[order(tG$pvalue),]
    tG <- rownames(tG)[1:numGenes]
    otherGenes <- unique(c(otherGenes,tG))
  }
  for(i in mySingles){
    tG <- VPA_Res[[celltype]][[i]]$SingleDonor$Results
    tG <- rownames(tG)[1:numGenes]
    tG <- rownames(tG)[1:numGenes]
    otherGenes <- unique(c(otherGenes,tG))
  }
  commonGenes <- mG[mG %in% sG]
  uniqueGenes <- commonGenes[commonGenes %ni% otherGenes]
  return(uniqueGenes)
}

myUniqueGenes2 <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  myUniqueGenes2[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(myUniqueGenes2[[celltype]]) <- rownames(myUniqueGenes2[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myKS[[celltype]])){
    for(j in colnames(myKS[[celltype]])){
      myUniqueGenes2[[celltype]][i,j] <- length(getUniqueGenes2(celltype=celltype, MultiDonor=i, SingleDonor=j, numGenes = 100))
    }
  }
}




mUG2 <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  mUG2[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(mUG2[[celltype]]) <- rownames(mUG2[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myKS[[celltype]])){
    for(j in colnames(myKS[[celltype]])){
      mUG2[[celltype]][i,j] <- length(getUniqueGenes2(celltype=celltype, MultiDonor=i, SingleDonor=j, numGenes = 500))
    }
  }
}
mUG2$Global

for(idx in seq(20,100,5)){
  mUG2 <- list()
  for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
    mUG2[[celltype]] <- matrix(ncol=5, nrow=5)
    colnames(mUG2[[celltype]]) <- rownames(mUG2[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
    for(i in rownames(myKS[[celltype]])){
      for(j in colnames(myKS[[celltype]])){
        mUG2[[celltype]][i,j] <- length(getUniqueGenes2(celltype=celltype, MultiDonor=i, SingleDonor=j, numGenes = idx))
      }
    }
  }
  print(idx)
  print(mUG2$Global)
}


mUG <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  mUG[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(mUG[[celltype]]) <- rownames(mUG[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myKS[[celltype]])){
    for(j in colnames(myKS[[celltype]])){
      mUG[[celltype]][i,j] <- length(getUniqueGenes(celltype=celltype, MultiDonor=i, SingleDonor=j, minFC = 0.5, maxPadj = 0.01))
    }
  }
}
mUG$Global



getOverlap2 <- function(Donor1="CW", Protocol1="MultiDonor",
                       Donor2="CW", Protocol2="SingleDonor", CellType="Global",
                       numGenes=200, myMetric = "pvalue", Decreasing=FALSE){
  tmp1 <- VPA_Res[[CellType]][[Donor1]][[Protocol1]]$Results
  tmp1 <- tmp1[order(abs(tmp1[,myMetric])),]
  upG1 <- rownames(tmp1)[1:numGenes]
  tmp2 <- VPA_Res[[CellType]][[Donor2]][[Protocol2]]$Results
  tmp2 <- tmp2[order(abs(tmp2[,myMetric])),]
  upG2 <- rownames(tmp2)[1:numGenes]
  return(sum(upG1 %in% upG2))
}


mOV <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  mOV[[celltype]] <- matrix(ncol=5, nrow=5)
  colnames(mOV[[celltype]]) <- rownames(mOV[[celltype]]) <- c("CW","H1","Mito210","PGP1","11a")
  for(i in rownames(myKS[[celltype]])){
    for(j in colnames(myKS[[celltype]])){
      mOV[[celltype]][i,j] <- getOverlap2(Donor1=i, Donor2=j, CellType=celltype, numGenes = 1000)
    }
  }
}
mOV$Global

printSigs <- function(){
  for(k in c("Global","aRG","oRG","IP","CFuPN","CPN")){
    for(j in c("SingleDonor","MultiDonor")){
      for(i in c("CW","H1","Mito210","PGP1","11a")){
        print(paste(k, j, i, ":", sum(VPA_Res[[k]][[i]][[j]]$Results$padj<0.1 &
                                        !is.na(VPA_Res[[k]][[i]][[j]]$Results$padj) &
                                        abs(VPA_Res[[k]][[i]][[j]]$Results$log2FoldChange)>0.5)))
      }
    }
  }
}

make_Comparison_Plot(t(mOV$Global), myTitle="New Overlap - Global")
make_Comparison_Plot(mOV$aRG, myTitle="New Overlap - aRG")
make_Comparison_Plot(mOV$oRG, myTitle="New Overlap - oRG")
make_Comparison_Plot(mOV$IP, myTitle="New Overlap - IP")
make_Comparison_Plot(mOV$CFuPN, myTitle="New Overlap - CFuPN")
make_Comparison_Plot(mOV$CPN, myTitle="New Overlap - CPN")




make_Comparison_Plot(myKS$Global, myTitle="KS Test\nGlobal - SD as Query")
make_Comparison_Plot(myKS$aRG, myTitle="KS Test\naRG - SD as Query")
make_Comparison_Plot(myKS$oRG, myTitle="KS Test\noRG - SD as Query")
make_Comparison_Plot(myKS$IP, myTitle="KS Test\nIP - SD as Query")
make_Comparison_Plot(myKS$CFuPN, myTitle="KS Test\nCFuPN - SD as Query")
make_Comparison_Plot(myKS$CPN, myTitle="KS Test\nCPN - SD as Query")

make_Comparison_Plot(myKSt$Global, myTitle="KS Test\nGlobal - MD as Query")
make_Comparison_Plot(myKSt$aRG, myTitle="KS Test\naRG - MD as Query")
make_Comparison_Plot(myKSt$oRG, myTitle="KS Test\noRG - MD as Query")
make_Comparison_Plot(myKSt$IP, myTitle="KS Test\nIP - MD as Query")
make_Comparison_Plot(myKSt$CFuPN, myTitle="KS Test\nCFuPN - MD as Query")
make_Comparison_Plot(myKSt$CPN, myTitle="KS Test\nCPN - MD as Query")




## Some plots to see it?
par(mfrow=c(5,1))
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$CW$SingleDonor$Results, rankMetric = "pvalue", Decreasing=FALSE)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$H1$SingleDonor$Results, rankMetric = "pvalue", Decreasing=FALSE)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$Mito210$SingleDonor$Results, rankMetric = "pvalue", Decreasing=FALSE)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$PGP1$SingleDonor$Results, rankMetric = "pvalue", Decreasing=FALSE)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN[["11a"]]$SingleDonor$Results, rankMetric = "pvalue", Decreasing=FALSE)

plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$CW$SingleDonor$Results)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$H1$SingleDonor$Results)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$Mito210$SingleDonor$Results)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN$PGP1$SingleDonor$Results)
plotKS(queryRes=VPA_Res$CPN[["11a"]]$MultiDonor$Results, refRes = VPA_Res$CPN[["11a"]]$SingleDonor$Results)





## Try calculating a z-score for the FC of each gene for each donor within a protocol?
FCmats <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  FCmats[[celltype]] <- list()
  sharedGenes <- rownames(VPA_Res[[celltype]]$CW$SingleDonor$Results)
  sharedGenes <- sharedGenes[sharedGenes %in% rownames(VPA_Res[[celltype]]$H1$SingleDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]]$Mito210$SingleDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]]$PGP1$SingleDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]][["11a"]]$SingleDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]]$H1$MultiDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]]$Mito210$MultiDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]]$PGP1$MultiDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]][["11a"]]$MultiDonor$Results) &
                               sharedGenes %in% rownames(VPA_Res[[celltype]][["CW"]]$MultiDonor$Results)]
  for(Protocol in c("SingleDonor","MultiDonor")){
    FCmats[[celltype]][[Protocol]] <- matrix(ncol=5, nrow=length(sharedGenes))
    rownames(FCmats[[celltype]][[Protocol]]) <- sharedGenes
    colnames(FCmats[[celltype]][[Protocol]]) <- c("CW","H1","Mito210","PGP1","11a")
    FCmats[[celltype]][[Protocol]][,"CW"] <- VPA_Res[[celltype]]$CW[[Protocol]]$Results[sharedGenes,"log2FoldChange"]
    FCmats[[celltype]][[Protocol]][,"H1"] <- VPA_Res[[celltype]]$H1[[Protocol]]$Results[sharedGenes,"log2FoldChange"]
    FCmats[[celltype]][[Protocol]][,"Mito210"] <- VPA_Res[[celltype]]$Mito210[[Protocol]]$Results[sharedGenes,"log2FoldChange"]
    FCmats[[celltype]][[Protocol]][,"PGP1"] <- VPA_Res[[celltype]]$PGP1[[Protocol]]$Results[sharedGenes,"log2FoldChange"]
    FCmats[[celltype]][[Protocol]][,"11a"] <- VPA_Res[[celltype]][["11a"]][[Protocol]]$Results[sharedGenes,"log2FoldChange"]
    tmpMat <- FCmats[[celltype]][[Protocol]]
    for(i in 1:nrow(tmpMat)){
      tmpMat[i,] <- (tmpMat[i,]-mean(tmpMat[i,]))/sd(tmpMat[i,])
    }
    FCmats[[celltype]][[Protocol]] <- tmpMat
  }
}

## Let's put a hold on comparisons with the Single-Donor, and just look at what makes each donor unique:
getUniqueGenes3 <- function(celltype="CFuPN", Donor="CW", Protocol="MultiDonor", minFC=0.5, maxPadj=0.1){
  myMultis <- c("CW","H1","Mito210","PGP1","11a")
  myMultis <- myMultis[myMultis != Donor]
  mG <- VPA_Res[[celltype]][[Donor]][[Protocol]]$Results
  mG <- mG[order(mG$log2FoldChange),]
  mG <- rownames(mG)[!is.na(mG$padj) & mG$padj < maxPadj & abs(mG$log2FoldChange) > minFC]
  otherGenes <- c()
  for(i in myMultis){
    tG <- VPA_Res[[celltype]][[i]][[Protocol]]$Results
    tG <- rownames(tG)[tG$padj < maxPadj & !is.na(tG$padj) & abs(tG$log2FoldChange) > minFC]
    otherGenes <- unique(c(otherGenes,tG))
  }
  uniqueGenes <- mG[mG %ni% otherGenes]
  return(uniqueGenes)
}

myUniverse <- rownames(pbulk)[apply(pbulk,1,sum)>100]
getUniqueGO <- function(celltype="CFuPN", Donor="CW", Protocol="MultiDonor", minFC=0.5, maxPadj=0.1, pvalueCutoff=0.1){
  myGenes <- getUniqueGenes3(celltype=celltype, Donor=Donor, Protocol=Protocol, minFC=minFC, maxPadj=maxPadj)
  myGO <- list()
  myGO$BP <- enrichGO(gene = myGenes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = myUniverse, pvalueCutoff=pvalueCutoff)
  myGO$MF <- enrichGO(gene = myGenes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", universe = myUniverse, pvalueCutoff=pvalueCutoff)
  return(myGO)
}

CW <- getUniqueGO(celltype="Global", Donor="CW")
H1 <- getUniqueGO(celltype="Global", Donor="H1")
Mito <- getUniqueGO(celltype="Global", Donor="Mito210")
PGP <- getUniqueGO(celltype="Global", Donor="PGP1")
a11 <- getUniqueGO(celltype="Global", Donor="11a")


BetweenDonors <- list()
myDonors <- c("CW","H1","Mito210","PGP1","11a")
for(i in 1:4){
  for(j in (i+1):5){
    BetweenDonors[[paste0(myDonors[i],"_",myDonors[j])]] <- GetDE_Results(Treatments=c("Control"),
                                                      IsGlobal = TRUE,
                                                      Genotypes = myDonors[c(i,j)],
                                                      Protocols=c("MultiDonor"),
                                                      myEffect = "Genotype")
  }
}
BetweenSingleDonors <- list()
for(i in 1:4){
  for(j in (i+1):5){
    BetweenSingleDonors[[paste0(myDonors[i],"_",myDonors[j])]] <- GetDE_Results(Treatments=c("Control"),
                                                                          IsGlobal = TRUE,
                                                                          Genotypes = myDonors[c(i,j)],
                                                                          Protocols=c("SingleDonor"),
                                                                          myEffect = "Genotype")
  }
}

myDonors <- c("H1","Mito210","PGP1")
BetweenNPCDonors <- list()
for(i in 1:2){
  for(j in (i+1):3){
    BetweenNPCDonors[[paste0(myDonors[i],"_",myDonors[j])]] <- GetDE_Results(Treatments=c("Control"),
                                                                                IsGlobal = TRUE,
                                                                                Genotypes = myDonors[c(i,j)],
                                                                                Protocols=c("NPC"),
                                                                                myEffect = "Genotype")
  }
}


tmp <- GetDE_Results(Treatments=c("Control"),
                     CellTypes = c("CFuPN"),
                     Genotypes = c("PGP1","Mito210"),
                     Protocols=c("MultiDonor"),
                     myEffect = "Genotype")

p1 <- makeVolcano(BetweenDonors[["CW_H1"]]$Results, myTitle="CW vs H1\nMultiDonor NSC; Global") + NoLegend()
p2 <- makeVolcano(BetweenDonors[["CW_Mito210"]]$Results, myTitle="CW vs Mito210\nMultiDonor NSC; Global") + NoLegend()
p3 <- makeVolcano(BetweenDonors[["CW_PGP1"]]$Results, myTitle="CW vs PGP1\nMultiDonor NSC; Global") + NoLegend()
p4 <- makeVolcano(BetweenDonors[["CW_11a"]]$Results, myTitle="CW vs 11a\nMultiDonor NSC; Global") + NoLegend()
p5 <- makeVolcano(BetweenDonors[["H1_Mito210"]]$Results, myTitle="H1 vs Mito210\nMultiDonor NSC; Global") + NoLegend()
p6 <- makeVolcano(BetweenDonors[["H1_PGP1"]]$Results, myTitle="H1 vs PGP1\nMultiDonor NSC; Global") + NoLegend()
p7 <- makeVolcano(BetweenDonors[["H1_11a"]]$Results, myTitle="H1 vs 11a\nMultiDonor NSC; Global") + NoLegend()
p8 <- makeVolcano(BetweenDonors[["Mito210_PGP1"]]$Results, myTitle="Mito210 vs PGP1\nMultiDonor NSC; Global") + NoLegend()
p9 <- makeVolcano(BetweenDonors[["Mito210_11a"]]$Results, myTitle="Mito210 vs 11a\nMultiDonor NSC; Global") + NoLegend()
p10 <- makeVolcano(BetweenDonors[["PGP1_11a"]]$Results, myTitle="PGP1 vs 11a\nMultiDonor NSC; Global") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_MultiDonor_BetweenDonors.svg", height=12, width=12)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_MultiDonor_BetweenDonors.jpeg", height=1000, width=1000)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()

p1 <- makeVolcano(BetweenSingleDonors[["CW_H1"]]$Results, myTitle="CW vs H1\nSingleDonor NSC; Global") + NoLegend()
p2 <- makeVolcano(BetweenSingleDonors[["CW_Mito210"]]$Results, myTitle="CW vs Mito210\nSingleDonor NSC; Global") + NoLegend()
p3 <- makeVolcano(BetweenSingleDonors[["CW_PGP1"]]$Results, myTitle="CW vs PGP1\nSingleDonor NSC; Global") + NoLegend()
p4 <- makeVolcano(BetweenSingleDonors[["CW_11a"]]$Results, myTitle="CW vs 11a\nSingleDonor NSC; Global") + NoLegend()
p5 <- makeVolcano(BetweenSingleDonors[["H1_Mito210"]]$Results, myTitle="H1 vs Mito210\nSingleDonor NSC; Global") + NoLegend()
p6 <- makeVolcano(BetweenSingleDonors[["H1_PGP1"]]$Results, myTitle="H1 vs PGP1\nSingleDonor NSC; Global") + NoLegend()
p7 <- makeVolcano(BetweenSingleDonors[["H1_11a"]]$Results, myTitle="H1 vs 11a\nSingleDonor NSC; Global") + NoLegend()
p8 <- makeVolcano(BetweenSingleDonors[["Mito210_PGP1"]]$Results, myTitle="Mito210 vs PGP1\nSingleDonor NSC; Global") + NoLegend()
p9 <- makeVolcano(BetweenSingleDonors[["Mito210_11a"]]$Results, myTitle="Mito210 vs 11a\nSingleDonor NSC; Global") + NoLegend()
p10 <- makeVolcano(BetweenSingleDonors[["PGP1_11a"]]$Results, myTitle="PGP1 vs 11a\nSingleDonor NSC; Global") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_SingleDonor_BetweenDonors.svg", height=12, width=12)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_SingleDonor_BetweenDonors.jpeg", height=1000, width=1000)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()

BD_SigGenes_MD <- list()
for(i in 1:4){
  for(j in (i+1):5){
    BD_SigGenes_MD[[paste(myDonors[i],myDonors[j],sep="-")]] <- getSigGenes(BetweenDonors[[paste(myDonors[i],myDonors[j],sep="_")]])
    tmp1 <- BetweenDonors[[paste(myDonors[i],myDonors[j],sep="_")]]$SimpUp@result
    tmp2 <- BetweenDonors[[paste(myDonors[i],myDonors[j],sep="_")]]$SimpDown@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
    tmp2$Direction <- rep("Down", nrow(tmp2))
    BD_SigGenes_MD[[paste0(myDonors[i],"-",myDonors[j],"_GO")]] <- rbind(tmp1,tmp2)
  }
}
write.xlsx(BD_SigGenes_MD, "/Users/tfaits/Documents/Chimera/Results/MultiDonor_NSC_DEGs_and_GO_Between_Donors.xlsx", rowNames=TRUE)
BD_SigGenes_SD <- list()
for(i in 1:4){
  for(j in (i+1):5){
    BD_SigGenes_SD[[paste(myDonors[i],myDonors[j],sep="-")]] <- getSigGenes(BetweenSingleDonors[[paste(myDonors[i],myDonors[j],sep="_")]])
    tmp1 <- BetweenSingleDonors[[paste(myDonors[i],myDonors[j],sep="_")]]$SimpUp@result
    tmp2 <- BetweenSingleDonors[[paste(myDonors[i],myDonors[j],sep="_")]]$SimpDown@result
    tmp1$Direction <- rep("Up", nrow(tmp1))
    tmp2$Direction <- rep("Down", nrow(tmp2))
    BD_SigGenes_SD[[paste0(myDonors[i],"-",myDonors[j],"_GO")]] <- rbind(tmp1,tmp2)
  }
}
write.xlsx(BD_SigGenes_SD, "/Users/tfaits/Documents/Chimera/Results/SingleDonor_NSC_DEGs_and_GO_Between_Donors.xlsx", rowNames=TRUE)





makeScatter <- function(inRes1, inRes2, myTitle="No Title", xLab="", yLab=""){
  shGens <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[shGens,]
  inRes2 <- inRes2[shGens,]
  inRes1$SecondLFC <- inRes2$log2FoldChange
  inRes1$SecondPadj <- inRes2$padj
  inRes1 <- inRes1[!is.na(inRes1$padj) & !is.na(inRes1$SecondPadj),]
  myCor <- data.frame(R2=paste("R^2 =",round(cor(inRes1$log2FoldChange, inRes1$SecondLFC), digits=3)), xpos=-Inf, ypos=Inf, hjustvar=-0.5, vjustvar=6)
  myPlot <- ggplot(inRes1, aes(x=log2FoldChange, y=SecondLFC, size=-log10(padj), color=-log10(SecondPadj))) +
    geom_point(alpha=0.6) +
    scale_color_gradientn(colors=c("#6b3b8b","#eddb1b","#fffc8c"), name="-log10(pvalue)\nin SingleDonors", limits=c(0,171), breaks=c(0, 50, 100, 150)) +
    scale_size_continuous(name="-log10(pvalue)\nin MultiDonors", range = c(0.3, 3), limits = c(0,307)) +
    ggtitle(myTitle) +
    xlim(c(-6,6)) + ylim(c(-6,6)) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    geom_text(data=myCor, aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=R2), inherit.aes = FALSE, fontface="bold") +
    xlab(xLab) + ylab(yLab)
  return(myPlot)
}

makeScatter <- function(inRes1, inRes2, myTitle="No Title", xLab="", yLab=""){
  shGens <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[shGens,]
  inRes2 <- inRes2[shGens,]
  inRes1$SecondLFC <- inRes2$log2FoldChange
  inRes1$SecondPadj <- inRes2$padj
  inRes1 <- inRes1[!is.na(inRes1$padj) & !is.na(inRes1$SecondPadj),]
  myCor <- data.frame(R2=paste("R^2 =",round(cor(inRes1$log2FoldChange, inRes1$SecondLFC), digits=3)), xpos=-Inf, ypos=Inf, hjustvar=-0.5, vjustvar=6)
  myPlot <- ggplot(inRes1, aes(x=log2FoldChange, y=SecondLFC, size=-log10(padj), color=-log10(SecondPadj))) +
    geom_point(alpha=0.8) +
    scale_color_gradientn(colors=c("#fffcac","#9b8b8c","#7b5b8c","#6b3b8b","#4b1b8b"), name="-log10(pvalue)\nin SingleDonors", limits=c(0,171), breaks=c(0, 50, 100, 150)) +
    scale_size_continuous(name="-log10(pvalue)\nin MultiDonors", range = c(0.3, 3), limits = c(0,307)) +
    ggtitle(myTitle) +
    xlim(c(-6,6)) + ylim(c(-6,6)) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    geom_text(data=myCor, aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=R2), inherit.aes = FALSE, fontface="bold") +
    xlab(xLab) + ylab(yLab)
  return(myPlot)
}
makeScatter5 <- function(inRes1, inRes2, myTitle="No Title", xLab="", yLab=""){
  shGens <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[shGens,]
  inRes2 <- inRes2[shGens,]
  inRes1$SecondLFC <- inRes2$log2FoldChange
  inRes1$SecondPadj <- inRes2$padj
  inRes1 <- inRes1[!is.na(inRes1$padj) & !is.na(inRes1$SecondPadj),]
  myCor <- data.frame(R2=paste("R^2 =",round(cor(inRes1$log2FoldChange, inRes1$SecondLFC), digits=3)), xpos=-Inf, ypos=Inf, hjustvar=-0.5, vjustvar=6)
  myPlot <- ggplot(inRes1, aes(x=log2FoldChange, y=SecondLFC, size=-log10(padj), color=-log10(SecondPadj))) +
    geom_point(alpha=1) +
    scale_color_gradient2(low="grey85", mid="red", high="red", midpoint=100, name="-log10(pvalue)\nin SingleDonors", limits=c(0,171), breaks=c(0, 50, 100, 150)) +
    scale_size_continuous(name="-log10(pvalue)\nin MultiDonors", range = c(0.3, 3), limits = c(0,307)) +
    ggtitle(myTitle) +
    xlim(c(-6,6)) + ylim(c(-6,6)) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    geom_text(data=myCor, aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=R2), inherit.aes = FALSE, fontface="bold") +
    xlab(xLab) + ylab(yLab)
  return(myPlot)
}

makeScatter6 <- function(inRes1, inRes2, myTitle="No Title", xLab="", yLab=""){
  shGens <- rownames(inRes1)[rownames(inRes1) %in% rownames(inRes2)]
  inRes1 <- inRes1[shGens,]
  inRes2 <- inRes2[shGens,]
  inRes1$SecondLFC <- inRes2$log2FoldChange
  inRes1$SecondPadj <- log10(inRes2$baseMean)
  inRes1 <- inRes1[!is.na(inRes1$padj) & !is.na(inRes1$SecondPadj),]
  myCor <- data.frame(R2=paste("R^2 =",round(cor(inRes1$log2FoldChange, inRes1$SecondLFC), digits=3)), xpos=-Inf, ypos=Inf, hjustvar=-0.5, vjustvar=6)
  myPlot <- ggplot(inRes1, aes(x=log2FoldChange, y=SecondLFC, size=log10(baseMean), color=SecondPadj)) +
    geom_point(alpha=0.6) +
    #scale_color_gradient2(low="grey85", mid="red", high="red", midpoint=100, name="-log10(pvalue)\nin SingleDonors", limits=c(0,171), breaks=c(0, 50, 100, 150)) +
    #scale_size_continuous(name="-log10(pvalue)\nin MultiDonors", range = c(0.3, 3), limits = c(0,307)) +
    ggtitle(myTitle) +
    xlim(c(-6,6)) + ylim(c(-6,6)) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    geom_text(data=myCor, aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=R2), inherit.aes = FALSE, fontface="bold") +
    xlab(xLab) + ylab(yLab)
  return(myPlot)
}

p1 <- makeScatter(BetweenDonors[["CW_H1"]]$Results, BetweenSingleDonors[["CW_H1"]]$Results, myTitle="CW vs H1", xLab="log2(FoldChange) in MultiDonor", yLab="log2(FoldChange) in SingleDonor") + NoLegend()
p2 <- makeScatter(BetweenDonors[["CW_Mito210"]]$Results, BetweenSingleDonors[["CW_Mito210"]]$Results,  myTitle="CW vs Mito210")+ NoLegend()
p3 <- makeScatter(BetweenDonors[["CW_PGP1"]]$Results, BetweenSingleDonors[["CW_PGP1"]]$Results,  myTitle="CW vs PGP1")+ NoLegend()
p4 <- makeScatter(BetweenDonors[["CW_11a"]]$Results, BetweenSingleDonors[["CW_11a"]]$Results,  myTitle="CW vs 11a")+ NoLegend()
p5 <- makeScatter(BetweenDonors[["H1_Mito210"]]$Results, BetweenSingleDonors[["H1_Mito210"]]$Results,  myTitle="H1 vs Mito210")+ NoLegend()
p6 <- makeScatter(BetweenDonors[["H1_PGP1"]]$Results, BetweenSingleDonors[["H1_PGP1"]]$Results,  myTitle="H1 vs PGP1")+ NoLegend()
p7 <- makeScatter(BetweenDonors[["H1_11a"]]$Results, BetweenSingleDonors[["H1_11a"]]$Results,  myTitle="H1 vs 11a")+ NoLegend()
p8 <- makeScatter(BetweenDonors[["Mito210_PGP1"]]$Results, BetweenSingleDonors[["Mito210_PGP1"]]$Results,  myTitle="Mito210 vs PGP1")+ NoLegend()
p9 <- makeScatter(BetweenDonors[["Mito210_11a"]]$Results, BetweenSingleDonors[["Mito210_11a"]]$Results,  myTitle="Mito210 vs 11a")+ NoLegend()
p10 <- makeScatter(BetweenDonors[["PGP1_11a"]]$Results, BetweenSingleDonors[["PGP1_11a"]]$Results,  myTitle="PGP1 vs 11a")+ NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Correlation_Of_BetweenDonors_SingleDonor_To_MultiDonor_Feb7_RevCol.svg", height=14, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Correlation_Of_BetweenDonors_SingleDonor_To_MultiDonor_Feb7_RevCol.jpeg", res=700, height=14, width=14, units = "in")
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()

p1 <- makeScatter(BetweenDonors[["CW_H1"]]$Results, BetweenSingleDonors[["CW_H1"]]$Results, myTitle="CW vs H1", xLab="log2(FoldChange) in MultiDonor", yLab="log2(FoldChange) in SingleDonor")
legend <- cowplot::get_legend(p1)
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Corr_Legend_RevCol.svg", height=4, width=2)
grid.newpage()
grid.draw(legend)
dev.off()

p1 <- makeScatter5(BetweenDonors[["CW_H1"]]$Results, BetweenSingleDonors[["CW_H1"]]$Results, myTitle="CW vs H1", xLab="log2(FoldChange) in MultiDonor", yLab="log2(FoldChange) in SingleDonor") + NoLegend()
p2 <- makeScatter5(BetweenDonors[["CW_Mito210"]]$Results, BetweenSingleDonors[["CW_Mito210"]]$Results,  myTitle="CW vs Mito210")+ NoLegend()
p3 <- makeScatter5(BetweenDonors[["CW_PGP1"]]$Results, BetweenSingleDonors[["CW_PGP1"]]$Results,  myTitle="CW vs PGP1")+ NoLegend()
p4 <- makeScatter5(BetweenDonors[["CW_11a"]]$Results, BetweenSingleDonors[["CW_11a"]]$Results,  myTitle="CW vs 11a")+ NoLegend()
p5 <- makeScatter5(BetweenDonors[["H1_Mito210"]]$Results, BetweenSingleDonors[["H1_Mito210"]]$Results,  myTitle="H1 vs Mito210")+ NoLegend()
p6 <- makeScatter5(BetweenDonors[["H1_PGP1"]]$Results, BetweenSingleDonors[["H1_PGP1"]]$Results,  myTitle="H1 vs PGP1")+ NoLegend()
p7 <- makeScatter5(BetweenDonors[["H1_11a"]]$Results, BetweenSingleDonors[["H1_11a"]]$Results,  myTitle="H1 vs 11a")+ NoLegend()
p8 <- makeScatter5(BetweenDonors[["Mito210_PGP1"]]$Results, BetweenSingleDonors[["Mito210_PGP1"]]$Results,  myTitle="Mito210 vs PGP1")+ NoLegend()
p9 <- makeScatter5(BetweenDonors[["Mito210_11a"]]$Results, BetweenSingleDonors[["Mito210_11a"]]$Results,  myTitle="Mito210 vs 11a")+ NoLegend()
p10 <- makeScatter5(BetweenDonors[["PGP1_11a"]]$Results, BetweenSingleDonors[["PGP1_11a"]]$Results,  myTitle="PGP1 vs 11a")+ NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Correlation_Of_BetweenDonors_SingleDonor_To_MultiDonor_Feb7_NewColor.svg", height=14, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Correlation_Of_BetweenDonors_SingleDonor_To_MultiDonor_Feb5.jpeg", height=1000, width=1000)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, layout_matrix=rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9),c(NA,NA,NA,10)))
dev.off()


makeScatter(BetweenDonors[["CW_H1"]]$Results, BetweenSingleDonors[["CW_H1"]]$Results, myTitle="CW vs H1", xLab="log2(FoldChange) in MultiDonor", yLab="log2(FoldChange) in SingleDonor")
ggsave("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Cor_Legend_Plot.svg", height=4, width=6)


numCW <- length(getUniqueGenes3(celltype="CFuPN",Donor = "CW", Protocol = "MultiDonor"))
numH1 <- length(getUniqueGenes3(celltype="CFuPN",Donor = "H1", Protocol = "MultiDonor"))
numMito <- length(getUniqueGenes3(celltype="CFuPN",Donor = "Mito210", Protocol = "MultiDonor"))
numPGP <- length(getUniqueGenes3(celltype="CFuPN",Donor = "PGP1", Protocol = "MultiDonor"))
num11a <- length(getUniqueGenes3(celltype="CFuPN",Donor = "11a", Protocol = "MultiDonor"))


UniqueGenes <- c(getUniqueGenes3(celltype="CFuPN",Donor = "CW", Protocol = "MultiDonor")[c(1:20, (numCW-19):numCW)],
                 getUniqueGenes3(celltype="CFuPN",Donor = "H1", Protocol = "MultiDonor")[c(1:20, (numH1-19):numH1)],
                 getUniqueGenes3(celltype="CFuPN",Donor = "Mito210", Protocol = "MultiDonor")[c(1:20, (numMito-19):numMito)],
                 getUniqueGenes3(celltype="CFuPN",Donor = "PGP1", Protocol = "MultiDonor")[c(1:20, (numPGP-19):numPGP)],
                 getUniqueGenes3(celltype="CFuPN",Donor = "11a", Protocol = "MultiDonor")[c(1:20, (num11a-19):num11a)])
UniqueGenesAll <- c(getUniqueGenes3(celltype="CFuPN",Donor = "CW", Protocol = "MultiDonor"),getUniqueGenes3(celltype="CFuPN",Donor = "H1", Protocol = "MultiDonor"), getUniqueGenes3(celltype="CFuPN",Donor = "Mito210", Protocol = "MultiDonor"), getUniqueGenes3(celltype="CFuPN",Donor = "PGP1", Protocol = "MultiDonor"), getUniqueGenes3(celltype="CFuPN",Donor = "11a", Protocol = "MultiDonor"))


myLFC <- data.frame(CW=VPA_Res$CFuPN$CW$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    H1=VPA_Res$CFuPN$H1$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    Mito210=VPA_Res$CFuPN$Mito210$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    PGP1=VPA_Res$CFuPN$PGP1$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    a11=VPA_Res$CFuPN[["11a"]]$MultiDonor$Results[UniqueGenes,"log2FoldChange"])
myLFC2 <- data.frame(CW=VPA_Res$CFuPN$CW$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                     CW_SD=VPA_Res$CFuPN$CW$SingleDonor$Results[UniqueGenes,"log2FoldChange"],
                     H1=VPA_Res$CFuPN$H1$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    H1_SD=VPA_Res$CFuPN$H1$SingleDonor$Results[UniqueGenes,"log2FoldChange"],
                    Mito210=VPA_Res$CFuPN$Mito210$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    Mito210_SD=VPA_Res$CFuPN$Mito210$SingleDonor$Results[UniqueGenes,"log2FoldChange"],
                    PGP1=VPA_Res$CFuPN$PGP1$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    PGP1_SD=VPA_Res$CFuPN$PGP1$SingleDonor$Results[UniqueGenes,"log2FoldChange"],
                    a11=VPA_Res$CFuPN[["11a"]]$MultiDonor$Results[UniqueGenes,"log2FoldChange"],
                    a11_SD=VPA_Res$CFuPN[["11a"]]$SingleDonor$Results[UniqueGenes,"log2FoldChange"])
myLFC3 <- data.frame(CW=VPA_Res$CFuPN$CW$MultiDonor$Results[UniqueGenesAll,"log2FoldChange"],
                     H1=VPA_Res$CFuPN$H1$MultiDonor$Results[UniqueGenesAll,"log2FoldChange"],
                     Mito210=VPA_Res$CFuPN$Mito210$MultiDonor$Results[UniqueGenesAll,"log2FoldChange"],
                     PGP1=VPA_Res$CFuPN$PGP1$MultiDonor$Results[UniqueGenesAll,"log2FoldChange"],
                     a11=VPA_Res$CFuPN[["11a"]]$MultiDonor$Results[UniqueGenesAll,"log2FoldChange"])
myLFC3[is.na(myLFC3)] <- 0

mySideAnno <- rowAnnotation(Genes_From_Donor=c(rep("CW",40),rep("H1",40),rep("Mito210",40),rep("PGP1",40),rep("11a",40)),
                            col=list(Genes_From_Donor=c("CW"="#0a9396", "GM"="#94d2bd", "H1"="#E9d8a6", "PGP1"="#ee9b00", "Mito210"="#ca6702","11a"="darkred")))

mySideAnno2 <- rowAnnotation(Donor=c(rep("CW",50),rep("H1",50),rep("Mito210",50),rep("PGP1",50),rep("11a",50)),
                            col=list(Donor=c("CW"="#0a9396", "GM"="#94d2bd", "H1"="#E9d8a6", "PGP1"="#ee9b00", "Mito210"="#ca6702","11a"="darkred")))
mySideAnno3 <- rowAnnotation(Donor=c(rep("CW",length(getUniqueGenes3(celltype="CFuPN",Donor = "CW", Protocol = "MultiDonor"))),
                                     rep("H1",length(getUniqueGenes3(celltype="CFuPN",Donor = "H1", Protocol = "MultiDonor"))),
                                     rep("Mito210",length(getUniqueGenes3(celltype="CFuPN",Donor = "Mito210", Protocol = "MultiDonor"))),
                                     rep("PGP1",length(getUniqueGenes3(celltype="CFuPN",Donor = "PGP1", Protocol = "MultiDonor"))),
                                     rep("11a",length(getUniqueGenes3(celltype="CFuPN",Donor = "11a", Protocol = "MultiDonor")))),
                             col=list(Donor=c("CW"="#0a9396", "GM"="#94d2bd", "H1"="#E9d8a6", "PGP1"="#ee9b00", "Mito210"="#ca6702","11a"="darkred")))

myLFC[is.na(myLFC)] <- 0
myLFC2[is.na(myLFC2)] <- 0
p1 <- Heatmap(myLFC, column_order = 1:5, left_annotation = mySideAnno,
              row_order=1:200, name="Log2(Fold Change)\nIn VPA vs Control",
              column_title = "Genes Uniquely DE in Each Donor\nVPA vs Control")
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Unique_Genes_For_Each_Donor_log2FC_in_VPA_Heatmap.svg", height=8, width=5)
grid.newpage()
p1
dev.off()

p2 <- Heatmap(myLFC2, column_order = 1:10, left_annotation = mySideAnno,
              row_order=1:200, name="Log2(Fold Change)\nIn VPA vs Control",
              column_title = "Genes Uniquely DE in Each Donor\nVPA vs Control")
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Unique_Genes_For_Each_Donor_log2FC_in_VPA_Heatmap_MD_AND_SD.svg", height=8, width=5)
grid.newpage()
p2
dev.off()

gseaplot()

Heatmap(abs(myLFC2), column_order = 1:5, left_annotation = mySideAnno2, row_order= 1:250)
Heatmap(t(scale(t(myLFC3))), column_order = 1:5, left_annotation = mySideAnno3, row_order=1:nrow(myLFC3))

p1 <- makeVolcano2(SingleMulti$aRG_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - aRGs") + NoLegend()
p2 <- makeVolcano2(SingleMulti$oRG_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - oRGs") + NoLegend()
p3 <- makeVolcano2(SingleMulti$IP_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - IPs") + NoLegend()
p4 <- makeVolcano2(SingleMulti$PN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - PNs") + NoLegend()
p5 <- makeVolcano2(SingleMulti$CFuPN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - CFuPNs") + NoLegend()
p6 <- makeVolcano2(SingleMulti$CPN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - CPNs") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_AllDonors_EachCellType_AdjustedP_Feb5.svg", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_AllDonors_EachCellType_AdjustedP_Feb5.jpeg", res=500, height=8, width=14, units="in")
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()


makeVolcano <- function(inTab, myTitle="No Title"){
  inTab$Sig <- "Not Sig"
  inTab$Sig[inTab$log2FoldChange > 0.5 & inTab$padj < 0.05] <- "Up in SD"
  inTab$Sig[inTab$log2FoldChange < -0.5 & inTab$padj < 0.05] <- "Down in SD"
  myLabs <- c(paste0("\n% Down: ",round(sum(inTab$Sig=="Down in SD")/nrow(inTab)*100,digits=1)),
              paste0("% Up: ",round(sum(inTab$Sig=="Up in SD")/nrow(inTab)*100,digits=1)),
              paste0("\n%Not Sig: ",round(sum(inTab$Sig == "Not Sig")/nrow(inTab)*100,digits=1)))
  ggplot(inTab, aes(x=log2FoldChange, y=-log10(pvalue), color=Sig)) +
    geom_point(size=1, alpha=0.7) +
    scale_color_manual(values=SigColors[unique(inTab$Sig)]) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    coord_cartesian(clip = 'off') +
    geom_abline(slope=0, intercept=min(-log10(inTab$pvalue[inTab$padj < 0.05 & !is.na(inTab$padj)])), linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_vline(xintercept=0.5, linetype="dashed") +
    annotate("text",  x=c(-Inf,Inf,0), y = Inf, label = myLabs, hjust=c(0,1,0.5), vjust=0) +
    ggtitle(myTitle)
}
makeVolcano(SingleMulti$aRG_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - aRGs") + NoLegend()
p1 <- makeVolcano(SingleMulti$aRG_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - aRGs") + NoLegend()
p2 <- makeVolcano(SingleMulti$oRG_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - oRGs") + NoLegend()
p3 <- makeVolcano(SingleMulti$IP_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - IPs") + NoLegend()
p4 <- makeVolcano(SingleMulti$PN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - PNs") + NoLegend()
p5 <- makeVolcano(SingleMulti$CFuPN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - CFuPNs") + NoLegend()
p6 <- makeVolcano(SingleMulti$CPN_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - CPNs") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_AllDonors_EachCellType_AdjustedP_Feb8.svg", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()
png("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_AllDonors_EachCellType_Feb14.png", height=8, width=14, units="in", res=300)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()


p1 <- makeVolcano2(SingleMulti$Global_SingleDonor_CW$Results, "Multi-Donor vs Single-Donor\nCW - Global") + NoLegend()
p2 <- makeVolcano2(SingleMulti$Global_SingleDonor_H1$Results, "Multi-Donor vs Single-Donor\nH1 - Global") + NoLegend()
p3 <- makeVolcano2(SingleMulti$Global_SingleDonor_Mito210$Results, "Multi-Donor vs Single-Donor\nMito210 - Global") + NoLegend()
p4 <- makeVolcano2(SingleMulti$Global_SingleDonor_PGP1$Results, "Multi-Donor vs Single-Donor\nPGP1 - Global") + NoLegend()
p5 <- makeVolcano2(SingleMulti$Global_SingleDonor_11a$Results, "Multi-Donor vs Single-Donor\n11a - Global") + NoLegend()
p6 <- makeVolcano2(SingleMulti$Global_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - Global") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_EachDonor_Global_AdjustedP.svg", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()
jpeg("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_EachDonor_Global_AdjustedP_Feb5.jpeg", res=500, units="in", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()

p1 <- makeVolcano(SingleMulti$Global_SingleDonor_CW$Results, "Multi-Donor vs Single-Donor\nCW - Global") + NoLegend()
p2 <- makeVolcano(SingleMulti$Global_SingleDonor_H1$Results, "Multi-Donor vs Single-Donor\nH1 - Global") + NoLegend()
p3 <- makeVolcano(SingleMulti$Global_SingleDonor_Mito210$Results, "Multi-Donor vs Single-Donor\nMito210 - Global") + NoLegend()
p4 <- makeVolcano(SingleMulti$Global_SingleDonor_PGP1$Results, "Multi-Donor vs Single-Donor\nPGP1 - Global") + NoLegend()
p5 <- makeVolcano(SingleMulti$Global_SingleDonor_11a$Results, "Multi-Donor vs Single-Donor\n11a - Global") + NoLegend()
p6 <- makeVolcano(SingleMulti$Global_SingleDonor_Global$Results, "Multi-Donor vs Single-Donor\nAll Donors - Global") + NoLegend()
svglite("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_EachDonor_Global_Feb8.svg", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()
png("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Volcano_EachDonor_Global_Feb14.png", res=300, units="in", height=8, width=14)
grid.newpage()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
dev.off()


makeVolcano2(BetweenDonors$)


CW <- getUniqueGO(celltype="CPN", Donor="CW")
H1 <- getUniqueGO(celltype="CPN", Donor="H1")
Mito <- getUniqueGO(celltype="CPN", Donor="Mito210")
PGP <- getUniqueGO(celltype="CPN", Donor="PGP1")
a11 <- getUniqueGO(celltype="CPN", Donor="11a")

makeDotPlot <- function(inGO){
  tmp <- simplify(inGO, cutoff=0.6)
  myP <- dotplot(tmp)
  return(myP)
}
dotplot(CW$BP)
dotplot(H1$BP)
makeDotPlot(Mito$BP)
makeDotPlot(PGP$BP)
makeDotPlot(a11$BP)

length(CW$BP@gene)
length(H1$BP@gene)
length(Mito$BP@gene)
length(PGP$BP@gene)
length(a11$BP@gene)


CW2 <- getUniqueGO(celltype="Global", Donor="CW", Protocol = "SingleDonor")
H12 <- getUniqueGO(celltype="Global", Donor="H1", Protocol = "SingleDonor")
Mito2 <- getUniqueGO(celltype="Global", Donor="Mito210", Protocol = "SingleDonor")
PGP2 <- getUniqueGO(celltype="Global", Donor="PGP1", Protocol = "SingleDonor")
a112 <- getUniqueGO(celltype="Global", Donor="11a", Protocol = "SingleDonor")
makeDotPlot(CW2$BP)
makeDotPlot(H12$BP)
makeDotPlot(Mito2$BP)
makeDotPlot(PGP2$BP)
makeDotPlot(a112$BP)

dotplot(CW$BP)
dotplot(H1$BP)
makeDotPlot(Mito$BP)
makeDotPlot(PGP$BP)
makeDotPlot(a11$BP)

length(CW2$BP@gene)
length(H12$BP@gene)
length(Mito2$BP@gene)
length(PGP2$BP@gene)
length(a112$BP@gene)

plotKS2 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE, myTitle=""){
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  queryRes <- queryRes[order((queryRes[,rankMetric]), decreasing=Decreasing),]
  queryRanks <- which(rownames(queryRes) %in% myGeneList)
  myP <- hist(queryRanks, breaks=nrow(queryRes), main=myTitle, xlab="", ylab="")
  return(myP)
}

runKS_Test2 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE){
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  queryRanks <- which(rownames(queryRes) %in% myGeneList)
  return(ks.test(queryRes[queryRanks,rankMetric], queryRes[,rankMetric]))
}

runKS_Test2(queryRes = VPA_Res$Global$CW$SingleDonor$Results, myGeneList=CW$Global$BP@gene, rankMetric="pvalue", Decreasing=FALSE)
runKS_Test2(queryRes = VPA_Res$Global$H1$SingleDonor$Results, myGeneList=H1$Global$BP@gene, rankMetric="pvalue", Decreasing=FALSE)
runKS_Test2(queryRes = VPA_Res$Global$Mito210$SingleDonor$Results, myGeneList=Mito$Global$BP@gene, rankMetric="pvalue", Decreasing=FALSE)
runKS_Test2(queryRes = VPA_Res$Global$PGP1$SingleDonor$Results, myGeneList=PGP$Global$BP@gene, rankMetric="pvalue", Decreasing=FALSE)
runKS_Test2(queryRes = VPA_Res$Global[["11a"]]$SingleDonor$Results, myGeneList=a11$Global$BP@gene, rankMetric="pvalue", Decreasing=FALSE)

par(mfrow=c(5,1))
plotKS2(queryRes = VPA_Res$Global$CW$SingleDonor$Results, myGeneList=CW$Global$BP@gene, myTitle="CW: KS-test pval < 1e-16", rankMetric="log2FoldChange", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$Global$H1$SingleDonor$Results, myGeneList=H1$Global$BP@gene, myTitle="H1: KS-test pval = 6.27e-8", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$Global$Mito210$SingleDonor$Results, myGeneList=Mito$Global$BP@gene, myTitle="Mito210: KS-test pval = 4.78e-10", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$Global$PGP1$SingleDonor$Results, myGeneList=PGP$Global$BP@gene, myTitle="PGP1: KS-test pval = 1.04e-13", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$Global[["11a"]]$SingleDonor$Results, myGeneList=a11$Global$BP@gene, myTitle="11a: KS-test pval < 1e-16", rankMetric="pvalue", Decreasing=FALSE)


runKS_Test2(queryRes = VPA_Res$CFuPN$CW$SingleDonor$Results, myGeneList=CW$CFuPN$BP@gene, rankMetric="pvalue", Decreasing=FALSE)$p.value
runKS_Test2(queryRes = VPA_Res$CFuPN$H1$SingleDonor$Results, myGeneList=H1$CFuPN$BP@gene, rankMetric="pvalue", Decreasing=FALSE)$p.value
runKS_Test2(queryRes = VPA_Res$CFuPN$Mito210$SingleDonor$Results, myGeneList=Mito$CFuPN$BP@gene, rankMetric="pvalue", Decreasing=FALSE)$p.value
runKS_Test2(queryRes = VPA_Res$CFuPN$PGP1$SingleDonor$Results, myGeneList=PGP$CFuPN$BP@gene, rankMetric="pvalue", Decreasing=FALSE)$p.value
runKS_Test2(queryRes = VPA_Res$CFuPN[["11a"]]$SingleDonor$Results, myGeneList=a11$CFuPN$BP@gene, rankMetric="pvalue", Decreasing=FALSE)$p.value

par(mfrow=c(5,1))
plotKS2(queryRes = VPA_Res$CFuPN$CW$SingleDonor$Results, myGeneList=CW$CFuPN$BP@gene, myTitle="CW: KS-test pval = 3.54e-8", rankMetric="log2FoldChange", Decreasing=TRUE)
plotKS2(queryRes = VPA_Res$CFuPN$H1$SingleDonor$Results, myGeneList=H1$CFuPN$BP@gene, myTitle="H1: KS-test pval = 5.21e-3", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$CFuPN$Mito210$SingleDonor$Results, myGeneList=Mito$CFuPN$BP@gene, myTitle="Mito210: KS-test pval = 3.26e-4", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$CFuPN$PGP1$SingleDonor$Results, myGeneList=PGP$CFuPN$BP@gene, myTitle="PGP1: KS-test pval = 5.25e-6", rankMetric="pvalue", Decreasing=FALSE)
plotKS2(queryRes = VPA_Res$CFuPN[["11a"]]$SingleDonor$Results, myGeneList=a11$CFuPN$BP@gene, myTitle="11a: KS-test pval < 1e-16", rankMetric="pvalue", Decreasing=FALSE)



plotKS3 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE, myTitle=""){
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  qGens <- rownames(queryRes) %in% myGeneList
  dVal <- c(0)
  isIn <- sum(qGens)
  notIn <- length(qGens)-sum(qGens)
  UpHeight <- 1/isIn
  DownHeight <- 1/notIn
  for(i in 1:length(qGens)){
    if(qGens[i]){
      dVal[i+1] <- dVal[i] + UpHeight
    }else{
      dVal[i+1] <- dVal[i] - DownHeight
    }
  }
  myP <- plot(1:length(dVal), dVal, main=myTitle, xlab="", ylab="", type="l", col="green", lwd=3, xaxt="n")
  return(myP)
}
par(mfrow=c(5,1))
plotKS3(queryRes = VPA_Res$Global$CW$SingleDonor$Results, myGeneList=CW$Global$BP@gene, myTitle="CW: KS-test pval = 3.54e-8", rankMetric="pvalue", Decreasing=FALSE)
plotKS3(queryRes = VPA_Res$Global$H1$SingleDonor$Results, myGeneList=H1$Global$BP@gene, myTitle="H1: KS-test pval = 5.21e-3", rankMetric="pvalue", Decreasing=FALSE)
plotKS3(queryRes = VPA_Res$Global$Mito210$SingleDonor$Results, myGeneList=Mito$Global$BP@gene, myTitle="Mito210: KS-test pval = 3.26e-4", rankMetric="pvalue", Decreasing=FALSE)
plotKS3(queryRes = VPA_Res$Global$PGP1$SingleDonor$Results, myGeneList=PGP$Global$BP@gene, myTitle="PGP1: KS-test pval = 5.25e-6", rankMetric="pvalue", Decreasing=FALSE)
plotKS3(queryRes = VPA_Res$Global[["11a"]]$SingleDonor$Results, myGeneList=a11$Global$BP@gene, myTitle="11a: KS-test pval < 1e-16", rankMetric="pvalue", Decreasing=FALSE)



for(i in 1:length(qGens)){
  if(qGens[i]){
    dVal[i+1] <- dVal[i] + UpHeight
  }else{
    dVal[i+1] <- dVal[i]
  }
}

CW <- list()
H1 <- list()
Mito <- list()
PGP <- list()
a11 <- list()
CWs <- list()
H1s <- list()
Mitos <- list()
PGPs <- list()
a11s <- list()
for(celltype in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  print(paste("Working on celltype:",celltype))
  CW[[celltype]] <- getUniqueGO(celltype=celltype, Donor="CW", Protocol = "MultiDonor", pvalueCutoff=0.1)
  H1[[celltype]] <- getUniqueGO(celltype=celltype, Donor="H1", Protocol = "MultiDonor", pvalueCutoff=0.1)
  Mito[[celltype]] <- getUniqueGO(celltype=celltype, Donor="Mito210", Protocol = "MultiDonor", pvalueCutoff=0.1)
  PGP[[celltype]] <- getUniqueGO(celltype=celltype, Donor="PGP1", Protocol = "MultiDonor", pvalueCutoff=0.1)
  a11[[celltype]] <- getUniqueGO(celltype=celltype, Donor="11a", Protocol = "MultiDonor", pvalueCutoff=0.1)
  # print("Starting the Single-Donors...")
  # CWs[[celltype]] <- getUniqueGO(celltype=celltype, Donor="CW", Protocol = "SingleDonor")
  # H1s[[celltype]] <- getUniqueGO(celltype=celltype, Donor="H1", Protocol = "SingleDonor")
  # Mitos[[celltype]] <- getUniqueGO(celltype=celltype, Donor="Mito210", Protocol = "SingleDonor")
  # PGPs[[celltype]] <- getUniqueGO(celltype=celltype, Donor="PGP1", Protocol = "SingleDonor")
  # a11s[[celltype]] <- getUniqueGO(celltype=celltype, Donor="11a", Protocol = "SingleDonor")
}


CountGOs <- function(inGO, ont="BP"){
  tmp <- inGO[[ont]]@result
  return(sum(tmp$qvalue < 0.1))
}

CountSimpleGOs <- function(inGO, ont="BP", cutoff=0.4){
  tmp <- inGO[[ont]]
  tmp <- simplify(tmp, cutoff=cutoff)
  return(sum(tmp@result$qvalue < 0.1))
}



myThing <- list(CW=CW, H1=H1, Mito=Mito, PGP=PGP, a11=a11)

for(j in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  print(j)
  for(i in names(myThing)){
    print(paste(i, CountGOs(myThing[[i]][[j]], ont="BP"),CountGOs(myThing[[i]][[j]], ont="MF")))
  }
  print("")
}


for(j in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  print(j)
  for(i in names(myThing)){
    print(paste(i, CountSimpleGOs(myThing[[i]][[j]], ont="BP"),CountSimpleGOs(myThing[[i]][[j]], ont="MF")))
  }
  print("")
}



for(j in c("Global","aRG","oRG","IP","CFuPN","CPN")){
  print(j)
  for(i in names(myThing)){
    print(paste(i, length(myThing[[i]][[j]]$BP@gene)))
  }
  print("")
}


plotKS4 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE, myTitle="", myColor="green"){
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  tmp <- runKS_Test2(queryRes=queryRes, myGeneList = myGeneList, rankMetric=rankMetric, Decreasing=Decreasing)
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  qGens <- rownames(queryRes) %in% myGeneList
  dVal <- c(0)
  isIn <- sum(qGens)
  notIn <- length(qGens)-sum(qGens)
  UpHeight <- 1/isIn
  DownHeight <- 1/notIn
  for(i in 1:length(qGens)){
    if(qGens[i]){
      dVal[i+1] <- dVal[i] + UpHeight
    }else{
      dVal[i+1] <- dVal[i] - DownHeight
    }
  }
  queryRes$Score <- dVal[2:length(dVal)]
  queryRes$Rank <- 1:nrow(queryRes)
  queryRes$SigGene <- queryRes$Rank
  queryRes$SigGene[rownames(queryRes) %ni% myGeneList] <- NA
  tmpRes <- queryRes[rownames(queryRes) %in% myGeneList,]
  maxY <- max(abs(dVal))/3
  myP <- ggplot(queryRes, aes(x=Rank)) +
    geom_segment(data=tmpRes, aes(x=Rank, xend=Rank, y=0, yend = maxY), size=0.15) +
    geom_line(data=queryRes,aes(y = Score), color=myColor,size=1) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5, size=8)) +
    ggtitle(myTitle) +
    ylab("Enrichment Score") +
    xlab(paste0("Rank of ", rankMetric)) +
    geom_text(aes(x=nrow(queryRes)*0.8, y=maxY*7/3, hjust=0.5, label=paste0("KS Test\np = ", format(tmp$p.value, scientific=TRUE, digits=3), "\nD = ", round(tmp$statistic, digits=4))),
              size=2, inherit.aes = FALSE, fontface="bold")
  return(myP)
}
plotKS6 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE, myTitle="", myColor="green"){
  ##Order the query by the <rankMetric>, then get the ranks of the <numGenes>:
  tmp <- runKS_Test2(queryRes=queryRes, myGeneList = myGeneList, rankMetric=rankMetric, Decreasing=Decreasing)
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  qGens <- rownames(queryRes) %in% myGeneList
  dVal <- c(0)
  isIn <- sum(qGens)
  notIn <- length(qGens)-sum(qGens)
  UpHeight <- 1/isIn
  DownHeight <- 1/notIn
  for(i in 1:length(qGens)){
    if(qGens[i]){
      dVal[i+1] <- dVal[i] + UpHeight
    }else{
      dVal[i+1] <- dVal[i] - DownHeight
    }
  }
  queryRes$Score <- dVal[2:length(dVal)]
  queryRes$Rank <- 1:nrow(queryRes)
  queryRes$SigGene <- queryRes$Rank
  queryRes$SigGene[rownames(queryRes) %ni% myGeneList] <- NA
  tmpRes <- queryRes[rownames(queryRes) %in% myGeneList,]
  maxY <- max(abs(dVal))/3
  myP <- ggplot(queryRes, aes(x=Rank)) +
    geom_segment(data=tmpRes, aes(x=Rank, xend=Rank, y=-maxY, yend = abs(log2FoldChange)/max(abs(log2FoldChange))*maxY*1.5), size=0.2) +
    geom_line(data=queryRes,aes(y = Score), color=myColor,size=1) +
    theme(panel.background = element_blank(),
          plot.title=element_text(face="bold", hjust=0.5, size=8)) +
    ggtitle(myTitle) +
    ylab("Enrichment Score") +
    xlab(paste0("Rank of ", rankMetric)) +
    geom_text(aes(x=nrow(queryRes)*0.8, y=maxY*7/3, hjust=0.5, label=paste0("KS Test\np = ", format(tmp$p.value, scientific=TRUE, digits=3), "\nD = ", round(tmp$statistic, digits=4))),
              size=2, inherit.aes = FALSE, fontface="bold")
  return(myP)
}
plotKS5 <- function(queryRes, myGeneList, rankMetric="log2FoldChange", Decreasing=TRUE, myTitle="", myColor="green"){
  queryRes <- queryRes[order(abs(queryRes[,rankMetric]), decreasing=Decreasing),]
  qGens <- rownames(queryRes) %in% myGeneList
  dVal <- c(0)
  isIn <- sum(qGens)
  notIn <- length(qGens)-sum(qGens)
  UpHeight <- 1/isIn
  DownHeight <- 1/notIn
  for(i in 1:length(qGens)){
    if(qGens[i]){
      dVal[i+1] <- dVal[i] + UpHeight
    }else{
      dVal[i+1] <- dVal[i] - DownHeight
    }
  }
  queryRes$Score <- dVal[2:length(dVal)]
  queryRes$Rank <- 1:nrow(queryRes)
  queryRes$SigGene <- queryRes$Rank
  queryRes$SigGene[rownames(queryRes) %ni% myGeneList] <- NA
  tmpRes <- queryRes[rownames(queryRes) %in% myGeneList,]
  maxY <- max(abs(dVal))/3
  myP <- ggplot(queryRes, aes(x=Rank)) +
    geom_segment(data=tmpRes, aes(x=Rank, xend=Rank, y=0, yend = maxY), size=0.15) +
    geom_line(data=queryRes,aes(y = Score), color=myColor,size=1) +
    theme(panel.background = element_blank(),
          plot.title=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()) +
    ggtitle(myTitle) +
    ylab("") +
    xlab("")
  return(myP)
}

## Make a bunch of KS plots:
for(i in myDonors){
  print(paste("Starting",i))
  plotKS4(queryRes = VPA_Res$Global[[i]]$SingleDonor$Results,
          myGeneList = getUniqueGenes3(celltype="Global", Donor=i, Protocol = "MultiDonor", maxPadj=0.1, minFC = 0.2),
          myTitle = paste0("Unique DEGs in MultiDonor ",i,"\nRanked in SingleDonor ",i), myColor=DonorCols[i])
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Unique_GenesFrom_MD_",i,"_Global_in_SD_Feb5.jpeg"),
         height=2.8, width=8)
  plotKS5(queryRes = VPA_Res$Global[[i]]$SingleDonor$Results,
          myGeneList = getUniqueGenes3(celltype="Global", Donor=i, Protocol = "MultiDonor", maxPadj=0.1, minFC = 0.2),
          myTitle = paste0("Unique DEGs in MultiDonor ",i,"\nRanked in SingleDonor ",i), myColor=DonorCols[i])
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Unique_GenesFrom_MD_",i,"_Global_in_SD_Feb5_NoLabel.jpeg"),
                height=2.5, width=8)
  print(paste("Done with",i))
}

for(i in myDonors){
  print(paste("Starting",i))
  plotKS4(queryRes = VPA_Res$Global[[i]]$SingleDonor$Results,
          myGeneList = rownames(getSigGenes(VPA_Res$Global[[i]]$MultiDonor)),
          myTitle = paste0("Unique DEGs in MultiDonor ",i,"\nRanked in SingleDonor ",i))
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Sig_GenesFrom_MD_",i,"_Global_in_SD.svg"),
         height=2.5, width=8)
  plotKS5(queryRes = VPA_Res$Global[[i]]$SingleDonor$Results,
          myGeneList = rownames(getSigGenes(VPA_Res$Global[[i]]$MultiDonor)),
          myTitle = paste0("Unique DEGs in MultiDonor ",i,"\nRanked in SingleDonor ",i))
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Sig_GenesFrom_MD_",i,"_Global_in_SD.jpeg"),
         height=2.5, width=8)
  print(paste("Done with",i))
}

for(i in myDonors){
  print(paste("Starting",i))
  plotKS4(queryRes = VPA_Res$Global[[i]]$MultiDonor$Results,
          myGeneList = getUniqueGenes3(celltype="Global", Donor=i, Protocol = "SingleDonor", maxPadj=0.1, minFC = 0.2),
          myTitle = paste0("Unique DEGs in SingleDonor ",i,"\nRanked in MultiDonor ",i), myColor = DonorCols[i])
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Unique_GenesFrom_SD_",i,"_Global_in_MD_Feb5.jpeg"),
         height=2.8, width=8)
  plotKS5(queryRes = VPA_Res$Global[[i]]$MultiDonor$Results,
          myGeneList = getUniqueGenes3(celltype="Global", Donor=i, Protocol = "SingleDonor", maxPadj=0.1, minFC = 0.2),
          myTitle = paste0("Unique DEGs in SingleDonor ",i,"\nRanked in MultiDonor ",i), myColor = DonorCols[i])
  ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/KS_Plot_Unique_GenesFrom_SD_",i,"_Global_in_MD_Feb5_NoLabel.jpeg"),
         height=2.5, width=8)
  print(paste("Done with",i))
}

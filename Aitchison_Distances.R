## Aitchison Distance
library(plyr)
library(robCompositions)
library(reshape2)
library(ggh4x)

treatCols <- c()
treatCols["EtOH"] <- "#B497D4"
treatCols["VPA"] <- "#F2CC54"
treatCols["Control"] <- "#EAEAEA"
TypeCols <- c("#8dd3c7", "#fff7bc", "#bebada", "#fcbba1", "#80b1d3", "#fdb462",
              "#fa9fb5", "#04cf18", "#a6d854", "#fccde5", "#d9d9d9")
names(TypeCols) <- c("aRG","oRG","IP","PN","CFuPN","CPN",
                     "IN","ChPL","CR","Cortical hem","Unknown")
DonorCols <- c("darkred","#0A9396","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702")
names(DonorCols) <- c("11a","CW","CIRM17","GM","H1","PGP1","Mito210")
myCols2 <- DonorCols

`%ni%` = Negate(`%in%`)
ProtCols <- c()
ProtCols["Single"] <- "#8cdb86"
ProtCols["Multi"] <- "#ce367c"
ProtCols["SingleDonor"] <- ProtCols["Single"]
ProtCols["MultiDonor"] <- ProtCols["Multi"]
ProtCols["SingleDonorNSC"] <- ProtCols["Single"]
ProtCols["MultiDonorNSC"] <- ProtCols["Multi"]
ProtCols["NPC"] <- "orange"
ProtCols["Atlas"] <- "#196099"

## Calculate some Aitchison Distances:

meta <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Metadata_AllMerged_Streamlined_November17_2023.RDS")
tmp <- table(meta$Individual, meta$MergedAnnotation) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)

myTab <- reshape2::melt(tmp)
colnames(myTab) <- c("Sample", "CellType", "Percent")
myTab$Sample <- as.character(myTab$Sample)
myTab$Batch <- sapply(myTab$Sample, function(x){strsplit(x,"_")[[1]][1]})
myTab$Mix <- sapply(myTab$Sample, function(x){strsplit(x,"_")[[1]][2]})
myTab$Treatment <- sapply(myTab$Sample, function(x){strsplit(x,"_")[[1]][3]})
myTab$Organoid <- sapply(myTab$Sample, function(x){paste(strsplit(x,"_")[[1]][1:4], collapse="_")})
myTab$Donor <- sapply(myTab$Sample, function(x){strsplit(x,"_")[[1]][5]})
myTab$CellType <- as.character(myTab$CellType)

Aitch <- aDist(tmp)

a2 <- as.matrix(Aitch)

aTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(a2)-1)){
  for(j in (i+1):ncol(a2)){
    aTab <- rbind(aTab, c(rownames(a2)[i], colnames(a2)[j], a2[i,j]))
  }
}

aTab <- aTab[-1,]
aTab$Treat1 <- sapply(aTab$Sample1, function(x){strsplit(x,"_")[[1]][3]})
aTab$Treat2 <- sapply(aTab$Sample2, function(x){strsplit(x,"_")[[1]][3]})
aTab$Donor1 <- sapply(aTab$Sample1, function(x){strsplit(x,"_")[[1]][5]})
aTab$Donor2 <- sapply(aTab$Sample2, function(x){strsplit(x,"_")[[1]][5]})
aTab$Batch1 <- sapply(aTab$Sample1, function(x){strsplit(x,"_")[[1]][1]})
aTab$Batch2 <- sapply(aTab$Sample2, function(x){strsplit(x,"_")[[1]][1]})
aTab$Mix1 <- sapply(aTab$Sample1, function(x){strsplit(x,"_")[[1]][2]})
aTab$Mix2 <- sapply(aTab$Sample2, function(x){strsplit(x,"_")[[1]][2]})
aTab$Organoid1 <- sapply(aTab$Sample1, function(x){paste(strsplit(x,"_")[[1]][1:4], collapse="_")})
aTab$Organoid2 <- sapply(aTab$Sample2, function(x){paste(strsplit(x,"_")[[1]][1:4], collapse="_")})
aTab$Protocol1 <- sapply(aTab$Sample1, function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})
aTab$Protocol2 <- sapply(aTab$Sample2, function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})
aTab$Distance <- as.numeric(aTab$Distance)
aTab$Within_Treatment <- aTab$Treat1 == aTab$Treat2
aTab$Within_Donor <- aTab$Donor1 == aTab$Donor2
aTab$Within_Batch <- aTab$Batch1 == aTab$Batch2
aTab$Within_Protocol <- aTab$Protocol1 == aTab$Protocol2

## Check whether treatment changes within-group distance
ggplot(aTab[aTab$Treat1 == aTab$Treat2 & aTab$Protocol1=="Single" & aTab$Protocol2=="Single",],
       aes(x=Treat1, y=Distance, fill=Treat1)) +
  geom_boxplot() +
  scale_fill_manual(values=treatCols[unique(aTab$Treat1)]) +
  ggtitle("Within-Treatment Aitchison Distance\nWithin Single-Donor") +
  theme(plot.title = element_text(face="bold", hjust=0.5))
## Check whether treatment changes within-group distance
ggplot(aTab[aTab$Treat1 == aTab$Treat2 & aTab$Protocol1=="Multi" & aTab$Protocol2=="Multi",],
       aes(x=Treat1, y=Distance, fill=Treat1)) +
  geom_boxplot() +
  scale_fill_manual(values=treatCols[unique(aTab$Treat1)]) +
  ggtitle("Within-Treatment Aitchison Distance\nWithin Multi-Donor") +
  theme(plot.title = element_text(face="bold", hjust=0.5))

## For each treatment, check the within-group distance between protocols:
ggplot(aTab[aTab$Treat1 == aTab$Treat2 & aTab$Protocol1 == aTab$Protocol2,],
       aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  scale_fill_manual(values=ProtCols[unique(aTab$Protocol1)]) +
  facet_grid(~Treat1) +
  ggtitle("Within-Treatment Aitchison Distance\nComparing Within-Single-Donor to Within-Multi-Donor") +
  theme(plot.title = element_text(face="bold", hjust=0.5))

## Within controls only, *within each donor*, compare within-group single-donor to multi-donor distances:
ggplot(aTab[aTab$Treat1 == "Control" & aTab$Treat2 == "Control" & aTab$Protocol1 == aTab$Protocol2 &
              aTab$Donor1 == aTab$Donor2,],
       aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  scale_fill_manual(values=ProtCols[unique(aTab$Protocol1)]) +
  facet_grid(~Donor1) +
  ggtitle("Within-Treatment Aitchison Distance\nWithin Each Donor - SD vs MD") +
  theme(plot.title = element_text(face="bold", hjust=0.5))

## Just for a moment, look only at the multidonor -
## Compare within-batch/within-donor distance to within-batch/between-donor
## And also to between-batch/within-donor and between-batch/between-donor:
aTab$Temporary <- aTab$Within_Batch + 2*aTab$Within_Donor
aTab$Temporary <- as.character(mapvalues(aTab$Temporary, from=c(3,2,1,0),
                                         to=c("Within-Batch\nWithin-Donor", "Between-Batch\nWithin-Donor",
                                              "Within-Batch\nBetween-Donor", "Between-Batch\nBetween-Donor")))
aTab$Temporary <- factor(aTab$Temporary, levels=c("Within-Batch\nWithin-Donor", "Between-Batch\nWithin-Donor",
                                                  "Within-Batch\nBetween-Donor", "Between-Batch\nBetween-Donor"))
ggplot(aTab[aTab$Treat1 == "Control" & aTab$Within_Treatment & aTab$Within_Protocol,],
       aes(x=Temporary, y=Distance, fill=Protocol1)) +
  geom_violin() +
  geom_boxplot(width=0.3) +
  facet_grid(~Protocol1) +
  theme(panel.background=element_rect(fill="white", color="black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey85", linewidth=0.2, linetype=2))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_MD_SD_WithinBetween_BatchAndDonor_",Sys.Date(),".png"), height=5, width=10)

ggplot(aTab[aTab$Treat1 == "Control" & aTab$Treat2 != "EtOH" & aTab$Within_Protocol,],
       aes(x=Temporary, y=Distance, fill=Treat2)) +
  geom_boxplot(width=0.5) +
  facet_grid(Protocol1~1) +
  scale_fill_manual(name="Within-Controls?", values=c(treatCols["Control"], treatCols["VPA"])) +
  theme(panel.background=element_rect(fill="white", color="black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey85", linewidth=0.2, linetype=2))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_MD_SD_WithinBetween_BatchAndDonor_WithinControl-BetweenConVPA_",Sys.Date(),".png"), height=8, width=6)

aTab$Donor1[aTab$Donor1=="CIRM17"] <- "CW"
aTab$Donor2[aTab$Donor2=="CIRM17"] <- "CW"


tmp <- aTab[aTab$Within_Treatment & aTab$Treat1 != "EtOH" &
              ((aTab$Within_Donor & aTab$Donor1=="PGP1" & aTab$Within_Protocol & aTab$Protocol1=="Multi") |
                 (aTab$Within_Protocol & aTab$Protocol1=="Multi" & !aTab$Within_Donor & aTab$Donor1 %in% c("PGP1","CW") & aTab$Donor2 %in% c("PGP1","CW")) |
                 (!aTab$Within_Protocol & aTab$Donor1=="PGP1" & aTab$Within_Donor) |
                 (!aTab$Within_Protocol & !aTab$Within_Donor & aTab$Donor1 %in% c("PGP1","CW") & aTab$Donor2 %in% c("PGP1","CW"))),]
tmp$Temporary <- tmp$Within_Donor + tmp$Within_Protocol*2
tmp$Temporary <- as.character(mapvalues(tmp$Temporary,from=c(3,2,1,0),
                                        to=c("PGP1-PGP1 Distance\nWithin MultiDonors",
                                             "PGP1-CW Distance\nWithin MultiDonors",
                                             "PGP1-PGP1 Distance\nMultiDonor-to-SingleDonor",
                                             "PGP1-CW Distance\nMultiDonor-to-SingleDonor")))
## Focus in, look at whether the protocol makes a bigger difference within a specific donor than the donor identity:
ggplot(tmp,
       aes(x=Temporary, y=Distance)) +
  geom_boxplot() +
  facet_grid(Treat1~1)



library(ComplexHeatmap)
topAnno <- HeatmapAnnotation(Treatment=as.character(sapply(colnames(a2), function(x){strsplit(x,"_")[[1]][3]})),
                             Donor=as.character(sapply(colnames(a2), function(x){strsplit(x,"_")[[1]][5]})),
                             Batch=as.character(sapply(colnames(a2), function(x){strsplit(x,"_")[[1]][1]})),
                             Protocol=as.character(sapply(colnames(a2), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                             col = list(Treatment=treatCols,
                                        Donor=myCols2,
                                        Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                        Protocol=ProtCols))
sideAnno <- rowAnnotation(Treatment=as.character(sapply(rownames(a2), function(x){strsplit(x,"_")[[1]][3]})),
                          Donor=as.character(sapply(rownames(a2), function(x){strsplit(x,"_")[[1]][5]})),
                          Batch=as.character(sapply(rownames(a2), function(x){strsplit(x,"_")[[1]][1]})),
                          Protocol=as.character(sapply(rownames(a2), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                          col = list(Treatment=treatCols,
                                     Donor=myCols2,
                                     Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                     Protocol=ProtCols))
Heatmap(a2, top_annotation = topAnno, left_annotation = sideAnno,
        column_order = c(1:ncol(a2)),
        row_order = c(1:nrow(a2)),
        show_column_names = FALSE, show_row_names = FALSE)
Heatmap(a2, top_annotation = topAnno, left_annotation = sideAnno,
        show_column_names = FALSE, show_row_names = FALSE)

NoEth <- a2[sapply(rownames(a2), function(x){strsplit(x,"_")[[1]][3]}) != "EtOH",
             sapply(colnames(a2), function(x){strsplit(x,"_")[[1]][3]}) != "EtOH"]
topAnno2 <- HeatmapAnnotation(Treatment=as.character(sapply(colnames(NoEth), function(x){strsplit(x,"_")[[1]][3]})),
                             Donor=as.character(sapply(colnames(NoEth), function(x){strsplit(x,"_")[[1]][5]})),
                             Batch=as.character(sapply(colnames(NoEth), function(x){strsplit(x,"_")[[1]][1]})),
                             Protocol=as.character(sapply(colnames(NoEth), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                             col = list(Treatment=treatCols,
                                        Donor=myCols2,
                                        Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                        Protocol=ProtCols))
sideAnno2 <- rowAnnotation(Treatment=as.character(sapply(rownames(NoEth), function(x){strsplit(x,"_")[[1]][3]})),
                          Donor=as.character(sapply(rownames(NoEth), function(x){strsplit(x,"_")[[1]][5]})),
                          Batch=as.character(sapply(rownames(NoEth), function(x){strsplit(x,"_")[[1]][1]})),
                          Protocol=as.character(sapply(rownames(NoEth), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                          col = list(Treatment=treatCols,
                                     Donor=myCols2,
                                     Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                     Protocol=ProtCols))
Heatmap(NoEth, top_annotation = topAnno2, left_annotation = sideAnno2,
        show_column_names = FALSE, show_row_names = FALSE)

Heatmap(NoEth, top_annotation = topAnno2, left_annotation = sideAnno2,
        column_order = order(sapply(colnames(NoEth), function(x){strsplit(x,"_")[[1]][3]}),sapply(colnames(NoEth), function(x){strsplit(x,"_")[[1]][5]})),
        row_order = order(sapply(rownames(NoEth), function(x){strsplit(x,"_")[[1]][3]}),sapply(rownames(NoEth), function(x){strsplit(x,"_")[[1]][5]})),
        show_column_names = FALSE, show_row_names = FALSE)

test <- table(meta$Individual, meta$Batch)

p1 <- ggplot(aTab[aTab$Within_Batch & aTab$Within_Treatment & aTab$Within_Donor & aTab$Treat1 != "EtOH",],
       aes(x=Donor1, y=Distance, fill=Donor1)) +
  geom_boxplot() +
  geom_jitter(pch=21, size=4) +
  facet_grid(Protocol1~Treat1) +
  scale_fill_manual(values=myCols2[unique(aTab$Donor1)]) +
  ylim(c(0,10)) +
  ggtitle("Within-batch distances")

p2 <- ggplot(aTab[!aTab$Within_Batch & aTab$Within_Treatment & aTab$Within_Donor
            & aTab$Within_Protocol & aTab$Treat1 != "EtOH",],
       aes(x=Donor1, y=Distance, fill=Donor1)) +
  geom_boxplot() +
  geom_jitter(pch=21, size=4) +
  facet_grid(Protocol1~Treat1) +
  scale_fill_manual(values=myCols2[unique(aTab$Donor1)]) +
  ylim(c(0,10)) +
  ggtitle("Between-batch distances")

cowplot::plot_grid(p1, p2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithinBatch-BetweenBatch_",Sys.Date(),".png"), height=8, width=14)


## Let's examine an example: Mito210-PGP1 differences:
## First, explore within-batch differences, compared across protocols (within controls)
ggplot(aTab[aTab$Within_Protocol & aTab$Within_Batch & aTab$Within_Treatment & aTab$Treat1!="EtOH" &
              ((aTab$Donor1=="Mito210" & aTab$Donor2=="PGP1") | (aTab$Donor1=="PGP1" & aTab$Donor2=="Mito210")),],
       aes(x=Batch1, y=Distance, fill=Batch1)) +
  facet_grid(Treat1~Protocol1, space="free", scales="free") +
  geom_boxplot() +
  geom_jitter(pch=21, size=3) +
  ggtitle("PGP1-Mito210 Aitchison Distances\nWithin-Batch") +
  theme(plot.title=element_text(face="bold", hjust=0.5))

## Compare to within-donor:
ggplot(aTab[aTab$Within_Protocol & aTab$Within_Batch & aTab$Within_Treatment & aTab$Treat1!="EtOH" &
              aTab$Donor1 %in% c("Mito210", "PGP1") & aTab$Donor2 %in% c("Mito210","PGP1"),],
       aes(x=Within_Donor, y=Distance)) +
  facet_grid(Treat1~Protocol1, space="free", scales="free") +
  geom_boxplot() +
  geom_jitter(pch=21, size=3) +
  ggtitle("PGP1-Mito210 Aitchison Distances\nWithin-Batch") +
  theme(plot.title=element_text(face="bold", hjust=0.5))

atlas <- readRDS("/Users/tfaits/Desktop/Documents/Atlas_RNA_data/ForLabelTransfer/Metadata_3mo.RDS")

tmp <- table(atlas$org, atlas$FinalName) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)

Aitch <- aDist(tmp)

atlas2 <- as.matrix(Aitch)

atlasTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(atlas2)-1)){
  for(j in (i+1):ncol(atlas2)){
    atlasTab <- rbind(atlasTab, c(rownames(atlas2)[i], colnames(atlas2)[j], atlas2[i,j]))
  }
}

atlasTab <- atlasTab[-1,]
atlasTab$Donor1 <- sapply(atlasTab$Sample1, function(x){strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1]})
atlasTab$Donor2 <- sapply(atlasTab$Sample2, function(x){strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1]})
atlasTab$Batch1 <- sapply(atlasTab$Sample1, function(x){paste(strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1:2], collapse="-")})
atlasTab$Batch2 <- sapply(atlasTab$Sample2, function(x){paste(strsplit(strsplit(x,"[.]")[[1]][2],"_")[[1]][1:2], collapse="-")})
atlasTab$Batch1[grep("Mito210_3m_PTEN", atlasTab$Sample1)] <- "Mit210-PTEN"
atlasTab$Batch2[grep("Mito210_3m_PTEN", atlasTab$Sample2)] <- "Mit210-PTEN"

atlasTab$Distance <- as.numeric(atlasTab$Distance)

atlasTab$Within_Donor <- atlasTab$Donor1 == atlasTab$Donor2
atlasTab$Within_Batch <- atlasTab$Batch1 == atlasTab$Batch2
atlasTab$Protocol1 <- "Atlas"
atlasTab$Protocol2 <- "Atlas"

#TempThing1 <- atlasTab[atlasTab$Within_Donor & atlasTab$Donor1 %in% c("PGP1","Mito210"),]
TempThing1 <- atlasTab
TempThing2 <- aTab[aTab$Within_Protocol & aTab$Within_Treatment & aTab$Treat1=="Control",colnames(TempThing1)]

TempThing <- rbind(TempThing1, TempThing2)
myCols2 <- c("#AE2012","darkred","#0A9396","#0A9396","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012","cyan")
names(myCols2) <- c("HUES66","11a","CIRM17","CW","CIRM","GM","H1","PGP1","Mito210","HUES66_AC2","Fetal")
ggplot(TempThing[TempThing$Within_Batch & TempThing$Within_Donor,], aes(x=Donor1, y=Distance, fill=Donor1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=myCols2[unique(TempThing$Donor1)]) +
  ggtitle("Within-Batch Aitchison Distances\nBy Donor") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlas_WithinBatch_",Sys.Date(),".png"), height=5, width=8)

ProtCols["Atlas"] <- "#196099"
ggplot(TempThing[!TempThing$Within_Batch,], aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=ProtCols[unique(TempThing$Protocol1)]) +
  ggtitle("Between-Batch Aitchison Distances") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlas_BetweenBatch_",Sys.Date(),".png"), height=5, width=5)

ggplot(TempThing[TempThing$Within_Batch & TempThing$Within_Donor,], aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=ProtCols[unique(TempThing$Protocol1)]) +
  ggtitle("Within-Batch (and Donor)\nAitchison Distances") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlas_WithinBatchDonor_",Sys.Date(),".png"), height=5, width=5)



czi <- readRDS("/Users/tfaits/Documents/Chimera/Seurats/Finalized_EtOH_VPA_Objects/CZI_Metadata.RDS")

tmp <- table(czi$orig.ident, czi$FinalName) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)

Aitch <- aDist(tmp)

czi2 <- as.matrix(Aitch)

cziTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(czi2)-1)){
  for(j in (i+1):ncol(czi2)){
    cziTab <- rbind(cziTab, c(rownames(czi2)[i], colnames(czi2)[j], czi2[i,j]))
  }
}

cziTab <- cziTab[-1,]
cziTab$Donor1 <- "Fetal"
cziTab$Donor2 <- "Fetal"
cziTab$Batch1 <- cziTab$Sample1
cziTab$Batch2 <- cziTab$Sample2

cziTab$Distance <- as.numeric(cziTab$Distance)

cziTab$Within_Donor <- TRUE
cziTab$Within_Batch <- TRUE
cziTab$Protocol1 <- "Fetal"
cziTab$Protocol2 <- "Fetal"

#TempThing1 <- atlasTab[atlasTab$Within_Donor & atlasTab$Donor1 %in% c("PGP1","Mito210"),]
TempThing3 <- cziTab
TempThing <- rbind(TempThing1, TempThing2, TempThing3)

ggplot(TempThing[TempThing$Within_Batch & TempThing$Within_Donor,], aes(x=Donor1, y=Distance, fill=Donor1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=myCols2[unique(TempThing$Donor1)]) +
  ggtitle("Within-Batch Aitchison Distances\nBy Donor") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlasCZI_WithinBatch_",Sys.Date(),".png"), height=5, width=9)

ggplot(TempThing[!TempThing$Within_Batch,], aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=ProtCols[unique(TempThing$Protocol1)]) +
  ggtitle("Between-Batch Aitchison Distances") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlas_BetweenBatch_",Sys.Date(),".png"), height=5, width=5)

ProtCols["Fetal"] <- "cyan"
ggplot(TempThing[TempThing$Within_Batch & TempThing$Within_Donor,], aes(x=Protocol1, y=Distance, fill=Protocol1)) +
  geom_boxplot() +
  facet_grid(~Protocol1, scales="free", space="free") +
  scale_fill_manual(values=ProtCols[unique(TempThing$Protocol1)]) +
  ggtitle("Within-Batch (and Donor)\nAitchison Distances") +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_WithAtlasCZI_WithinBatchDonor_",Sys.Date(),".png"), height=5, width=6)



Trevino <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Trevino_Metadata.RDS")
tmp <- table(Trevino$Sample.ID, Trevino$CellTypes) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)
Aitch <- aDist(tmp)
Trevino2 <- as.matrix(Aitch)
TrevinoTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(Trevino2)-1)){
  for(j in (i+1):ncol(Trevino2)){
    TrevinoTab <- rbind(TrevinoTab, c(rownames(Trevino2)[i], colnames(Trevino2)[j], Trevino2[i,j]))
  }
}
TrevinoTab <- TrevinoTab[-1,]
TrevinoTab$Donor1 <- "Fetal"
TrevinoTab$Donor2 <- "Fetal"
TrevinoTab$Batch1 <- TrevinoTab$Sample1
TrevinoTab$Batch2 <- TrevinoTab$Sample2
TrevinoTab$Distance <- as.numeric(TrevinoTab$Distance)
TrevinoTab$Within_Donor <- TRUE
TrevinoTab$Within_Batch <- TRUE
TrevinoTab$Protocol1 <- "Fetal"
TrevinoTab$Protocol2 <- "Fetal"

Polioudakis <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Polioudakis_Metadata.RDS")
tmp <- table(Polioudakis$Donor, Polioudakis$Cluster) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)
Aitch <- aDist(tmp)
Polioudakis2 <- as.matrix(Aitch)
PolioudakisTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(Polioudakis2)-1)){
  for(j in (i+1):ncol(Polioudakis2)){
    PolioudakisTab <- rbind(PolioudakisTab, c(rownames(Polioudakis2)[i], colnames(Polioudakis2)[j], Polioudakis2[i,j]))
  }
}
PolioudakisTab <- PolioudakisTab[-1,]
PolioudakisTab$Donor1 <- "Fetal"
PolioudakisTab$Donor2 <- "Fetal"
PolioudakisTab$Batch1 <- PolioudakisTab$Sample1
PolioudakisTab$Batch2 <- PolioudakisTab$Sample2
PolioudakisTab$Distance <- as.numeric(PolioudakisTab$Distance)
PolioudakisTab$Within_Donor <- TRUE
PolioudakisTab$Within_Batch <- TRUE
PolioudakisTab$Protocol1 <- "Fetal"
PolioudakisTab$Protocol2 <- "Fetal"




npc <- readRDS("/Users/tfaits/Documents/Chimera/Seurats/Finalized_EtOH_VPA_Objects/D23_D15_M3M4_Trimmed_Metadata.RDS")
#npc <- npc[npc$Genotype %in% c("PGP1","Mito210") & npc$orig.ident %in% c("August2022_Mix1_D23_1","August2022_Mix1_D23_2"),]
npc <- npc[npc$orig.ident %in% c("August2022_Mix1_D23_1","August2022_Mix1_D23_2"),]

npc$Individual <- paste(npc$Sample, npc$Genotype, sep="-")
tmp <- table(npc$Individual, npc$CellType) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)

Aitch <- aDist(tmp)

npc2 <- as.matrix(Aitch)

npcTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(npc2)-1)){
  for(j in (i+1):ncol(npc2)){
    npcTab <- rbind(npcTab, c(rownames(npc2)[i], colnames(npc2)[j], npc2[i,j]))
  }
}
npcTab <- npcTab[-1,]
npcTab$Sample1 <- as.character(npcTab$Sample1)
npcTab$Sample2 <- as.character(npcTab$Sample2)
npcTab$Donor1 <- sapply(npcTab$Sample1, function(x){strsplit(x,"-")[[1]][2]})
npcTab$Donor2 <- sapply(npcTab$Sample2, function(x){strsplit(x,"-")[[1]][2]})
npcTab$Batch1 <- "NPC_A"
npcTab$Batch2 <- "NPC_A"

npcTab$Distance <- as.numeric(npcTab$Distance)

npcTab$Within_Donor <- npcTab$Donor1==npcTab$Donor2
npcTab$Within_Batch <- TRUE
npcTab$Protocol1 <- "NPC"
npcTab$Protocol2 <- "NPC"

#TempThing1 <- atlasTab[atlasTab$Within_Donor & atlasTab$Donor1 %in% c("PGP1","Mito210"),]
TempThing <- rbind(TempThing1, TempThing2, npcTab)

TempThing <- TempThing[TempThing$Donor1 %in% c("PGP1", "Mito210") & TempThing$Donor2 %in% c("PGP1","Mito210"),]
TempThing$Donors <- paste(TempThing$Donor1, TempThing$Donor2, sep="-")
TempThing$Donors[TempThing$Donors=="Mito210-PGP1"] <- "PGP1-Mito210"
TempThing$Batches <- ifelse(TempThing$Within_Batch, "Within Batch", "Between Batches")
TempThing$Donors <- factor(TempThing$Donors, levels=c("PGP1-PGP1","Mito210-Mito210","PGP1-Mito210"))
TempThing$Batches <- factor(TempThing$Batches, levels=c("Within Batch","Between Batches"))
TempThing$Protocol1 <- factor(TempThing$Protocol1, levels=c("Multi","Single","NPC","Atlas"))
ProtoCol <- c("#bbe4ea", "#a0cee7", "#a0b1db", "#9484be")
names(ProtoCol) <- c("Atlas","Single","Multi","NPC")
ggplot(TempThing, aes(x=Donors, y=Distance, fill=Protocol1)) +
  geom_boxplot(outlier.shape=NA, outlier.color="green") +
  facet_nested(~Batches + Protocol1) +
  scale_fill_manual(values=ProtCols[levels(TempThing$Protocol1)]) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(color="grey75", linewidth=0.3)) +
  geom_jitter(width=0.15, aes(color=Protocol1)) +
  scale_color_manual(values=ProtoCol[levels(TempThing$Protocol1)])
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_AllProtocols_PGP1Mito210_WithinBetweenBatch_",Sys.Date(),".svg"),
       height=5, width=8)


TT2 <- TempThing
PolioudakisTab$Donors <- "Polioudakis et al"
PolioudakisTab$Batches <- "Fetal Datasets"
cziTab$Donors <- "Uzquiano et al"
cziTab$Batches <- "Fetal Datasets"
TrevinoTab$Donors <- "Trevino et al"
TrevinoTab$Batches <- "Fetal Datasets"
TT2 <- rbind(TT2, PolioudakisTab, cziTab, TrevinoTab)

ProtCols <- c("#bbe4ea", "#a0cee7", "#a0b1db", "#9484be", "cyan")
names(ProtCols) <- c("Atlas","Single","Multi","NPC", "Fetal")
levels(TT2$Protocol1)
ggplot(TT2, aes(x=Donors, y=Distance, fill=Protocol1)) +
  geom_boxplot(outlier.shape=NA, outlier.color="green") +
  facet_nested(~Batches + Protocol1, scales="free", space="free") +
  scale_fill_manual(values=ProtCols[levels(TT2$Protocol1)]) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(color="grey75", linewidth=0.3)) +
  geom_jitter(width=0.15, color="black", alpha=0.8)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_AllProtocols_PGP1Mito210_WithinBetweenBatch_",Sys.Date(),".svg"),
       height=5, width=8)





# npc <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/M3M4_NPC_NSC_UpdatedAnnotations_Metadata.RDS")
# npc <- npc[npc$orig.ident %in% c("August2022_Mix1_D23_1","August2022_Mix1_D23_2"),]
# npc$Individual <- paste(npc$orig.ident, npc$Genotype, sep="-")
# npc$CellType <- npc$MergedAnnotation
# tmp <- table(npc$Individual, npc$CellType) + 0.1
# tmp <- 100*tmp/apply(tmp,1,sum)
# 
# Aitch <- aDist(tmp)
# 
# npc2 <- as.matrix(Aitch)
# 
# npcTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
# for(i in 1:(nrow(npc2)-1)){
#   for(j in (i+1):ncol(npc2)){
#     npcTab <- rbind(npcTab, c(rownames(npc2)[i], colnames(npc2)[j], npc2[i,j]))
#   }
# }
# npcTab <- npcTab[-1,]
# npcTab$Sample1 <- as.character(npcTab$Sample1)
# npcTab$Sample2 <- as.character(npcTab$Sample2)
# npcTab$Donor1 <- sapply(npcTab$Sample1, function(x){strsplit(x,"-")[[1]][2]})
# npcTab$Donor2 <- sapply(npcTab$Sample2, function(x){strsplit(x,"-")[[1]][2]})
# npcTab$Batch1 <- "NPC_A"
# npcTab$Batch2 <- "NPC_A"
# 
# npcTab$Distance <- as.numeric(npcTab$Distance)
# 
# npcTab$Within_Donor <- npcTab$Donor1==npcTab$Donor2
# npcTab$Within_Batch <- TRUE
# npcTab$Protocol1 <- "NPC"
# npcTab$Protocol2 <- "NPC"

## Points - jitter plot showing the within-batch/within-donor Aitchison distance (should be a direct replacement for Fig 3 F)
TempThing <- rbind(TempThing1, TempThing2, npcTab)
TempThing <- TempThing[TempThing$Donor1==TempThing$Donor2 & TempThing$Batch1==TempThing$Batch2,]
TempThing$Protocol1 <- factor(TempThing$Protocol1, levels=c("Multi","Single","NPC","Atlas"))
DonorCols["HUES66"] <- "red"
ggplot(TempThing, aes(x=Protocol1, y=Distance, fill=Donor1)) +
  geom_jitter(pch=21, width=0.2, size=2) +
  scale_fill_manual(values=DonorCols[unique(TempThing$Donor1)]) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(color="grey75", linewidth=0.3))
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_AllProtocols_WithinBatchAndDonor_",Sys.Date(),".svg"),
       height=4, width=4)

mean(cziTab$Distance)


ggplot(TempThing, aes(x=Protocol1, y=Distance, fill=Donor1)) +
  geom_jitter(pch=21, width=0.2, size=2) +
  scale_fill_manual(values=DonorCols[unique(TempThing$Donor1)], name="Donor") +
  ylim(c(0,9)) +
  xlab("Protocol") +
  ylab("Aitchison Distance") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(color="grey75", linewidth=0.3)) +
  geom_abline(intercept=mean(cziTab$Distance), slope=0, linetype="dashed",
              color="#696969", linewidth=1) +
  geom_abline(intercept=mean(TrevinoTab$Distance), slope=0, linetype="dashed",
              color="#3B3B3B", linewidth=1) +
  geom_abline(intercept=mean(PolioudakisTab$Distance), slope=0, linetype="dashed",
              color="#aeaeae", linewidth=1) +
  annotate(geom="text", label=c("CZI","Trevino","Polioudakis"), x=4.5,
           y=c(mean(cziTab$Distance), mean(TrevinoTab$Distance), mean(PolioudakisTab$Distance)),
           vjust=-1, size=3, hjust=1)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_AllProtocols_WithinBatchAndDonor_WithFetals_",Sys.Date(),".svg"),
       height=4, width=4)

## Points - jitter plot showing the within-donor (ignore batches?) Aitchison distance (May replace Fig 3 F)
TempThing <- rbind(TempThing1, TempThing2, npcTab)
TempThing <- TempThing[TempThing$Donor1==TempThing$Donor2 & TempThing$Protocol1==TempThing$Protocol2,]
TempThing$Protocol1 <- factor(TempThing$Protocol1, levels=c("Multi","Single","NPC","Atlas"))
DonorCols["HUES66"] <- "red"
ggplot(TempThing, aes(x=Protocol1, y=Distance, fill=Donor1)) +
  geom_jitter(pch=21, width=0.2, size=2) +
  xlab("Protocol") + ylab("Aitchison Distance\nBetween Replicates") +
  scale_fill_manual(values=DonorCols[unique(TempThing$Donor1)], name="Donor") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background = element_rect(fill="white"),
          panel.grid.major.y = element_line(color="grey75", linewidth=0.3)) +
  geom_abline(intercept=mean(cziTab$Distance), slope=0, linetype="dashed", color="#696969", linewidth=1) +
  geom_abline(intercept=mean(TrevinoTab$Distance), slope=0, linetype="dashed", color="#3B3B3B", linewidth=1) +
  geom_abline(intercept=mean(PolioudakisTab$Distance), slope=0, linetype="dashed", color="#aeaeae", linewidth=1) +
  annotate(geom="text", label=c("CZI","Trevino","Polioudakis"), x=4.5,
           y=c(mean(cziTab$Distance), mean(TrevinoTab$Distance), mean(PolioudakisTab$Distance)),
           vjust=-1, size=3, hjust=1)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_Boxplot_AllProtocols_WithinDonorOnly_WithFetal_",Sys.Date(),".svg"),
         height=4, width=4)



meta$Individual <- paste(meta$Individual, meta$Genotype, sep="-")
NSCs <- meta[,c("Individual","MergedAnnotation","Protocol")]
npc$MergedAnnotation <- npc$CellType
NPCs <- npc[,c("Individual","MergedAnnotation")]
NPCs$Protocol <- "NPC"
atlas$Individual <- mapvalues(atlas$org,
                              from=c("DorsalKadoshima.GM_3mon_1","DorsalKadoshima.GM_3mon_2","DorsalKadoshima.GM_3mon_3",
                                     "DorsalKadoshima.HUES66_3mon_1","DorsalKadoshima.HUES66_3mon_2","DorsalKadoshima.HUES66_3mon_3",    
                                     "DorsalKadosima.HUES66_3m_d106_1","DorsalKadosima.HUES66_3m_d106_2","DorsalKadosima.HUES66_3m_d106_3",  
                                     "DorsalKadoshima.Mito210_3m_1","DorsalKadoshima.Mito210_3m_2","DorsalKadoshima.Mito210_3m_3",     
                                     "DorsalKadoshima.Mito210_3m_PTEN_1","DorsalKadoshima.Mito210_3m_PTEN_2","DorsalKadoshima.Mito210_3m_PTEN_3",
                                     "DorsalKadoshima.PGP1_b1_3mon_1","DorsalKadoshima.PGP1_b1_3mon_2","DorsalKadoshima.PGP1_b1_3mon_3",   
                                     "DorsalKadoshima.PGP1_b2_3mon_1","DorsalKadoshima.PGP1_b2_3mon_2","DorsalKadoshima.PGP1_b2_3mon_3"),
                              to=c("GM_1-GM","GM_2-GM","GM_3-GM",
                                   "HUES66_a1-HUES66","HUES66_a2-HUES66","HUES66_a3-HUES66",
                                   "HUES66_b1-HUES66","HUES66_b2-HUES66","HUES66_b3-HUES66",
                                   "Mito210_a1-Mito210","Mito210_a2-Mito210","Mito210_a3-Mito210",
                                   "Mito210_b1-Mito210","Mito210_b2-Mito210","Mito210_b3-Mito210",
                                   "PGP1_a1-PGP1","PGP1_a2-PGP1","PGP1_a3-PGP1",
                                   "PGP1_b1-PGP1","PGP1_b2-PGP1","PGP1_b3-PGP1"))
atlas$MergedAnnotation <- atlas$FinalName
VELASCO <- atlas[,c("Individual","MergedAnnotation")]
VELASCO$Protocol <- "Velasco"
VELASCO$MergedAnnotation <- as.character(VELASCO$MergedAnnotation)
VELASCO$MergedAnnotation[VELASCO$MergedAnnotation == "Cajal Retzius"] <- "CR"
VELASCO$MergedAnnotation[VELASCO$MergedAnnotation == "Immature IN"] <- "IN"
FULL <- rbind(NSCs, NPCs, VELASCO)
tmp <- table(FULL$Individual, FULL$MergedAnnotation) + 0.1
tmp <- 100*tmp/apply(tmp,1,sum)
Aitch <- aDist(tmp)
full2 <- as.matrix(Aitch)
fullTab <- data.frame(Sample1="Fake", Sample2="AlsoFake", Distance=99)
for(i in 1:(nrow(full2)-1)){
  for(j in (i+1):ncol(full2)){
    fullTab <- rbind(fullTab, c(rownames(full2)[i], colnames(full2)[j], full2[i,j]))
  }
}
fullTab <- fullTab[-1,]
fullTab$Sample1 <- as.character(fullTab$Sample1)
fullTab$Sample2 <- as.character(fullTab$Sample2)
fullTab$Donor1 <- sapply(fullTab$Sample1, function(x){strsplit(x,"-")[[1]][2]})
fullTab$Donor2 <- sapply(fullTab$Sample2, function(x){strsplit(x,"-")[[1]][2]})
fullTab$Protocol1 <- sapply(fullTab$Sample1, function(x){FULL$Protocol[FULL$Individual==x][1]})
fullTab$Protocol2 <- sapply(fullTab$Sample2, function(x){FULL$Protocol[FULL$Individual==x][1]})

fullTab$Within_Protocol <- fullTab$Protocol1 == fullTab$Protocol2
fullTab$Within_Donor <- fullTab$Donor1 == fullTab$Donor2
fullTab$InvolvesMultidonor <- fullTab$Protocol1=="MultiDonorNSC" | fullTab$Protocol2=="MultiDonorNSC"
fullTab$OtherProt <- "NA"
for(i in which(fullTab$InvolvesMultidonor)){
  fullTab$OtherProt[i] <- ifelse(fullTab$Protocol1[i]=="MultiDonorNSC", fullTab$Protocol2[i], fullTab$Protocol1[i])
}
fullTab$OtherProt[fullTab$Within_Protocol] <- fullTab$Protocol1[fullTab$Within_Protocol]
fullTab$Distance <- as.numeric(fullTab$Distance)
fullTab$Treatment1 <- "Control"
fullTab$Treatment2 <- "Control"
fullTab$Treatment1[grep("VPA",fullTab$Sample1)] <- "VPA"
fullTab$Treatment2[grep("VPA",fullTab$Sample2)] <- "VPA"
fullTab$Treatment1[grep("EtOH",fullTab$Sample1)] <- "EtOH"
fullTab$Treatment2[grep("EtOH",fullTab$Sample2)] <- "EtOH"
fullTab$Within_Treatment <- fullTab$Treatment1==fullTab$Treatment2
fullTab$Control_Only <- fullTab$Treatment1=="Control" & fullTab$Treatment2 =="Control"

fullTab$OtherProt[fullTab$OtherProt=="SingleDonorNSC"] <- "NSC\nSingle\nDonor"
fullTab$OtherProt[fullTab$OtherProt=="MultiDonorNSC"] <- "NSC\nMultidonor"
fullTab$OtherProt <- factor(fullTab$OtherProt, levels=c("Velasco","NSC\nSingle\nDonor","NSC\nMultidonor","NPC"))
myLab <- c("Between-Protocols:\nMultiDonor NSC vs","Within-Protocol")
names(myLab) <- c("FALSE","TRUE")
ggplot(fullTab[fullTab$Within_Donor & (fullTab$InvolvesMultidonor | fullTab$Within_Protocol) & fullTab$Control_Only,],
       aes(x=OtherProt, y=Distance, fill=Donor1)) +
  geom_jitter(pch=21, width=0.3) +
  facet_grid(~Within_Protocol, space="free", scales="free", labeller=labeller(Within_Protocol=myLab)) +
  scale_fill_manual(values=DonorCols[unique(fullTab$Donor1)], name="Donor") +
  theme(panel.background = element_blank()) +
  xlab("") +
  geom_abline(intercept=mean(cziTab$Distance), slope=0, linetype="dashed", color="#696969", linewidth=1) +
  geom_abline(intercept=mean(TrevinoTab$Distance), slope=0, linetype="dashed", color="#3B3B3B", linewidth=1) +
  geom_abline(intercept=mean(PolioudakisTab$Distance), slope=0, linetype="dashed", color="#aeaeae", linewidth=1)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_BetweenProtocolsAndWithinProtocols_WithinDonor_",Sys.Date(),".svg"),
       height=4, width=6)

fT2a <- fullTab
fT2a$Within_Protocol[fT2a$Within_Protocol] <- "Within-Protocol"
fT2a$Within_Protocol[fT2a$Within_Protocol==FALSE] <- "Between-Protocols:\nMultiDonor NSC vs"
PT <- PolioudakisTab
PT$Batch1 <- NULL
PT$Batch2 <- NULL
PT$Within_Donor <- NULL
PT$Within_Batch <- NULL
PT$Donors <- NULL
PT$Batches <- NULL
PT$Within_Protocol <- "Fetal"
PT$Within_Donor <- TRUE
PT$InvolvesMultidonor <- FALSE
PT$OtherProt <- "Polioudakis"
PT$Treatment1 <- "Control"
PT$Treatment2 <- "Control"
PT$Within_Treatment <- TRUE
PT$Control_Only <- TRUE
cT <- cziTab
cT$Batch1 <- NULL
cT$Batch2 <- NULL
cT$Within_Donor <- NULL
cT$Within_Batch <- NULL
cT$Donors <- NULL
cT$Batches <- NULL
cT$Within_Protocol <- "Fetal"
cT$Within_Donor <- TRUE
cT$InvolvesMultidonor <- FALSE
cT$OtherProt <- "CZI"
cT$Treatment1 <- "Control"
cT$Treatment2 <- "Control"
cT$Within_Treatment <- TRUE
cT$Control_Only <- TRUE
TT <- TrevinoTab
TT$Batch1 <- NULL
TT$Batch2 <- NULL
TT$Within_Donor <- NULL
TT$Within_Batch <- NULL
TT$Donors <- NULL
TT$Batches <- NULL
TT$Within_Protocol <- "Fetal"
TT$Within_Donor <- TRUE
TT$InvolvesMultidonor <- FALSE
TT$OtherProt <- "Trevino"
TT$Treatment1 <- "Control"
TT$Treatment2 <- "Control"
TT$Within_Treatment <- TRUE
TT$Control_Only <- TRUE
TT$Donor1 <- "Trevino"
cT$Donor1 <- "CZI"
PT$Donor1 <- "Polioudakis"

fT2 <- rbind(fT2a, PT, cT, TT)
table(fT2$Within_Protocol)
fT2$Within_Protocol <- factor(fT2$Within_Protocol, levels=c("Between-Protocols:\nMultiDonor NSC vs", "Within-Protocol", "Fetal"))
DonorCols["Polioudakis"] <- DonorCols["CZI"] <- DonorCols["Trevino"] <- "cyan"
ggplot(fT2[fT2$Within_Donor & (fT2$InvolvesMultidonor | fT2$Within_Protocol %in% c("Within-Protocol","Fetal")) & fT2$Control_Only,],
       aes(x=OtherProt, y=Distance, fill=Donor1)) +
  geom_boxplot(width=0.4, aes(fill="grey"), outlier.shape = NA) +
  geom_abline(intercept=mean(cziTab$Distance), slope=0, linetype="dashed", color="#696969", linewidth=1) +
  geom_abline(intercept=mean(TrevinoTab$Distance), slope=0, linetype="dashed", color="#3B3B3B", linewidth=1) +
  geom_abline(intercept=mean(PolioudakisTab$Distance), slope=0, linetype="dashed", color="#aeaeae", linewidth=1) +
  geom_jitter(pch=21, width=0.2) +
  facet_grid(~Within_Protocol, space="free", scales="free") +
  scale_fill_manual(values=DonorCols[unique(fT2$Donor1)], name="Donor") +
  theme(panel.background = element_blank()) +
  xlab("")
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aitchison_BetweenProtocolsAndWithinProtocols_WithinDonor_",Sys.Date(),".svg"),
       height=4, width=7)



ggplot(fullTab[fullTab$InvolvesMultidonor | fullTab$Within_Protocol,],
       aes(x=OtherProt, y=Distance)) +
  geom_violin(width=0.5) +
  geom_jitter(width=0.05, alpha=0.4) +
  facet_grid(~Within_Protocol, space="free", scales="free", labeller=labeller(Within_Protocol=myLab)) +
  theme(panel.background = element_blank()) +
  xlab("")
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Plots/Aithison_BetweenProtocolsAndWithinProtocols_DonorAgnostic_",Sys.Date(),".png"),
       height=4, width=6)



VelascoIn <- fullTab$Distance[fullTab$Protocol1=="Velasco" & fullTab$Within_Protocol & fullTab$Within_Donor & fullTab$Control_Only]
NPCIn <- fullTab$Distance[fullTab$Protocol1=="NPC" & fullTab$Within_Protocol & fullTab$Within_Donor & fullTab$Control_Only]
SDIn <- fullTab$Distance[fullTab$Protocol1=="SingleDonorNSC" & fullTab$Within_Protocol & fullTab$Within_Donor & fullTab$Control_Only]
MDIn <- fullTab$Distance[fullTab$Protocol1=="MultiDonorNSC" & fullTab$Within_Protocol & fullTab$Within_Donor & fullTab$Control_Only]
VelascoMD <- fullTab$Distance[fullTab$OtherProt=="Velasco" & fullTab$InvolvesMultidonor & fullTab$Within_Donor & fullTab$Control_Only]
NPCMD <- fullTab$Distance[fullTab$OtherProt=="NPC" & fullTab$InvolvesMultidonor & fullTab$Within_Donor & fullTab$Control_Only]
SDMD <- fullTab$Distance[fullTab$OtherProt=="NSC\nSingle\nDonor" & fullTab$InvolvesMultidonor & fullTab$Within_Donor & fullTab$Control_Only]


var.test(x=VelascoIn, y=VelascoMD, alternative="two.sided")$p.value
t.test(VelascoIn, VelascoMD)$p.value
var.test(MDIn, SDMD, alternative="two.sided")$p.value
t.test(MDIn, SDMD)$p.value
var.test(SDIn, SDMD, alternative="two.sided")$p.value
t.test(SDIn, SDMD)$p.value
var.test(SDMD, NPCMD, alternative="two.sided")$p.value
t.test(SDMD, NPCMD)$p.value
var.test(VelascoIn, SDMD, alternative="two.sided")$p.value
t.test(VelascoIn, SDMD)$p.value
var.test(VelascoIn, VelascoMD, alternative="two.sided")$p.value
t.test(VelascoIn, VelascoMD)$p.value
ThingMat <- matrix(ncol=2, nrow=7)
ThingMat[,1] <- c(mean(VelascoIn), mean(NPCIn), mean(SDIn), mean(MDIn), mean(VelascoMD), mean(NPCMD), mean(SDMD))
ThingMat[,2] <- c(sd(VelascoIn), sd(NPCIn), sd(SDIn), sd(MDIn), sd(VelascoMD), sd(NPCMD), sd(SDMD))
rownames(ThingMat) <- c("Velasco_Internal","NPC_Internal","SD_Internal","MD_Internal","Velasco_MD","NPC_MD","SD_MD")
ThingMat
library(lsr)
cohensD(SDMD, MDIn)
cohensD(SDMD, SDIn)
cohensD(VelascoMD, MDIn)
cohensD(VelascoMD, VelascoIn)
cohensD(SDMD, MDIn)

t.test(VelascoMD, c(VelascoIn,MDIn))
cohensD(VelascoMD, c(VelascoIn,MDIn))

t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)
t.test(VelascoIn, VelascoMD)



t.test(VelascoIn, SDIn)$p.value
var.test(VelascoIn, SDIn, alternative="two.sided")
wilcox.test(VelascoIn, SDIn, alternative="two.sided")$p.value

t.test(VelascoIn, MDIn)$p.value
var.test(VelascoIn, MDIn, alternative="two.sided")
wilcox.test(VelascoIn, MDIn, alternative="two.sided")$p.value

t.test(MDIn, SDIn)$p.value
var.test(MDIn, SDIn, alternative="two.sided")
wilcox.test(MDIn, SDIn, alternative="two.sided")$p.value

t.test(NPCIn, SDIn)$p.value
var.test(NPCIn, SDIn, alternative="two.sided")
wilcox.test(NPCIn, SDIn, alternative="two.sided")

t.test(NPCIn, MDIn)$p.value
var.test(NPCIn, MDIn, alternative="two.sided")
wilcox.test(NPCIn, MDIn, alternative="two.sided")$p.value

t.test(NPCIn, VelascoIn)$p.value
var.test(NPCIn, VelascoIn, alternative="two.sided")
wilcox.test(NPCIn, VelascoIn, alternative="two.sided")$p.value


OnlyThing <- NoEth[grepl("Mito210", rownames(NoEth)) | grepl("PGP1", rownames(NoEth)),
                   grepl("Mito210", colnames(NoEth)) | grepl("PGP1", colnames(NoEth))]
OnlyThing <- OnlyThing[grep("Control", rownames(OnlyThing)), grep("Control",colnames(OnlyThing))]
topAnno2 <- HeatmapAnnotation(Treatment=as.character(sapply(colnames(OnlyThing), function(x){strsplit(x,"_")[[1]][3]})),
                              Donor=as.character(sapply(colnames(OnlyThing), function(x){strsplit(x,"_")[[1]][5]})),
                              Batch=as.character(sapply(colnames(OnlyThing), function(x){strsplit(x,"_")[[1]][1]})),
                              Protocol=as.character(sapply(colnames(OnlyThing), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                              col = list(Treatment=treatCols,
                                         Donor=myCols2,
                                         Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                         Protocol=ProtCols))
sideAnno2 <- rowAnnotation(Treatment=as.character(sapply(rownames(OnlyThing), function(x){strsplit(x,"_")[[1]][3]})),
                           Donor=as.character(sapply(rownames(OnlyThing), function(x){strsplit(x,"_")[[1]][5]})),
                           Batch=as.character(sapply(rownames(OnlyThing), function(x){strsplit(x,"_")[[1]][1]})),
                           Protocol=as.character(sapply(rownames(OnlyThing), function(x){ifelse(grepl("SingleDonor",strsplit(x,"_")[[1]][1]),"Single","Multi")})),
                           col = list(Treatment=treatCols,
                                      Donor=myCols2,
                                      Batch=c("BatchA"="pink","BatchB"="skyblue","SingleDonorA"="darkred","SingleDonorB"="darkblue","BatchC"="orange"),
                                      Protocol=ProtCols))
Heatmap(OnlyThing, top_annotation = topAnno2, left_annotation = sideAnno2,
        row_order = order(grepl("SingleDonor",rownames(OnlyThing)),sapply(rownames(OnlyThing), function(x){strsplit(x,"_")[[1]][5]})),
        column_order = order(grepl("SingleDonor",colnames(OnlyThing)),sapply(colnames(OnlyThing), function(x){strsplit(x,"_")[[1]][5]})),
        show_column_names = FALSE, show_row_names = FALSE)


aTab$Protocols <- "Between Protocols"
aTab$Protocols[aTab$Protocol1==aTab$Protocol2] <- aTab$Protocol1[aTab$Protocol1==aTab$Protocol2]
aTab$Protocols <- factor(aTab$Protocols, levels=c("Single", "Multi", "Between Protocols"))
ggplot(aTab[aTab$Protocol1%in% c("Multi","Single") & aTab$Protocol2 %in% c("Multi","Single") &
              aTab$Within_Donor & aTab$Within_Treatment & aTab$Treat1=="Control",],
       aes(x=Protocols, y=Distance, fill=Donor1)) + 
  geom_boxplot() +
  facet_grid(~Donor1) +
  scale_fill_manual(values=DonorCols[unique(aTab$Donor1)]) +
  theme(axis.text.x=element_text(angle=90),
        panel.background = element_blank())







  
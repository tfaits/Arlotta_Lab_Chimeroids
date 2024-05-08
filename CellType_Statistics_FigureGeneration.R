## Cell type statistics and some figure generation
library(Seurat)
library(ggplot2)
library(plyr)
library(reshape2)
library(ggh4x)
library(lme4)
library(edgeR)

seurMeta <- readRDS("Metadata_NSC_SD-MD_Control-VPA-EtOH.rds")
multi <- seurMeta[seurMeta$Protocol=="MultiDonorNSC",]
mock <- seurMeta[seurMeta$Protocol=="SingleDonorNSC",]
myCols2 <- c("#AE2012","darkred","#0A9396","#0A9396","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012")
names(myCols2) <- c("HUES66","11a","CIRM17","CW","CIRM","GM","H1","PGP1","Mito210","HUES66_AC2")
TypeCols <- c("#8dd3c7", "#bebada", "#fccde5", "#a6d854","#a6d854", "#fb9a99",
                       "#08519c", "#bf812d","#dfc27d",
                       "#d9d9d9","#d9d9d9",'#fff7bc', '#80b1d3', '#fcbba1', "#fdb462", "#fa9fb5",'#fa9fb5','#fa9fb5',
                       '#dd3497', '#d9f0a3' ,'#a1d99b', '#fee391', "#02818a",
                       "#1ef7e9","#1ef7e9", "#04cf18","#04cf18","#04cf18","#b6c6b4",
                       "#d376ab","#c17cc1","#cccc8c","#66a556","#555555")
names(TypeCols) <- c("aRG", "IP", "Cortical hem", "Cajal Retzius", "CR", "Newborn PN",
                                            "Newborn DL PN", "Subcortical progenitors", "Subcortical neurons",
                                            "Unknown","Unknown/Ribo","oRG", "CFuPN", "PN", "CPN", "Immature IN","IN","GABA IN",
                                            "IN progenitors", "Astroglia", "oRG/Astroglia", "oRG II", "Glial precursors",
                                            "Cycling","Cycling Progenitors","Choroid_Plexus","ChPl","ChPL","aRG/PN",
                                            "GABA INP", "Cycling GABA INP","Mito","Dorsal pallium","Low Quality")
treatCols <- c()
treatCols["EtOH"] <- "#B497D4"
treatCols["VPA"] <- "#F2CC54"
treatCols["CNT"] <- "#EAEAEA"
treatCols["Control"] <- treatCols["CNT"]

mock$Replicate <- sapply(mock$Sample, function(x){strsplit(x,"_")[[1]][length(strsplit(x,"_")[[1]])]})
full <- rbind(mock[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb",
                      "seurat_clusters","Batch", "Mix","Genotype","Treatment","Sample","Organoid","Individual",
                      "Protocol","MergedAnnotation","UMAP_1","UMAP_2","WorkingAnnotation")],
              multi[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb",
                       "seurat_clusters","Batch", "Mix","Genotype","Treatment","Sample","Organoid","Individual",
                       "Protocol","MergedAnnotation","UMAP_1","UMAP_2","WorkingAnnotation")])


full$Genotype[full$Genotype=="CIRM17"] <- "CW"
full$GenoFactor <- factor(full$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))

### Figure 1
#### D
ControlMultidonors <- full[full$Treat=="Control" & full$Protocol=="MultiDonorNSC",]
myMat <- table(ControlMultidonors$Organoid, ControlMultidonors$Genotype)
myMat <- myMat/apply(myMat,1,sum)*100
myMat <- melt(myMat)
colnames(myMat) <- c("Sample","Genotype","Percent")
myMat$Sample <- as.character(myMat$Sample)
myMat$Mix <- as.character(sapply(myMat$Sample, function(x){strsplit(x,"_")[[1]][2]}))
myMat$Batch <- as.character(sapply(myMat$Sample, function(x){strsplit(x,"_")[[1]][1]}))
myMat$Replicate <- as.character(sapply(myMat$Sample, function(x){strsplit(x,"_")[[1]][4]}))
myMat$Genotype <- factor(myMat$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))
ggplot(myMat, aes(x=Replicate, y=Percent, fill=Genotype)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=myCols2[levels(myMat$Genotype)]) +
  facet_grid(Mix~Batch, space="free", scales="free") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color="grey90"))

full$WAfactor <- factor(full$MergedAnnotation,
                        levels=c("aRG","oRG","IP","PN","CFuPN","CPN",
                                 "IN","ChPL","CR","Unknown"))

ggplot(full[full$Treatment=="Control" & full$Protocol=="MultiDonorNSC",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  scale_color_manual(values=TypeCols[levels(full$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none") +
  ggtitle(paste0("All donors [",floor(sum(full$Treatment=="Control" & full$Protocol=="MultiDonorNSC")/1000),",",
                 sum(full$Treatment=="Control" & full$Protocol=="MultiDonorNSC")%%1000," cells]"))

ggplot(full[full$Treatment=="Control" & full$Protocol=="MultiDonorNSC",], aes(x=UMAP_1, y=UMAP_2, color=GenoFactor)) +
  geom_point(size=0.15) +
  scale_color_manual(values=myCols2[levels(full$GenoFactor)], name="Donor") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none") +
  ggtitle("Donor distribution")

tmp <- table(paste(multi$Treatment, multi$Mix, multi$Genotype, sep="-"))
multi$NumCells1 <- paste0(multi$Genotype,"\n[",tmp[paste(multi$Treatment, multi$Mix, multi$Genotype, sep="-")],"]")
multi$NumCells1 <- factor(multi$NumCells1,
                          levels=c("CW\n[16077]","CW\n[15897]","CW\n[25299]","CW\n[21647]","CW\n[5813]","CW\n[2705]",
                                   "H1\n[10660]","H1\n[6275]","H1\n[6307],H1\n[7382]","H1\n[6825]","H1\n[4452]",
                                   "Mito210\n[14409]","Mito210\n[11280]","Mito210\n[9441]","Mito210\n[10651]","Mito210\n[13755]","Mito210\n[9835]",
                                   "PGP1\n[14282]","PGP1\n[10230]","PGP1\n[13003]","PGP1\n[11377]","PGP1\n[7720]","PGP1\n[5785]",
                                   "11a\n[4492]","11a\n[6436]","11a\n[3192]"))
ggplot(multi[multi$Treatment=="Control" & multi$Mix=="Mix1",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells1) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

ggplot(multi[multi$Treatment=="Control" & multi$Mix=="Mix2",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells1) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      Qup = quantile(x[[col]], probs=0.75, na.rm=TRUE)[[1]],
      Qdown = quantile(x[[col]], probs=0.25, na.rm=TRUE)[[1]],
      median = median(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum$Percent <- data_sum$mean
  return(data_sum)
}
propTable <- table(multi$MergedAnnotation[multi$Treatment=="Control"], multi$Individual[multi$Treatment=="Control"])
propTable <- t(t(propTable)/apply(propTable,2,sum)*100)
propTable <- melt(propTable)
colnames(propTable) <- c("CellType","PseudoSample","value")
propTable$PseudoSample <- as.character(propTable$PseudoSample)
propTable$Line <- sapply(propTable$PseudoSample, function(x){strsplit(x,"_")[[1]][5]})
propTable$Treatment <- sapply(propTable$PseudoSample, function(x){strsplit(x,"_")[[1]][3]})
propTable$Sample <- sapply(propTable$PseudoSample, function(x){paste(strsplit(x,"_")[[1]][c(1,2,4)], collapse="_")})
propTable$Batch <- sapply(propTable$PseudoSample, function(x){strsplit(x,"_")[[1]][1]})
propTable$Mix <- sapply(propTable$PseudoSample, function(x){strsplit(x,"_")[[1]][2]})
propTable$PseudoTreat <- paste(propTable$Treatment, propTable$Line, sep="-")

pt2 <- data_summary(propTable, varname="value", groupnames=c("Line","CellType"))
pt2$CellType <- factor(pt2$CellType, levels=levels(multi$WAfactor))
pt2$Line[pt2$Line=="CIRM17"] <- "CW"
pt2$Line <- factor(pt2$Line, levels=levels(multi$GenoFactor))

ggplot(pt2, aes(x=CellType, y=median, fill=CellType)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=TypeCols[levels(pt2$CellType)]) +
  facet_nested(~Line, scales="free", space="free")+
  xlab("") + ylab("Percent abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("[Mix1 + Mix2] Cell type distribution by donor") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_blank(),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12),
        axis.ticks.x = element_blank())

propTable$Line[propTable$Line == "CIRM17"] <- "CW"
propTable$Line <- factor(propTable$Line, levels=c("CW","H1","Mito210","PGP1","11a"))
ggplot(pt2, aes(x=CellType, y=median, fill=CellType)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=TypeCols[levels(pt2$CellType)]) +
  facet_nested(~Line, scales="free", space="free")+
  xlab("") + ylab("Percent abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("[Mix1 + Mix2] Cell type distribution by donor") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_blank(),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12),
        axis.ticks.x = element_blank()) +
  geom_jitter(data=propTable, aes(x=CellType, y=value, fill=CellType),
              color="black", pch=21, alpha=0.5, size=1, width=0.25) +
  geom_errorbar(data=pt2, aes(ymin=Qdown, ymax=Qup), width=.2, position=position_dodge(.9))


myMat <- table(multi$GenoFactor[multi$Treatment=="Control"], multi$WAfactor[multi$Treatment=="Control"])
myMat <- myMat/apply(myMat,1,sum)*100
myMat <- melt(myMat)
ggplot(myMat, aes(x=Var1, y=value, fill=Var2)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values=TypeCols[levels(myMat$Var2)], name="Cell type") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color="grey90"))+
  xlab("") + ylab("Percent of cells")


## Let's begin some statistics!
##
tmp <- table(multi$Treatment)
multi$NumCells_Treatment <- paste0(multi$Treatment, "\n[", tmp[multi$Treatment],"]")
multi$NumCells_Treatment <- mapvalues(multi$NumCells_Treatment,
                                      from=c("Control\n[103602]", "EtOH\n[60082]", "VPA\n[111543]"),
                                      to=c("Control\n[103,602]", "EtOH\n[60,082]", "VPA\n[111,543]"))
multi$NumCells_Treatment <- factor(multi$NumCells_Treatment,
                                   levels=c("Control\n[103,602]", "EtOH\n[60,082]", "VPA\n[111,543]"))

ggplot(multi, aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells_Treatment) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

tmp <- table(multi$Organoid)
multi$NumCells_Org <- paste0("Batch.", multi$Batch, " ", sapply(multi$Organoid, function(x){strsplit(x,"_")[[1]][4]}), "\n[",
                             floor(tmp[multi$Organoid]/1000), ",", tmp[multi$Organoid]%%1000,"]")

ggplot(multi[multi$Treatment=="EtOH" & multi$Mix=="Mix1",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells_Org) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

ggplot(multi[multi$Treatment=="EtOH" & multi$Mix=="Mix2",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells_Org) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

ggplot(multi[multi$Treatment=="VPA" & multi$Mix=="Mix1",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells_Org) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

ggplot(multi[multi$Treatment=="VPA" & multi$Mix=="Mix2",], aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  facet_grid(~NumCells_Org) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        axis.text=element_blank(),
        axis.title=element_blank())

## Barplots for Figure 4
propTable <- table(multi$MergedAnnotation, multi$Individual)
propTable <- t(t(propTable)/apply(propTable,2,sum)*100)
propTable <- melt(propTable)
colnames(propTable) <- c("CellType","Sample_Genotype","value")
propTable$Sample <- as.character(sapply(as.character(propTable$Sample_Genotype), function(x){paste(strsplit(x,"_")[[1]][1:4],collapse="_")}))
propTable$Treatment <- sapply(propTable$Sample, function(x){strsplit(x,"_")[[1]][3]})
propTable$Genotype <- as.character(sapply(as.character(propTable$Sample_Genotype), function(x){strsplit(x,"_")[[1]][5]}))
pt2 <- data_summary(propTable, varname="value", groupnames=c("Treatment","CellType"))


## Statistics time! Run a linear model on each dataset (multidonor/singledonor)
myTab <- unclass(table(multi$Individual, multi$MergedAnnotation))
myTab <- data.frame(myTab)
myTab$LibSize <- apply(myTab,1,sum)
myTab$Treatment <- sapply(rownames(myTab), function(x){strsplit(x,"_")[[1]][3]})
myTab$Genotype  <- sapply(rownames(myTab), function(x){strsplit(x,"_")[[1]][5]})
myTab$Sample    <- sapply(rownames(myTab), function(x){paste(strsplit(x,"_")[[1]][1:4],collapse="_")})
myTab$Batch <- sapply(rownames(myTab), function(x){strsplit(x,"_")[[1]][1]})
myTab$Mix <- sapply(rownames(myTab), function(x){strsplit(x,"_")[[1]][2]})
RandMod <- list()
InterMod <- list()
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR")){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="EtOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="EtOH",])
}
vpaP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ
RandMod <- list()
InterMod <- list()
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL")){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
}
etoP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  etoP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
etoQ <- p.adjust(etoP, method="BH")
etoQ

pt2$Significant <- ""
for(i in names(vpaQ)[vpaQ<0.1]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "*"
}
for(i in names(vpaQ)[vpaQ<0.01]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "**"
}
for(i in names(vpaQ)[vpaQ<0.001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "***"
}
for(i in names(vpaQ)[vpaQ<0.0001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "****"
}
for(i in names(etoQ)[etoQ<0.1]){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$CellType==i] <- "*"
}
for(i in names(etoQ)[etoQ<0.01]){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$CellType==i] <- "**"
}
for(i in names(etoQ)[etoQ<0.001]){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$CellType==i] <- "***"
}
for(i in names(etoQ)[etoQ<0.0001]){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$CellType==i] <- "****"
}
pt2$CellType <- as.character(pt2$CellType)
pt2$CellType <- factor(pt2$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN", "ChPL","CR"))
pt2 <- pt2[!is.na(pt2$CellType),]
ggplot(pt2[pt2$CellType != "Unknown",], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

## Single-donor:
myTabm <- unclass(table(mock$Individual, mock$MergedAnnotation))
myTabm <- data.frame(myTabm)
myTabm$LibSize <- apply(myTabm,1,sum)
myTabm$Treatment <- sapply(rownames(myTabm), function(x){strsplit(x,"_")[[1]][3]})
myTabm$Genotype  <- sapply(rownames(myTabm), function(x){strsplit(x,"_")[[1]][5]})
myTabm$Sample    <- sapply(rownames(myTabm), function(x){paste(strsplit(x,"_")[[1]][1:4],collapse="_")})
myTabm$Batch <- sapply(rownames(myTabm), function(x){strsplit(x,"_")[[1]][1]})
RandMod <- list()
InterMod <- list()
for(i in unique(mock$MergedAnnotation)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + (1|Genotype), data=myTabm[myTabm$Treatment!="EtOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + (1|Genotype), data=myTabm[myTabm$Treatment!="EtOH",])
}
vpaPm <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaPm[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQm <- p.adjust(vpaPm, method="BH")
vpaQm
RandMod <- list()
InterMod <- list()
for(i in unique(mock$MergedAnnotation)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + (1|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + (1|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
}
etoPm <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  etoPm[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
etoQm <- p.adjust(etoPm, method="BH")
etoQm

propTableM <- table(mock$MergedAnnotation, mock$Individual)
propTableM <- t(t(propTableM)/apply(propTableM,2,sum)*100)
propTableM <- melt(propTableM)
colnames(propTableM) <- c("CellType","Sample_Genotype","value")
propTableM$Sample <- as.character(sapply(as.character(propTableM$Sample_Genotype), function(x){paste(strsplit(x,"_")[[1]][1:4],collapse="_")}))
propTableM$Treatment <- sapply(propTableM$Sample, function(x){strsplit(x,"_")[[1]][3]})
propTableM$Genotype <- as.character(sapply(as.character(propTableM$Sample_Genotype), function(x){strsplit(x,"_")[[1]][5]}))
ptM <- data_summary(propTableM, varname="value", groupnames=c("Treatment","CellType"))

ptM$Significant <- ""
for(i in names(vpaQm)[vpaQm<0.1]){
  ptM$Significant[ptM$Treatment=="VPA" & ptM$CellType==i] <- "*"
}
for(i in names(vpaQm)[vpaQm<0.01]){
  ptM$Significant[ptM$Treatment=="VPA" & ptM$CellType==i] <- "**"
}
for(i in names(vpaQm)[vpaQm<0.001]){
  ptM$Significant[ptM$Treatment=="VPA" & ptM$CellType==i] <- "***"
}
for(i in names(vpaQm)[vpaQm<0.0001]){
  ptM$Significant[ptM$Treatment=="VPA" & ptM$CellType==i] <- "****"
}
for(i in names(etoQm)[etoQm<0.1]){
  ptM$Significant[ptM$Treatment=="EtOH" & ptM$CellType==i] <- "*"
}
for(i in names(etoQm)[etoQm<0.01]){
  ptM$Significant[ptM$Treatment=="EtOH" & ptM$CellType==i] <- "**"
}
for(i in names(etoQm)[etoQm<0.001]){
  ptM$Significant[ptM$Treatment=="EtOH" & ptM$CellType==i] <- "***"
}
for(i in names(etoQm)[etoQm<0.0001]){
  ptM$Significant[ptM$Treatment=="EtOH" & ptM$CellType==i] <- "****"
}
ptM$CellType <- as.character(ptM$CellType)
ptM$CellType <- factor(ptM$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN", "ChPL","CR"))
ptM <- ptM[!is.na(ptM$CellType),]
ggplot(ptM[ptM$CellType != "Unknown",], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(ptM$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=ptM,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

## Splitting the above barplots up into "high", "mid", and "low" abundance cell types:
ggplot(ptM[ptM$CellType %in% c("aRG","IN","CR"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(ptM$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12))  +
  geom_jitter(data=propTableM[propTableM$CellType %in% c("aRG","IN","CR"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=ptM[ptM$CellType %in% c("aRG","IN","CR"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
#Mid SD:
ggplot(ptM[ptM$CellType %in% c("oRG","IP","ChPL"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(ptM$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12))  +
  geom_jitter(data=propTableM[propTableM$CellType %in% c("oRG","IP","ChPL"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=ptM[ptM$CellType %in% c("oRG","IP","ChPL"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
#High SD:
ggplot(ptM[ptM$CellType %in% c("PN","CFuPN","CPN"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(ptM$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12))  +
  geom_jitter(data=propTableM[propTableM$CellType %in% c("PN","CFuPN","CPN"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=ptM[ptM$CellType %in% c("PN","CFuPN","CPN"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## LowProp MD:
ggplot(pt2[pt2$CellType %in% c("aRG","IN","CR"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("aRG","IN","CR"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=pt2[pt2$CellType %in% c("aRG","IN","CR"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## MD Mid
ggplot(pt2[pt2$CellType %in% c("oRG","IP","ChPL"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("oRG","IP","ChPL"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=pt2[pt2$CellType %in% c("oRG","IP","ChPL"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## MD High
ggplot(pt2[pt2$CellType %in% c("PN","CFuPN","CPN"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("PN","CFuPN","CPN"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
  geom_text(data=pt2[pt2$CellType %in% c("PN","CFuPN","CPN"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

## UMAPs for Figure 5:
tmp <- table(paste(multi$Genotype, multi$Treatment, sep="_"))
multi$NumCells_Donor <- paste0(paste(multi$Genotype, multi$Treatment, sep="_"),
                               "\n[",floor(tmp[paste(multi$Genotype, multi$Treatment, sep="_")]/1000),",",
                               tmp[paste(multi$Genotype, multi$Treatment, sep="_")]%%1000, "]")
multi$NumCells_Donor <- factor(multi$NumCells_Donor,
                               levels=c("CW_Control\n[31,974]","H1_Control\n[16,935]","Mito210_Control\n[25,689]",
                                        "PGP1_Control\n[24,512]","11a_Control\n[4,492]",
                                        "CW_EtOH\n[8,518]","H1_EtOH\n[11,277]","Mito210_EtOH\n[23,590]",
                                        "PGP1_EtOH\n[13,505]","11a_EtOH\n[3,192]",
                                        "CW_VPA\n[46,946]","H1_VPA\n[13,689]","Mito210_VPA\n[20,92]",
                                        "PGP1_VPA\n[24,380]","11a_VPA\n[6,436]"))
ggplot(multi, aes(x=UMAP_1, y=UMAP_2, color=WAfactor)) +
  geom_point(size=0.15) +
  scale_color_manual(values=TypeCols[levels(multi$WAfactor)], name="Cell Type") +
  theme(panel.background = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none") +
  facet_wrap(~NumCells_Donor, ncol=5)


## Barplot showing donor proportions in mix1&mix2?
## Try a joint model with mix as a covariate!
## VPA BothMix
abundances <- table(multi$Genotype[multi$Treatment != "EtOH"], multi$Sample[multi$Treatment!="EtOH"]) 
abundances <- unclass(abundances)
extra.info <- multi[match(colnames(abundances), multi$Sample),c("Sample","Treatment","Mix","Batch")]
vpaWholeChim3 <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Mix)  + factor(Treatment), vpaWholeChim3$samples)
vpaWholeChim3 <- calcNormFactors(vpaWholeChim3, method="TMM")
vpaWholeChim3 <- estimateDisp(vpaWholeChim3, design4, trend="none")
vpaWholeChim3 <- glmQLFit(vpaWholeChim3, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim3 <- glmQLFTest(vpaWholeChim3, coef=ncol(design4))
rvwc1 <- topTags(resVPAWholeChim3)
topTags(resVPAWholeChim3)

## ETOH M1
abundances <- table(multi$Genotype[multi$Treatment != "VPA"], multi$Sample[multi$Treatment!="VPA"]) 
abundances <- unclass(abundances)
extra.info <- multi[match(colnames(abundances), multi$Sample),c("Sample","Treatment","Mix","Batch")]
vpaWholeChim3 <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Mix) + factor(Treatment), vpaWholeChim3$samples)
vpaWholeChim3 <- calcNormFactors(vpaWholeChim3, method="TMM")
vpaWholeChim3 <- estimateDisp(vpaWholeChim3, design4, trend="none")
vpaWholeChim3 <- glmQLFit(vpaWholeChim3, design4, robust=TRUE, abundance.trend=FALSE)
resETOHWholeChim3 <- glmQLFTest(vpaWholeChim3, coef=ncol(design4))
rewc1 <- topTags(resETOHWholeChim3)
topTags(resETOHWholeChim3)

## Plot the results with significance markers
meta <- multi
propTable <- table(meta$Genotype, meta$Sample)
propTable <- t(t(propTable)/apply(propTable,2,sum)*100)
propTable <- melt(propTable)
colnames(propTable) <- c("Genotype","Sample","value")
propTable$Treatment <- sapply(propTable$Sample, function(x){strsplit(as.character(x),"_")[[1]][3]})
propTable$Mix <- as.character(sapply(as.character(propTable$Sample), function(x){strsplit(as.character(x),"_")[[1]][2]}))
propTable$ExpectedValue <- 25
propTable$ExpectedValue[propTable$Mix=="Mix2"] <- 20
propTable$LogFC <- log2(propTable$value) - log2(propTable$ExpectedValue)
propTable <- propTable[!(propTable$Mix=="Mix1" & propTable$Genotype=="11a"),]
propTable$DifferenceFromExpected <- propTable$value - propTable$ExpectedValue
propTable$ProportionalDeviation <- propTable$value/propTable$ExpectedValue
pt2 <- data_summary(propTable, varname="value", groupnames=c("Treatment","Genotype"))
pt2$Significant <- ""
for(i in which(rvwc1$table$FDR<0.1)){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$Genotype==rownames(rvwc1$table)[i]] <- "*"
}
for(i in which(rvwc1$table$FDR<0.01)){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$Genotype==rownames(rvwc1$table)[i]] <- "**"
}
for(i in which(rvwc1$table$FDR<0.001)){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$Genotype==rownames(rvwc1$table)[i]] <- "***"
}
for(i in which(rvwc1$table$FDR<0.0001)){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$Genotype==rownames(rvwc1$table)[i]] <- "****"
}
treatCols["Control"] <- treatCols["CNT"]
#ETOH M1
for(i in which(rewc1$table$FDR<0.1)){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$Genotype==rownames(rewc1$table)[i]] <- "*"
}
for(i in which(rewc1$table$FDR<0.01)){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$Genotype==rownames(rewc1$table)[i]] <- "**"
}
for(i in which(rewc1$table$FDR<0.001)){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$Genotype==rownames(rewc1$table)[i]] <- "***"
}
for(i in which(rewc1$table$FDR<0.0001)){
  pt2$Significant[pt2$Treatment=="EtOH" & pt2$Genotype==rownames(rewc1$table)[i]] <- "****"
}
pt2$Genotype <- as.character(pt2$Genotype)
pt2$Genotype <- factor(pt2$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))
ggplot(pt2, aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~Genotype)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_blank(),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="white"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12),
        axis.ticks.x = element_blank()) +
  geom_text(data=pt2,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

ggplot(pt2, aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_jitter(data=propTable, aes(x=Treatment, y=value, fill=Treatment), alpha=0.5, width=0.25, color="black", pch=21, size=1) +
  geom_errorbar(data=pt2, aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~Genotype)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_blank(),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="white"),
        panel.grid.major.y = element_line(color="grey96"),
        panel.grid.major.x = element_blank(),
        strip.text=element_text(face="bold",size=12),
        axis.ticks.x = element_blank()) +
  geom_text(data=pt2,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")


## Okay, moving on the the larger mixed-effects models that look for...
## ...differences in response to treatment between the donor lines!
#Treatment + Batch + Mix + (Treatment|Genotype) + (1|Sample)
#Treatment + Batch + Mix + (1|Genotype) + (1|Sample)
RandMod <- list()
InterMod <- list()
for(i in unique(multi$MergedAnnotation)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="EtOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (Treatment|Genotype), data=myTab[myTab$Treatment!="EtOH",])
}
vpaP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ[order(vpaQ)]

RandModm <- list()
InterModm <- list()
for(i in unique(mock$MergedAnnotation)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandModm[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + (1|Genotype), data=myTabm[myTabm$Treatment!="EtOH",])
  InterModm[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + (Treatment|Genotype), data=myTabm[myTabm$Treatment!="EtOH",])
}
vpaPm <- c()
for(i in names(RandModm)){
  print(i)
  print(anova(RandModm[[i]], InterModm[[i]])[2,8])
  vpaPm[i] <- anova(RandModm[[i]], InterModm[[i]])[2,8]
}
vpaQm <- p.adjust(vpaPm, method="BH")
vpaQm[order(vpaQm)]

RandMod2 <- list()
InterMod2 <- list()
for(i in unique(multi$MergedAnnotation)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod2[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
  InterMod2[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (Treatment|Genotype), data=myTab[myTab$Treatment!="VPA",])
}
etoP <- c()
for(i in names(InterMod2)){
  print(i)
  print(anova(RandMod2[[i]], InterMod2[[i]])[2,8])
  etoP[i] <- anova(RandMod2[[i]], InterMod2[[i]])[2,8]
}
etoQ <- p.adjust(etoP[c("aRG","oRG","IP","PN","CFuPN","CPN")], method="BH")
etoQ

RandModm <- list()
InterModm <- list()
for(i in unique(mock$CellType)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandModm[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
  InterModm[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (Treatment|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
}
etoPm <- c()
for(i in unique(mock$CellType)){
  print(i)
  print(anova(RandModm[[i]], InterModm[[i]])[2,8])
  etoPm[i] <- anova(RandModm[[i]], InterModm[[i]])[2,8]
}
etoQm <- p.adjust(etoPm, method="BH")
etoQm


#### These results give us just one change that itself changes with donor: IN (in VPA)
## Three others have a raw p-value < 0.1 (and FDR < 0.2): CFuPN, aRG, and CPN
## Let's work on the pairwise-comparisons within these things
myTab$Genotype[myTab$Genotype=="CIRM17"] <- "CW"
myGens <- c("CW","Mito210","H1","PGP1","11a")
PairwiseRes <- list()
for(celltype in c("IN","CFuPN","aRG","CPN")){
  PairwiseRes[[celltype]] <- list()
  for(i in 1:4){
    geni <- myGens[i]
    for(j in (i+1):5){
      genj <- myGens[j]
      PairwiseRes[[celltype]][[paste(geni,genj, sep="-")]] <- list()
      myTab$TypeOfInterest <- myTab[[celltype]]
      myMod1 <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype),
                         data=myTab[myTab$Treatment!="ETOH" & myTab$Genotype %in% c(geni, genj),])
      myMod2 <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (Treatment|Genotype),
                         data=myTab[myTab$Treatment!="ETOH" & myTab$Genotype %in% c(geni, genj),])
      PairwiseRes[[celltype]][[paste(geni,genj, sep="-")]][["Mod1"]] <- myMod1
      PairwiseRes[[celltype]][[paste(geni,genj, sep="-")]][["Mod2"]] <- myMod2
      PairwiseRes[[celltype]][[paste(geni,genj, sep="-")]][["Pval"]] <- anova(myMod1, myMod2)[2,8]
    }
  }
}
PWR <- list()
for(i in names(PairwiseRes)){
  PWR[[i]] <- c()
  for(j in names(PairwiseRes[[i]])){
    PWR[[i]] <- c(PWR[[i]], PairwiseRes[[i]][[j]]$Pval)
  }
}
PWR
PWR_Q <- list()
for(i in names(PWR)){
  PWR_Q[[i]] <- p.adjust(PWR[[i]], method="BH")
}


## Look just at untreated PGP1 samples. Compare across protocols: NSC MD, NSC SD, NPC, Velasco (atlas):
nsc1 <- full
nsc1 <- nsc1[nsc1$Treatment=="Control",]
nsc1$Batch[nsc1$Batch=="SingleDonorA"] <- "NSC_D"
nsc1$Batch[nsc1$Batch=="SingleDonorB"] <- "NSC_C"
nsc1$Batch[nsc1$Batch=="B"] <- "NSC_B"
nsc1$Batch[nsc1$Batch=="C"] <- "NSC_C"
nsc1$CellType <- nsc1$MergedAnnotation
atlas <- readRDS("Metadata_Atlas_3mo.RDS")
atlas <- atlas[atlas$dataset %in% c("4","5","6","7"),]
atlas$Protocol <- "Velasco"
atlas$Batch <- as.character(mapvalues(atlas$dataset, from=c("4","5","6","7"),
                                      to=c("Velasco_A", "Velasco_B","Velasco_C","Velasco_D")))
atlas$CellType <- atlas$FinalName
atlas$Sample <- atlas$org
atlas$Genotype <- as.character(mapvalues(atlas$dataset, from=c("4","5","6","7"),
                                         to=c("Mito210","Mito210","PGP1","PGP1")))

npc <- readRDS("Metadata_M3M4_NPC_NSC.RDS")
npc <- npc[npc$orig.ident %in% c("August2022_Mix1_D23_1","August2022_Mix1_D23_2"),]
npc$Protocol <- "NPC"
npc$Batch <- "NPC_A"
npc$Sample <- as.character(mapvalues(npc$orig.ident,
                                     from=c("August2022_Mix1_D23_1","August2022_Mix1_D23_2"),
                                     to=c("NPC_1","NPC_2")))
npc$CellType <- npc$MergedAnnotation

tmp <- rbind(nsc1[,c("Sample","Batch","CellType", "Protocol","Genotype")],
             atlas[,c("Sample","Batch","CellType", "Protocol","Genotype")],
             npc[,c("Sample","Batch","CellType", "Protocol","Genotype")])
tmp$AllMeta <- paste(tmp$Sample, tmp$Batch, tmp$Protocol, tmp$Genotype, sep="<")

tmp$CellType[tmp$CellType=="Immature IN"] <- "IN"

myTab <- unclass(table(tmp$AllMeta, tmp$CellType))
myTab <- data.frame(myTab)
myTab$LibSize <- apply(myTab,1,sum)
myTab$Sample <- sapply(rownames(myTab), function(x){strsplit(x,"<")[[1]][1]})
myTab$Batch  <- sapply(rownames(myTab), function(x){strsplit(x,"<")[[1]][2]})
myTab$Protocol <- sapply(rownames(myTab), function(x){strsplit(x,"<")[[1]][3]})
myTab$Genotype <- sapply(rownames(myTab), function(x){strsplit(x,"<")[[1]][4]})
myTab$ID <- rownames(myTab)

## We are most interested in comparing each other protocol to the NSC MD. Do sequentially:
## We don't really need a MEM here. We can simply use edgeR's negative-binomial model with library-size as an offset
abundances <- t(myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),11:16]
PGP1_SDMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designPGP <- model.matrix(~factor(Batch) + factor(Protocol), PGP1_SDMD$samples)
PGP1_SDMD <- calcNormFactors(PGP1_SDMD, method="TMM")
PGP1_SDMD <- estimateDisp(PGP1_SDMD, designPGP, trend="none")
PGP1_SDMD <- glmQLFit(PGP1_SDMD, designPGP, robust=TRUE, abundance.trend=FALSE)
resPGP_sdmd <- glmQLFTest(PGP1_SDMD, coef=ncol(designPGP))
rpsdmd <- topTags(resPGP_sdmd)
topTags(resPGP_sdmd)

abundances <- t(myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),11:16]
Mito210_SDMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designMito <- model.matrix(~factor(Batch) + factor(Protocol), Mito210_SDMD$samples)
Mito210_SDMD <- calcNormFactors(Mito210_SDMD, method="TMM")
Mito210_SDMD <- estimateDisp(Mito210_SDMD, designMito, trend="none")
Mito210_SDMD <- glmQLFit(Mito210_SDMD, designMito, robust=TRUE, abundance.trend=FALSE)
resMito_sdmd <- glmQLFTest(Mito210_SDMD, coef=ncol(designMito))
rmsdmd <- topTags(resMito_sdmd)
topTags(resMito_sdmd)

abundances <- t(myTab[myTab$Genotype=="H1" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="H1" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),11:16]
H1_SDMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designH <- model.matrix(~factor(Batch) + factor(Protocol), H1_SDMD$samples)
H1_SDMD <- calcNormFactors(H1_SDMD, method="TMM")
H1_SDMD <- estimateDisp(H1_SDMD, designH, trend="none")
H1_SDMD <- glmQLFit(H1_SDMD, designH, robust=TRUE, abundance.trend=FALSE)
resH_sdmd <- glmQLFTest(H1_SDMD, coef=ncol(designH))
rhsdmd <- topTags(resH_sdmd)
topTags(resH_sdmd)

abundances <- t(myTab[myTab$Genotype=="CW" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="CW" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),11:16]
CW1_SDMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designCW <- model.matrix(~factor(Batch) + factor(Protocol), CW1_SDMD$samples)
CW1_SDMD <- calcNormFactors(CW1_SDMD, method="TMM")
CW1_SDMD <- estimateDisp(CW1_SDMD, designCW)
CW1_SDMD <- glmQLFit(CW1_SDMD, designCW, robust=TRUE, abundance.trend=FALSE)
resCW_sdmd <- glmQLFTest(CW1_SDMD, coef=ncol(designCW))
rcsdmd <- topTags(resCW_sdmd)
topTags(resCW_sdmd)

abundances <- t(myTab[myTab$Genotype=="11a" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="11a" & myTab$Protocol %in% c("SingleDonorNSC","MultiDonorNSC"),11:16]
a11_SDMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design11a <- model.matrix(~factor(Batch) + factor(Protocol), a11_SDMD$samples)
a11_SDMD <- calcNormFactors(a11_SDMD, method="TMM")
a11_SDMD <- estimateDisp(a11_SDMD, design11a)
a11_SDMD <- glmQLFit(a11_SDMD, design11a, robust=TRUE, abundance.trend=FALSE)
res11_sdmd <- glmQLFTest(a11_SDMD, coef=ncol(design11a))
r1sdmd <- topTags(res11_sdmd)
topTags(res11_sdmd)

abundances <- t(myTab[myTab$Genotype=="H1" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="H1" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),11:16]
H1_NPCMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designH <- model.matrix(~factor(Protocol), H1_NPCMD$samples)
H1_NPCMD <- calcNormFactors(H1_NPCMD, method="TMM")
H1_NPCMD <- estimateDisp(H1_NPCMD, designH)
H1_NPCMD <- glmQLFit(H1_NPCMD, designH, robust=TRUE, abundance.trend=FALSE)
resH_NPCmd <- glmQLFTest(H1_NPCMD, coef=ncol(designH))
rhNPCmd <- topTags(resH_NPCmd)
topTags(resH_NPCmd)

abundances <- t(myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("Velasco","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("Velasco","MultiDonorNSC"),11:16]
PGP1_VelascoMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designPGP <- model.matrix(~factor(Protocol), PGP1_VelascoMD$samples)
PGP1_VelascoMD <- calcNormFactors(PGP1_VelascoMD, method="TMM")
PGP1_VelascoMD <- estimateDisp(PGP1_VelascoMD, designPGP, trend="none")
PGP1_VelascoMD <- glmQLFit(PGP1_VelascoMD, designPGP, robust=TRUE, abundance.trend=FALSE)
resPGP_Velascomd <- glmQLFTest(PGP1_VelascoMD, coef=ncol(designPGP))
rpVelascomd <- topTags(resPGP_Velascomd)
topTags(resPGP_Velascomd)

abundances <- t(myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("Velasco","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("Velasco","MultiDonorNSC"),11:16]
Mito210_VelascoMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designMito <- model.matrix(~factor(Protocol), Mito210_VelascoMD$samples)
Mito210_VelascoMD <- calcNormFactors(Mito210_VelascoMD, method="TMM")
Mito210_VelascoMD <- estimateDisp(Mito210_VelascoMD, designMito, trend="none")
Mito210_VelascoMD <- glmQLFit(Mito210_VelascoMD, designMito, robust=TRUE, abundance.trend=FALSE)
resMito_Velascomd <- glmQLFTest(Mito210_VelascoMD, coef=ncol(designMito))
rmVelascomd <- topTags(resMito_Velascomd)
topTags(resMito_Velascomd)

abundances <- t(myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="PGP1" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),11:16]
PGP1_NPCMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designPGP <- model.matrix(~factor(Protocol), PGP1_NPCMD$samples)
PGP1_NPCMD <- calcNormFactors(PGP1_NPCMD, method="TMM")
PGP1_NPCMD <- estimateDisp(PGP1_NPCMD, designPGP, trend="none")
PGP1_NPCMD <- glmQLFit(PGP1_NPCMD, designPGP, robust=TRUE, abundance.trend=FALSE)
resPGP_NPCmd <- glmQLFTest(PGP1_NPCMD, coef=ncol(designPGP))
rpNPCmd <- topTags(resPGP_NPCmd)
topTags(resPGP_NPCmd)

abundances <- t(myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),1:10])
abundances <- unclass(abundances)
extra.info <- myTab[myTab$Genotype=="Mito210" & myTab$Protocol %in% c("NPC","MultiDonorNSC"),11:16]
Mito210_NPCMD <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
designMito <- model.matrix(~factor(Protocol), Mito210_NPCMD$samples)
Mito210_NPCMD <- calcNormFactors(Mito210_NPCMD, method="TMM")
Mito210_NPCMD <- estimateDisp(Mito210_NPCMD, designMito, trend="none")
Mito210_NPCMD <- glmQLFit(Mito210_NPCMD, designMito, robust=TRUE, abundance.trend=FALSE)
resMito_NPCmd <- glmQLFTest(Mito210_NPCMD, coef=ncol(designMito))
rmNPCmd <- topTags(resMito_NPCmd)
topTags(resMito_NPCmd)

abundances <- myTab[,1:10]
abundances <- 100*abundances/apply(abundances,1,sum)

daMat <- data.frame(Protocol="Other", CellType="Perfect", Genotype="HW", Percent=100, Pval=1)
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown")){
  daMat <- rbind(daMat, c("MultiDonor", i, "PGP1", mean(abundances[myTab$Genotype=="PGP1" & myTab$Protocol=="MultiDonorNSC",i]), 1))
  daMat <- rbind(daMat, c("MultiDonor", i, "H1", mean(abundances[myTab$Genotype=="H1" & myTab$Protocol=="MultiDonorNSC",i]), 1))
  daMat <- rbind(daMat, c("MultiDonor", i, "CW", mean(abundances[myTab$Genotype=="CW" & myTab$Protocol=="MultiDonorNSC",i]), 1))
  daMat <- rbind(daMat, c("MultiDonor", i, "11a", mean(abundances[myTab$Genotype=="11a" & myTab$Protocol=="MultiDonorNSC",i]), 1))
  daMat <- rbind(daMat, c("SingleDonor", i, "PGP1", mean(abundances[myTab$Genotype=="PGP1" & myTab$Protocol=="SingleDonorNSC",i]), as.data.frame(rpsdmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("Velasco", i, "PGP1", mean(abundances[myTab$Genotype=="PGP1" & myTab$Protocol=="Velasco",i]), as.data.frame(rpVelascomd)[i,"FDR"]))
  daMat <- rbind(daMat, c("NPC", i, "PGP1", mean(abundances[myTab$Genotype=="PGP1" & myTab$Protocol=="NPC",i]), as.data.frame(rpNPCmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("MultiDonor", i, "Mito210", mean(abundances[myTab$Genotype=="Mito210" & myTab$Protocol=="MultiDonorNSC",i]), 1))
  daMat <- rbind(daMat, c("SingleDonor", i, "Mito210", mean(abundances[myTab$Genotype=="Mito210" & myTab$Protocol=="SingleDonorNSC",i]), as.data.frame(rmsdmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("Velasco", i, "Mito210", mean(abundances[myTab$Genotype=="Mito210" & myTab$Protocol=="Velasco",i]), as.data.frame(rmVelascomd)[i,"FDR"]))
  daMat <- rbind(daMat, c("NPC", i, "Mito210", mean(abundances[myTab$Genotype=="Mito210" & myTab$Protocol=="NPC",i]), as.data.frame(rmNPCmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("SingleDonor", i, "H1", mean(abundances[myTab$Genotype=="H1" & myTab$Protocol=="SingleDonorNSC",i]), as.data.frame(rhsdmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("SingleDonor", i, "CW", mean(abundances[myTab$Genotype=="CW" & myTab$Protocol=="SingleDonorNSC",i]), as.data.frame(rcsdmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("SingleDonor", i, "11a", mean(abundances[myTab$Genotype=="11a" & myTab$Protocol=="SingleDonorNSC",i]), as.data.frame(r1sdmd)[i,"FDR"]))
  daMat <- rbind(daMat, c("NPC", i, "H1", mean(abundances[myTab$Genotype=="H1" & myTab$Protocol=="NPC",i]), as.data.frame(rhNPCmd)[i,"FDR"]))
}
daMat$Pval <- as.numeric(daMat$Pval)
daMat$Percent <- as.numeric(daMat$Percent)
daMat <- daMat[-1,]
daMat$Sig <- ""
daMat$Sig[daMat$Pval < 0.05] <- paste(round(daMat$Pval[daMat$Pval < 0.05], digits=3), "*")
daMat$Sig[daMat$Pval < 0.005] <- paste(round(daMat$Pval[daMat$Pval < 0.005], digits=4), "**")
daMat$Sig[daMat$Pval < 0.0005] <- paste(round(daMat$Pval[daMat$Pval < 0.0005], digits=5), "***")
daMat$Protocol[daMat$Protocol=="MultiDonor"] <- "MD\nNSC"
daMat$Protocol[daMat$Protocol=="SingleDonor"] <- "SD\nNSC"
daMat$Protocol <- factor(daMat$Protocol, levels=c("MD\nNSC","SD\nNSC","NPC","Velasco"))
daMat$CellType <- factor(daMat$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"))
daMat$Genotype <- factor(daMat$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))
ggplot(daMat, aes(x=Protocol, y=Percent, fill=CellType, label=Sig)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=TypeCols[levels(daMat$CellType)]) +
  facet_grid(~Genotype, space="free", scales="free") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme(panel.background=element_blank())

daMat$Sig[daMat$Percent < 1] <- ""
ggplot(daMat, aes(x=Protocol, y=Percent, fill=CellType, label=Sig)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=TypeCols[levels(daMat$CellType)]) +
  facet_grid(~Genotype, space="free", scales="free") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme(panel.background=element_blank())

RandMod <- list()
InterMod <- list()
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown")){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Batch), data=myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC"),])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Protocol + offset(log(LibSize)) + (1|Batch), data=myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC"),])
}
vpaP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ

## We are most interested in comparing each other protocol to the NSC MD. Do sequentially:
RandMod <- list()
InterMod <- list()
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown")){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Batch), data=myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC") & myTab$Genotype=="Mito210",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Protocol + offset(log(LibSize)) + (1|Batch), data=myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC") & myTab$Genotype=="Mito210",])
}
vpaP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ


mT <- table(paste(tmp$Protocol, tmp$Genotype, sep="<"), tmp$CellType)
mT <- mT/apply(mT,1,sum)*100
mT <- melt(mT)
colnames(mT) <- c("Group","CellType","Percent")
mT$Group <- as.character(mT$Group)
mT$Protocol <- sapply(mT$Group, function(x){strsplit(x,"<")[[1]][1]})
mT$Genotype <- sapply(mT$Group, function(x){strsplit(x,"<")[[1]][2]})
mT$CellType <- factor(mT$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"))
mT$Protocol <- factor(mT$Protocol, levels=c("MultiDonorNSC","SingleDonorNSC","NPC","Velasco"))
ggplot(mT, aes(x=Protocol, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=TypeCols[levels(mT$CellType)]) +
  facet_grid(~Genotype) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


mT2 <- melt(cbind(myTab[,1:10]/apply(myTab[,1:10],1,sum)*100, ID=myTab$ID))
colnames(mT2) <- c("Sample","CellType","Percent")
mT2$Sample <- as.character(mT2$Sample)
mT2$Protocol <- sapply(mT2$Sample, function(x){strsplit(x,"<")[[1]][3]})
mT2$Genotype <- sapply(mT2$Sample, function(x){strsplit(x,"<")[[1]][4]})
mT2$Batch <- sapply(mT2$Sample, function(x){strsplit(x,"<")[[1]][2]})
mT2$Sample <- sapply(mT2$Sample, function(x){strsplit(x,"<")[[1]][1]})
mT2$CellType <- factor(mT2$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"))
mT2$Protocol <- factor(mT2$Protocol, levels=c("MultiDonorNSC","SingleDonorNSC","NPC","Velasco"))
pgp1 <- ggplot(mT2[mT2$Genotype=="PGP1",], aes(x=Sample, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=TypeCols[levels(mT$CellType)]) +
  facet_grid(~Protocol, space="free", scales="free") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=0.2),
        axis.title.x=element_blank(),
        panel.background = element_rect(color="black",fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="grey70", linewidth=0.2)) +
  ggtitle("PGP1")
mito210 <- ggplot(mT2[mT2$Genotype=="Mito210",], aes(x=Sample, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=TypeCols[levels(mT$CellType)]) +
  facet_grid(~Protocol, space="free", scales="free") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=0.2),
        axis.title.x=element_blank(),
        panel.background = element_rect(color="black",fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="grey70", linewidth=0.2)) +
  ggtitle("Mito210")
cowplot::plot_grid(pgp1, mito210, ncol=1)

ggplot(mT2[mT2$Genotype=="PGP1",], aes(x=Sample, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=TypeCols[levels(mT$CellType)]) +
  facet_nested(~Protocol + Batch, space="free", scales="free") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=0.2),
        axis.title.x=element_blank(),
        panel.background = element_rect(color="black",fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="grey70", linewidth=0.2)) +
  ggtitle("PGP1")



abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC") & myTab$Genotype == "Mito210",1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)

abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC") & myTab$Genotype == "PGP1",1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)

abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","SingleDonorNSC"),1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Genotype) + factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)


abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","Velasco") & myTab$Genotype == "Mito210",1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)

abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","Velasco") & myTab$Genotype == "PGP1",1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)

abundances <- t(myTab[myTab$Protocol %in% c("MultiDonorNSC","Velasco"),1:10])
head(abundances)
extra.info <- myTab[match(colnames(abundances), rownames(myTab)),c("Sample","Protocol","Batch","Genotype")]
vpaWholeChim <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Genotype) + factor(Protocol), vpaWholeChim$samples)
vpaWholeChim <- calcNormFactors(vpaWholeChim, method="TMM")
vpaWholeChim <- estimateDisp(vpaWholeChim, design4, trend="none")
fit.vpaWholeChim <- glmQLFit(vpaWholeChim, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim <- glmQLFTest(fit.vpaWholeChim, coef=ncol(design4))
rvwc <- topTags(resVPAWholeChim)
topTags(resVPAWholeChim)



                             
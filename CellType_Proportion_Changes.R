require(plyr)
require(lme4)

treatCols <- c()
treatCols["ETOH"] <- "#B497D4"
treatCols["VPA"] <- "#F2CC54"
treatCols["Control"] <- "#EAEAEA"
treatCols["CNT"] <- "#EAEAEA"
`%ni%` = Negate(`%in%`)

############################################
####
####
############################################
#### Run some stats
## Let's look at changes in cell type distribution *by donor*
## Ok, a big one: Look at cell-type-by-line distribution changes.
## Essentially, do the different lines have different changes in cell type distributions?
seur <- readRDS("Completed_Metadata_All_Annotation.RDS")
seur$CellType_Full[seur$seurat_clusters=="25"] <- "ChPL"
seur <- seur[seur$Genotype!="DBL" & seur$CellType_Full != "Low Quality",]
seur$CellType <- as.character(seur$CellType_Full)
seur$Batch[seur$Batch == "B6"] <- "B7"
seur$Genotype[seur$Genotype=="CIRM17"] <- "CW"
seur$CellType[seur$CellType=="Immature IN"] <- "Immature_IN"
seur$CellType[seur$CellType=="Cajal Retzius"] <- "Cajal_Retzius"
seur$Sample_Genotype <- paste(seur$Sample, seur$Genotype, sep="-")
myTab <- unclass(table(seur$Sample_Genotype, seur$CellType))
myTab <- data.frame(myTab)
myTab$LibSize <- apply(myTab,1,sum)
myTab$Treatment <- sapply(rownames(myTab), function(x){strsplit(x,"_")[[1]][1]})
myTab$Genotype  <- sapply(rownames(myTab), function(x){strsplit(x,"-")[[1]][2]})
myTab$Sample    <- sapply(rownames(myTab), function(x){strsplit(x,"-")[[1]][1]})
myTab$Batch <- seur[match(myTab$Sample, seur$Sample),"Batch"]
myTab$Mix <- seur[match(myTab$Sample, seur$Sample),"Mix"]
myTab$Batch[myTab$Treatment=="CNT" & myTab$Batch=="B2"] <- "B1"

## Load in the mock data, treat similarly:
mock <- readRDS("Multi_Merged_Metadata.RDS")
mock$CellType <- as.character(mapvalues(mock$CellType, from=c("Immature IN", "Cajal Retzius"), to=c("IN","CR")))
myTabm <- unclass(table(mock$Sample_Name, mock$CellType))
myTabm <- data.frame(myTabm)
myTabm$LibSize <- apply(myTabm,1,sum)
myTabm$Treatment <- sapply(rownames(myTabm), function(x){strsplit(x," ")[[1]][1]})
myTabm$Genotype  <- sapply(rownames(myTabm), function(x){strsplit(x," ")[[1]][2]})
myTabm$Sample    <- rownames(myTabm)

## before we get into random effects, check a mixed model to see how it compares to the edgeR results
RandMod <- list()
InterMod <- list()
for(i in unique(seur$CellType)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="ETOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="ETOH",])
}
vpaP <- c()
for(i in unique(seur$CellType)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ
RandMod <- list()
InterMod <- list()
for(i in unique(seur$CellType)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
}
vpaP <- c()
for(i in unique(seur$CellType)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
etoQ <- p.adjust(vpaP, method="BH")
etoQ
## Make some plots of the Chimeroid Random Effects results (not looking at random slopes yet - just checking treatment effect)
seur$Sample_Genotype <- paste(seur$Sample, seur$Genotype, sep="-")
meta <- seur
propTable <- table(meta$CellType, meta$Sample_Genotype)
propTable <- t(t(propTable)/apply(propTable,2,sum)*100)
propTable <- melt(propTable)
colnames(propTable) <- c("CellType","Sample_Genotype","value")
propTable$Sample <- as.character(sapply(as.character(propTable$Sample), function(x){strsplit(x,"-")[[1]][1]}))
propTable$Treatment <- sapply(propTable$Sample, function(x){strsplit(x,"_")[[1]][1]})
propTable$Genotype <- as.character(sapply(as.character(propTable$Sample_Genotype), function(x){strsplit(x,"-")[[1]][2]}))
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      Qup = quantile(x[[col]], 0.75, na.rm=TRUE)[[1]],
      Qdown = quantile(x[[col]], 0.25, na.rm=TRUE)[[1]],
      median = median(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum$Percent <- data_sum$mean
  return(data_sum)
}
pt2 <- data_summary(propTable, varname="value", groupnames=c("Treatment","CellType"))
pt2$Significant <- ""
pt2$CellType2 <- sapply(pt2$CellType, function(x){paste(strsplit(as.character(x)," ")[[1]],collapse="_")})
for(i in names(vpaQ)[vpaQ<0.1]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType2==i] <- "*"
}
for(i in names(vpaQ)[vpaQ<0.01]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType2==i] <- "**"
}
for(i in names(vpaQ)[vpaQ<0.001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType2==i] <- "***"
}
for(i in names(vpaQ)[vpaQ<0.0001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType2==i] <- "****"
}
for(i in names(etoQ)[etoQ<0.1]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType2==i] <- "*"
}
for(i in names(etoQ)[etoQ<0.01]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType2==i] <- "**"
}
for(i in names(etoQ)[etoQ<0.001]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType2==i] <- "***"
}
for(i in names(etoQ)[etoQ<0.0001]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType2==i] <- "****"
}
pt2$CellType <- as.character(pt2$CellType)
pt2$CellType[pt2$CellType=="Immature_IN"] <- "IN"
pt2$CellType[pt2$CellType=="Cajal_Retzius"] <- "CR"
pt2$CellType <- factor(pt2$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN", "ChPL","CR","Unknown"))
ggplot(pt2, aes(x=Treatment, y=median, fill=Treatment)) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
ggsave("/Users/tfaits/Documents/Chimera/Figures/CellType_Proportions_Whole_Barplot_0531_RandomNegBin.png", width=12,height=4)
ggsave("/Users/tfaits/Documents/Chimera/Figures/CellType_Proportions_Whole_Barplot_0531_RandomNegBin.svg", width=12,height=4)
## High-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),], aes(x=Treatment, y=median, fill=Treatment)) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## Middle-proportion (for splitting up the scale):
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("oRG","IP","ChPL"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## Low-proportion (for splitting up the scale):
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("aRG","IN","CR"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

## Add points to these plots
## High-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("CPN","PN","CFuPN"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 6, color = "black", fontface = "bold")
## Middle-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("oRG","IP","ChPL"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("oRG","IP","ChPL"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("oRG","IP","ChPL"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## Low-proportion (for splitting up the scale):
propTable$CellType <- as.character(propTable$CellType)
propTable$CellType[propTable$CellType=="Immature_IN"] <- "IN"
propTable$CellType[propTable$CellType=="Cajal_Retzius"] <- "CR"
ggplot(pt2[pt2$CellType %in% c("aRG","IN","CR"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("aRG","IN","CR"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("aRG","IN","CR"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

## before we get into random effects, check a mixed model for the Single-Donor Mocks:
RandMod <- list()
InterMod <- list()
for(i in unique(mock$CellType)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="ETOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="ETOH",])
}
vpaP <- c()
for(i in unique(mock$CellType)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQM <- p.adjust(vpaP, method="BH")
vpaQM
RandMod <- list()
InterMod <- list()
for(i in unique(mock$CellType)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="VPA",])
}
vpaP <- c()
for(i in unique(mock$CellType)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
etoQM <- p.adjust(vpaP, method="BH")
etoQM
## Plot the mocks with the random effects results:
mock$Sample_Genotype <- mock$Sample_Name
meta <- mock
propTable <- table(meta$CellType, meta$Sample_Name)
propTable <- t(t(propTable)/apply(propTable,2,sum)*100)
propTable <- melt(propTable)
colnames(propTable) <- c("CellType","Sample_Genotype","value")
propTable$Sample <- as.character(propTable$Sample_Genotype)
propTable$Treatment <- sapply(propTable$Sample, function(x){strsplit(x," ")[[1]][1]})
propTable$Genotype <- as.character(sapply(as.character(propTable$Sample_Genotype), function(x){strsplit(x," ")[[1]][2]}))
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      Qup = quantile(x[[col]], 0.75, na.rm=TRUE)[[1]],
      Qdown = quantile(x[[col]], 0.25, na.rm=TRUE)[[1]],
      median = median(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum$Percent <- data_sum$mean
  return(data_sum)
}
pt2 <- data_summary(propTable, varname="value", groupnames=c("Treatment","CellType"))
pt2$Significant <- ""
pt2$CellType2 <- sapply(pt2$CellType, function(x){paste(strsplit(as.character(x)," ")[[1]],collapse="_")})
for(i in names(vpaQM)[vpaQM<0.1]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "*"
}
for(i in names(vpaQM)[vpaQM<0.01]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "**"
}
for(i in names(vpaQM)[vpaQM<0.001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "***"
}
for(i in names(vpaQM)[vpaQM<0.0001]){
  pt2$Significant[pt2$Treatment=="VPA" & pt2$CellType==i] <- "****"
}
for(i in names(etoQM)[etoQM<0.1]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType==i] <- "*"
}
for(i in names(etoQM)[etoQM<0.01]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType==i] <- "**"
}
for(i in names(etoQM)[etoQM<0.001]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType==i] <- "***"
}
for(i in names(etoQM)[etoQM<0.0001]){
  pt2$Significant[pt2$Treatment=="ETOH" & pt2$CellType==i] <- "****"
}
treatCols[["Control"]] <- "#EAEAEA"
pt2$CellType <- as.character(pt2$CellType)
pt2$CellType <- factor(pt2$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN", "ChPL","CR","Unknown"))
ggplot(pt2, aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=treatCols[unique(pt2$Treatment)]) +
  facet_grid(~CellType)+
  xlab("") + ylab("% of Sequenced Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("(Single-Donor Mocks) Cell Type Proportions") +
  theme(axis.text=element_text(face="bold",size=14),
        axis.text.x=element_text(angle=60),
        title = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.background = element_rect(fill="white", color="black"),
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2,aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## High-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("CPN","PN","CFuPN"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("CPN","PN","CFuPN"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## Middle-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("oRG","IP","ChPL"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("oRG","IP","ChPL"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("oRG","IP","ChPL"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")
## Low-proportion (for splitting up the scale):
ggplot(pt2[pt2$CellType %in% c("aRG","IN","CR"),], aes(x=Treatment, y=median, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=propTable[propTable$CellType %in% c("aRG","IN","CR"),],
              aes(y=value), pch=21, alpha=0.5, width=0.25, size=1) +
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
        panel.grid.major = element_line(color="grey96"),
        strip.text=element_text(face="bold",size=12)) +
  geom_text(data=pt2[pt2$CellType %in% c("aRG","IN","CR"),],aes(x=Treatment, y=Qup+1, label=Significant), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

##
## While I've got myTab loaded in, let's make some stacked-barplots for each donor within the Multidonor Chimeroids
## Collapse the replicates
myTab2 <- myTab
myTab2[,1:10] <- myTab2[,1:10]/myTab2$LibSize*100
myTab2 <- myTab2[myTab2$Treatment=="CNT",]
myMelt <- melt(myTab2[,c(1:10, 13)])
colnames(myMelt) <- c("Genotype","CellType","Percent")
# Divide by number of reps (for each line) so we'll plot the average values:
myMelt$Percent[myMelt$Genotype=="11a"] <- myMelt$Percent[myMelt$Genotype=="11a"]/3
myMelt$Percent[myMelt$Genotype!="11a"] <- myMelt$Percent[myMelt$Genotype!="11a"]/7
myMelt$Genotype <- factor(myMelt$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))
myMelt$CellType <- factor(myMelt$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"))
myMelt2 <- data.frame(Genotype=c(rep("CW",10),rep("H1",10),rep("Mito210",10),rep("PGP1",10),rep("11a",10)),
                      CellType=rep(c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"),5),
                      Percent=0)
for(i in c("CW","H1","Mito210","PGP1","11a")){
  for(j in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown")){
    myMelt2$Percent[myMelt2$Genotype==i & myMelt2$CellType==j] <- sum(myMelt$Percent[myMelt$Genotype==i & myMelt$CellType==j])
  }
}
myMelt2$Genotype <- factor(myMelt2$Genotype, levels=c("CW","H1","Mito210","PGP1","11a"))
myMelt2$CellType <- factor(myMelt2$CellType, levels=c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR","Unknown"))
ggplot(myMelt2, aes(x=Genotype, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=TypeCols[levels(myMelt$CellType)]) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        panel.background = element_rect(fill="white",color="grey85"))+
  ggtitle("Cell Type Stacked Barplot for Mix1+Mix2\nGrouped By Donor")

#Treatment + Batch + Mix + (Treatment|Genotype) + (1|Sample)
#Treatment + Batch + Mix + (1|Genotype) + (1|Sample)
RandMod <- list()
InterMod <- list()
for(i in unique(seur$CellType)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (1|Genotype), data=myTab[myTab$Treatment!="ETOH",])
  InterMod[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Sample) + (Treatment|Genotype), data=myTab[myTab$Treatment!="ETOH",])
}
vpaP <- c()
for(i in unique(seur$CellType)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ[order(vpaQ)]

RandModm <- list()
InterModm <- list()
for(i in unique(mock$CellType)){
  myTabm$TypeOfInterest <- myTabm[[i]]
  RandModm[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data=myTabm[myTabm$Treatment!="ETOH",])
  InterModm[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (Treatment|Genotype), data=myTabm[myTabm$Treatment!="ETOH",])
}
vpaPm <- c()
for(i in unique(mock$CellType)){
  print(i)
  print(anova(RandModm[[i]], InterModm[[i]])[2,8])
  vpaPm[i] <- anova(RandModm[[i]], InterModm[[i]])[2,8]
}
vpaQm <- p.adjust(vpaPm, method="BH")
vpaQm

RandMod2 <- list()
InterMod2 <- list()
for(i in unique(seur$CellType)){
  myTab$TypeOfInterest <- myTab[[i]]
  RandMod2[[i]]  <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (1|Genotype), data=myTab[myTab$Treatment!="VPA",])
  InterMod2[[i]] <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + Batch + Mix + (Treatment|Genotype), data=myTab[myTab$Treatment!="VPA",])
}
etoP <- c()
for(i in unique(seur$CellType)){
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

## These results give us just three changes that themselves change with donor: oRGs, IPs, and CPNs (in VPA)
## Let's work on the pairwise-comparisons within these things
myTab$Genotype[myTab$Genotype=="CIRM17"] <- "CW"
myGens <- c("CW","Mito210","H1","PGP1","11a")
PairwiseRes <- list()
for(celltype in c("oRG","IP","CPN")){
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
for(i in c("oRG","IP","CPN")){
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

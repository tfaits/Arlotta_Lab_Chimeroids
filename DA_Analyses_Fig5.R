## Some new Differential Cell Type Abundance (DA) analyses for Figure 5?
library(ggplot2)
library(reshape2)
library(lme4)
library(edgeR)

meta <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/Metadata_AllMerged_Streamlined_November17_2023.RDS")

## First, let's pull out IN proportions, and separate each donor?
## eg: look for how much INs are induced in 11a, compare that change between multidonor and single donor?

## First actual step is to create a long dataframe with counts for each cell type, including library size (for DA offset)
fullTab <- table(meta$Individual, meta$MergedAnnotation)
rownames(fullTab) <- paste(rownames(fullTab), apply(fullTab,1,sum), sep="_")
myTab <- melt(fullTab)
colnames(myTab) <- c("Individual","CellType","Percent")
myTab$Individual <- as.character(myTab$Individual)
myTab$Batch <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][1]}))
myTab$Mix <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][2]}))
myTab$Treatment <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][3]}))
myTab$Replicate <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][4]}))
myTab$Donor <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][5]}))
myTab$LibSize <- as.numeric(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][6]}))

fTab <- data.frame(Individual=rownames(fullTab))
for(i in colnames(fullTab)){
  fTab[[i]] <- fullTab[,i]
}
fTab$Batch <- as.character(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][1]}))
fTab$Mix <- as.character(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][2]}))
fTab$Treatment <- as.character(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][3]}))
fTab$Replicate <- as.character(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][4]}))
fTab$Donor <- as.character(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][5]}))
fTab$LibSize <- as.numeric(sapply(fTab$Individual, function(x){strsplit(x,"_")[[1]][6]}))


##Example code for starting a more standard model (using edgeR)
abundances <- table(meta$MergedAnnotation, meta$Individual)
abundances <- unclass(abundances)
extra.info <- meta[match(colnames(abundances), meta$Individual),c("Genotype","Individual","Treatment","Mix","Batch")]
extra.info$Protocol <- ifelse(extra.info$Batch %in% c("A","B","C"), "MultiDonor","SingleDonor")
extra.info$Batch[extra.info$Batch %in% c("C","SingleDonorB")] <- "Nov23"

runEdgeR <- function(Protocol="MultiDonor",
                     Donor="CW"){
  TMP <- DGEList(counts=abundances[,extra.info$Genotype==Donor &
                                     extra.info$Protocol == Protocol &
                                     extra.info$Treatment %in% c("Control","VPA")],
                 samples=extra.info[extra.info$Genotype==Donor &
                                      extra.info$Protocol == Protocol &
                                      extra.info$Treatment %in% c("Control","VPA"),],
                 lib.size = colSums(abundances[,extra.info$Genotype==Donor &
                                                 extra.info$Protocol == Protocol &
                                                 extra.info$Treatment %in% c("Control","VPA")]))
  myDesign <- model.matrix(~factor(Batch) + factor(Treatment), TMP$samples)
  TMP <- calcNormFactors(TMP, method="TMM")
  TMP <- estimateDisp(TMP, myDesign)
  TMP <- glmQLFit(TMP, myDesign, robust=TRUE, abundance.trend=FALSE)
  resTMP <- glmQLFTest(TMP, coef=ncol(myDesign))
  tmp <- topTags(resTMP)
  return(as.data.frame(tmp))
}


md_cw <- runEdgeR()
sd_cw <- runEdgeR(Protocol = "SingleDonor")
md_pgp <- runEdgeR(Donor = "PGP1")
sd_pgp <- runEdgeR(Protocol = "SingleDonor",Donor = "PGP1")
md_mito <- runEdgeR(Donor = "Mito210")
sd_mito <- runEdgeR(Protocol = "SingleDonor",Donor = "Mito210")
md_h1 <- runEdgeR(Donor = "H1")
sd_h1 <- runEdgeR(Protocol = "SingleDonor",Donor = "H1")
md_11a <- runEdgeR(Donor = "11a")
sd_11a <- runEdgeR(Protocol = "SingleDonor",Donor = "11a")

md_cw <- md_cw[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
sd_cw <- sd_cw[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
md_pgp <- md_pgp[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
sd_pgp <- sd_pgp[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
md_mito <- md_mito[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
sd_mito <- sd_mito[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
md_h1 <- md_h1[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
sd_h1 <- sd_h1[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
md_11a <- md_11a[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]
sd_11a <- sd_11a[c("aRG","oRG","IP","PN","CFuPN","CPN","IN","CR","ChPL"),]

plot(sd_cw[,"logFC"], md_cw[,"logFC"])
cor(sd_cw[,"logFC"], md_cw[,"logFC"])
plot(sd_pgp[,"logFC"], md_pgp[,"logFC"])
cor(sd_pgp[,"logFC"], md_pgp[,"logFC"])
plot(sd_mito[,"logFC"], md_mito[,"logFC"])
cor(sd_mito[,"logFC"], md_mito[,"logFC"])
plot(sd_h1[,"logFC"], md_h1[,"logFC"])
cor(sd_h1[,"logFC"], md_h1[,"logFC"])
plot(sd_11a[,"logFC"], md_11a[,"logFC"])
cor(sd_11a[,"logFC"], md_11a[,"logFC"])

thing <- list(md_cw, md_pgp, md_mito, md_h1, md_11a,
              sd_cw, sd_pgp, sd_mito, sd_h1, sd_11a)
names(thing) <- c("md_cw", "md_pgp", "md_mito", "md_h1", "md_11a",
                  "sd_cw", "sd_pgp", "sd_mito", "sd_h1", "sd_11a")
corMat <- matrix(nrow=length(thing), ncol=length(thing), data=NA)
colnames(corMat) <- names(thing)
rownames(corMat) <- names(thing)
for(i in names(thing)){
  for(j in names(thing)){
    corMat[i, j] <- cor(thing[[i]]$logFC, thing[[j]]$logFC)
  }
}
corMat2 <- matrix(nrow=length(thing), ncol=length(thing), data=NA)
colnames(corMat2) <- names(thing)
rownames(corMat2) <- names(thing)
for(i in names(thing)){
  for(j in names(thing)){
    corMat2[i, j] <- cor(thing[[i]]$FDR, thing[[j]]$FDR)
  }
}
library(ComplexHeatmap)
Heatmap(corMat, column_order = 1:10, row_order = 1:10)

meltMat <- melt(corMat)
meltMat2 <- melt(corMat2)
meltMat$FDR_Cor <- meltMat2$value
ggplot(meltMat[grepl("md", meltMat$Var1) & grepl("sd", meltMat$Var2),], aes(x=Var1, y=Var2, size=FDR_Cor, fill=value)) +
  geom_point(pch=21) +
  scale_size_continuous(name="P-val Correlation", range = c(2,15)) +
  scale_fill_gradient2(low="darkblue", mid="white", high="red", midpoint=0.7, name="FoldChange Correlation") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        panel.background = element_blank())

md_cw["IN",]
sd_cw["IN",]
md_pgp["IN",]
sd_pgp["IN",]
md_mito["IN",]
sd_mito["IN",]
md_h1["IN",]
sd_h1["IN",]
md_11a["IN",]
sd_11a["IN",]


temp <- table(meta$Individual, meta$MergedAnnotation)
temp <- temp/apply(temp,1,sum)*100
myTab <- melt(temp)
colnames(myTab) <- c("Individual","CellType","Percent")
myTab$Individual <- as.character(myTab$Individual)
myTab$Batch <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][1]}))
myTab$Mix <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][2]}))
myTab$Treatment <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][3]}))
myTab$Replicate <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][4]}))
myTab$Donor <- as.character(sapply(myTab$Individual, function(x){strsplit(x,"_")[[1]][5]}))

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
myTab$Protocol <- ifelse(myTab$Batch %in% c("BatchA","BatchB","BatchC"), "MultiDonor","SingleDonor")
pt2 <- data_summary(myTab, varname="Percent", groupnames=c("Treatment","CellType","Protocol","Donor"))
pt2$Donor[pt2$Donor=="CIRM17"] <- "CW"
pt2$FDR <- NA
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR")){
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="CW" & pt2$Protocol=="MultiDonor"] <- paste0("p ~ ",round(md_cw[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="CW" & pt2$Protocol=="SingleDonor"] <- paste0("p ~ ",round(sd_cw[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="PGP1" & pt2$Protocol=="MultiDonor"] <- paste0("p ~ ",round(md_pgp[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="PGP1" & pt2$Protocol=="SingleDonor"] <- paste0("p ~ ",round(sd_pgp[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="Mito210" & pt2$Protocol=="MultiDonor"] <- paste0("p ~ ",round(md_mito[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="Mito210" & pt2$Protocol=="SingleDonor"] <- paste0("p ~ ",round(sd_mito[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="H1" & pt2$Protocol=="MultiDonor"] <- paste0("p ~ ",round(md_h1[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="H1" & pt2$Protocol=="SingleDonor"] <- paste0("p ~ ",round(sd_h1[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="11a" & pt2$Protocol=="MultiDonor"] <- paste0("p ~ ",round(md_11a[i,"FDR"],digits=4))
  pt2$FDR[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="11a" & pt2$Protocol=="SingleDonor"] <- paste0("p ~ ",round(sd_11a[i,"FDR"],digits=4))
}

pt2$LogFC <- NA
for(i in c("aRG","oRG","IP","PN","CFuPN","CPN","IN","ChPL","CR")){
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="CW" & pt2$Protocol=="MultiDonor"] <- paste0("logFC = ",round(md_cw[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="CW" & pt2$Protocol=="SingleDonor"] <- paste0("logFC = ",round(sd_cw[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="PGP1" & pt2$Protocol=="MultiDonor"] <- paste0("logFC = ",round(md_pgp[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="PGP1" & pt2$Protocol=="SingleDonor"] <- paste0("logFC = ",round(sd_pgp[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="Mito210" & pt2$Protocol=="MultiDonor"] <- paste0("logFC = ",round(md_mito[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="Mito210" & pt2$Protocol=="SingleDonor"] <- paste0("logFC = ",round(sd_mito[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="H1" & pt2$Protocol=="MultiDonor"] <- paste0("logFC = ",round(md_h1[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="H1" & pt2$Protocol=="SingleDonor"] <- paste0("logFC = ",round(sd_h1[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="11a" & pt2$Protocol=="MultiDonor"] <- paste0("logFC = ",round(md_11a[i,"logFC"],digits=2))
  pt2$LogFC[pt2$Treatment=="VPA" & pt2$CellType== i & pt2$Donor=="11a" & pt2$Protocol=="SingleDonor"] <- paste0("logFC = ",round(sd_11a[i,"logFC"],digits=2))
}

p1 <- ggplot(pt2[pt2$CellType=="IN" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=median, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=Qdown, ymax=Qup), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

p2 <- ggplot(pt2[pt2$CellType=="IN" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold")

cowplot::plot_grid(p1, p2)

ggplot(pt2[pt2$CellType=="IN" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold") +
  ggtitle("IN proportion\n(median, upper/lower quartile)")
ggplot(pt2[pt2$CellType=="IN" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=LogFC), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold") +
  ggtitle("IN proportion\n(median, upper/lower quartile)")

ggplot(pt2[pt2$CellType=="IP" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold") +
  ggtitle("IP proportion\n(median, upper/lower quartile)")

ggplot(pt2[pt2$CellType=="CPN" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold") +
  ggtitle("CPN proportion\n(median, upper/lower quartile)")

ggplot(pt2[pt2$CellType=="ChPL" & pt2$Treatment %in% c("Control","VPA"),], aes(x=Treatment, y=mean, fill=Donor)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=DonorCols[unique(pt2$Donor)]) +
  facet_grid(Protocol~Donor) +
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=Treatment, y=mean+sd+1, label=FDR), inherit.aes = FALSE,
            size = 4, color = "black", fontface = "bold") +
  ggtitle("ChPL proportion\n(median, upper/lower quartile)")


## What even is this?
myTab$TypeOfInterest <- myTab[[i]]
RandMod  <- glmer.nb(IN ~ Treatment + offset(log(LibSize)) + Batch, data=fTab[fTab$Donor =="PGP1" & fTab$Batch %in% c("BatchA","BatchB","BatchC"),])
InterMod <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Batch) + (1|Mix) + (1|Sample) + (Treatment|Genotype), data=myTab[myTab$Treatment!="EtOH",])
vpaP <- c()
for(i in names(RandMod)){
  print(i)
  print(anova(RandMod[[i]], InterMod[[i]])[2,8])
  vpaP[i] <- anova(RandMod[[i]], InterMod[[i]])[2,8]
}
vpaQ <- p.adjust(vpaP, method="BH")
vpaQ[order(vpaQ)]


library(Seurat)
seur <- readRDS("/n/arlotta_lab2/Users/Tyler/PublicData/Human/Krieg2023/local.rds")
seur <- subset(seur, subset=region %in% c("BA10", "BA22", "BA24", "BA46", "BA8", "BA9", "FIC", "INS", "PFC") & 
                 age_range %in% c("0-1 years","1-2 years","2-4 years","4-10 years","10-20 years"))
saveRDS(seur, "/n/arlotta_lab2/Users/Tyler/PublicData/Human/KriegForKaren/ChildCortex.RDS")







multi <- meta[meta$Protocol=="MultiDonorNSC",]



## Example code for starting a mixed model (using lme4)
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

##Example code for starting a more standard model (using edgeR)
abundances <- table(seur$Genotype[seur$Treatment != "ETOH" & seur$Mix=="M1"], seur$Sample[seur$Treatment!="ETOH"& seur$Mix=="M1"]) 
abundances <- unclass(abundances)
extra.info <- seur[match(colnames(abundances), seur$Sample),c("Genotype","Sample","Treatment","Mix","Batch")]
vpaWholeChim3 <- DGEList(counts=abundances, samples=extra.info, lib.size = colSums(abundances))
design4 <- model.matrix(~factor(Batch) + factor(Treatment), vpaWholeChim3$samples)
vpaWholeChim3 <- calcNormFactors(vpaWholeChim3, method="TMM")
vpaWholeChim3 <- estimateDisp(vpaWholeChim3, design4, trend="none")
vpaWholeChim3 <- glmQLFit(vpaWholeChim3, design4, robust=TRUE, abundance.trend=FALSE)
resVPAWholeChim1 <- glmQLFTest(vpaWholeChim3, coef=ncol(design4))
rvwc1 <- topTags(resVPAWholeChim1)
topTags(resVPAWholeChim1)



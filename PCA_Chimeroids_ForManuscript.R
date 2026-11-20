library(multcomp)
library(emmeans)
library(multcompView)
library(ggplot2)
library(plyr)
library(reshape2)

`%ni%` = Negate(`%in%`)
## Figure 1:
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

## PCA analyses on pseudobulked data
chimP <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PseudoBulk_CellType_UsingAnnotationFromNovember17b.RDS")
chimM <- readRDS("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PseudoBulk_CellType_Metadata_UsingAnnotationFromNovember17b.RDS")
mB <- chimM

head(chimM)
chimM$Dataset <- chimM$Protocol
myMeta <- chimM[,c("Sample","Treatment","CellType","Genotype","Batch","Mix","Dataset")]
myMeta$Genotype[myMeta$Genotype %in% c("CIRM","CIRM17")] <- "CW"
myPseudo <- chimP[,myMeta$Dataset=="MultiDonor"]
myMeta <- myMeta[myMeta$Dataset=="MultiDonor",]
myPseudo <- myPseudo[,!grepl("NPC", myMeta$Sample)]
myMeta <- myMeta[!grepl("NPC", myMeta$Sample),]
myCPM <- t(t(myPseudo)/apply(myPseudo,2,sum)*1000000)
myLogCPM <- log(myCPM+1)


## Make a function to run this PCA stuff on subgroups of samples; eg, run PCA on just the VPA+control IPs
makeMyPCA_Chim <- function(celltype="IP", treatment=c("Control","VPA")){
  tmpM <- myMeta[myMeta$CellType==celltype & myMeta$Treatment %in% treatment, c("Sample","Treatment","CellType","Genotype","Batch","Mix","Dataset")]
  tmpP <- myLogCPM[,myMeta$CellType==celltype & myMeta$Treatment %in% treatment]
  tmpP <- tmpP[apply(tmpP,1,mean)>1,]
  varGenes <- apply(tmpP,1,function(x){var(x)/mean(x)})
  topVarGenes <- rownames(tmpP)[(log(varGenes) + apply(tmpP,1,mean))>1]
  tmpPCA <- prcomp(t(tmpP[topVarGenes,]))
  ## Run a filter for outliers:
  # tmpM$Distance <- apply(tmpPCA$x, 1, function(x){sqrt(sum(x*x))})
  # cutoffB <- quantile(tmpM$Distance)[[4]] + 5*(quantile(tmpM$Distance)[[4]]-quantile(tmpM$Distance)[[2]])
  # print(paste("Failing:", paste(rownames(tmpM)[tmpM$Distance >= cutoffB], collapse=", ")))
  # tmpP <- tmpP[,tmpM$Distance < cutoffB]
  # tmpM <- tmpM[tmpM$Distance < cutoffB,]
  # varGenes <- apply(tmpP,1,function(x){var(x)/mean(x)})
  # topVarGenes <- rownames(tmpP)[(log(varGenes) + apply(tmpP,1,mean))>1]
  # print(paste("Number of genes used =", length(topVarGenes)))
  # tmpPCA <- prcomp(t(tmpP[topVarGenes,]))
  tmpM$Distance <- apply(tmpPCA$x, 1, function(x){sqrt(sum(x*x))})
  tmpM <- cbind(tmpM, tmpPCA$x)
  return(list(tmpM,tmpPCA))
}
VPA_IP2 <- makeMyPCA_Chim()
VPA_oRG2 <- makeMyPCA_Chim(celltype = "oRG")
VPA_aRG2 <- makeMyPCA_Chim(celltype = "aRG")
VPA_CPN2 <- makeMyPCA_Chim(celltype = "CPN")
VPA_CFuPN2 <- makeMyPCA_Chim(celltype = "CFuPN")
VPA_PN2 <- makeMyPCA_Chim(celltype = "PN")
VPA_ChPL2 <- makeMyPCA_Chim(celltype = "ChPL")
VPA_IN2 <- makeMyPCA_Chim(celltype = "IN")
VPA_CR2 <- makeMyPCA_Chim(celltype="CR")

makePCAPlot <- function(inList, PCa="PC1", PCb="PC2"){
  inList[[1]][["PCa"]] <- inList[[1]][[PCa]]
  inList[[1]][["PCb"]] <- inList[[1]][[PCb]]
  myTitle <- paste(unique(inList[[1]]$Treatment)[2],inList[[1]]$CellType[1],sep="-treated ")
  percVar <- (inList[[2]]$sdev^2)/sum(inList[[2]]$sdev^2)*100
  names(percVar) <- colnames(inList[[2]]$x)
  ggplot(inList[[1]], aes(x=PCa, y=PCb, color=Genotype, shape=Treatment)) +
    geom_point(size=3) +
    scale_color_manual(values=myCols2[unique(inList[[1]][["Genotype"]])]) +
    ggtitle(myTitle) +
    theme(axis.title = element_text(face="bold"),
          strip.text = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust=0.5)) +
    xlab(paste0(PCa," (", round(percVar[[PCa]], digits=1), "% of variance)"))+
    ylab(paste0(PCb," (", round(percVar[[PCb]], digits=1), "% of variance)"))
}
makePCAPlot(inList=VPA_IP2, PCa="PC2", PCb="PC3")
makePCAPlot(inList=VPA_IP2, PCa="PC1", PCb="PC2")

for(i in c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")){
  print(paste0(i, ": Batch: ",round(cor(VPA_IP2[[1]][[i]],as.numeric(factor(VPA_IP2[[1]]$Batch))), digits=3),
               "    Donor: ",round(cor(VPA_IP2[[1]][[i]],as.numeric(factor(VPA_IP2[[1]]$Genotype))), digits=3),
               "    Treatment: ",round(cor(VPA_IP2[[1]][[i]],as.numeric(factor(VPA_IP2[[1]]$Treatment))), digits=3)))
}


##Ok, go through all of that, but with the whole-organoid data instead of the split-by-celltype:
cp2 <- readRDS("PseudoBulk_Global_UsingAnnotationFromNovember17b.RDS")
cm2 <- readRDS("PseudoBulk_Global_Metadata_UsingAnnotationFromNovember17b.RDS")
head(cm2)
cm2$Dataset <- cm2$Protocol
myMeta2 <- cm2[,c("Sample","Treatment","Genotype","Batch","Mix","Dataset")]
myMeta2$Genotype[myMeta2$Genotype %in% c("CIRM","CIRM17")] <- "CW"
myPseudo2 <- cp2[,myMeta2$Dataset=="MultiDonor"]
myMeta2 <- myMeta2[myMeta2$Dataset=="MultiDonor",]
myCPM2 <- t(t(myPseudo2)/apply(myPseudo2,2,sum)*1000000)
myLogCPM2 <- log(myCPM2+1)

makeMyPCA2_Chim <- function(treatment=c("Control","VPA")){
  tmpM <- myMeta2[myMeta2$Treatment %in% treatment, c("Sample","Treatment","Genotype","Batch","Mix","Dataset")]
  tmpP <- myLogCPM2[,myMeta2$Treatment %in% treatment]
  tmpP <- tmpP[apply(tmpP,1,mean)>1,]
  varGenes <- apply(tmpP,1,function(x){var(x)/mean(x)})
  topVarGenes <- rownames(tmpP)[(log(varGenes) + apply(tmpP,1,mean))>1]
  tmpPCA <- prcomp(t(tmpP[topVarGenes,]))
  ## Run a filter for outliers:
  tmpM$Distance <- apply(tmpPCA$x, 1, function(x){sqrt(sum(x*x))})
  cutoffB <- quantile(tmpM$Distance)[[4]] + 2.5*(quantile(tmpM$Distance)[[4]]-quantile(tmpM$Distance)[[2]])
  print(paste("Passing:", paste(rownames(tmpM)[tmpM$Distance < cutoffB], collapse=", ")))
  print(paste("Failing:", paste(rownames(tmpM)[tmpM$Distance >= cutoffB], collapse=", ")))
  print(paste0("Cutoff = ", cutoffB, "; Highest passing = ", max(tmpM$Distance[tmpM$Distance < cutoffB]),
               "; Lowest failing = ", min(tmpM$Distance[tmpM$Distance >= cutoffB])))
  tmpP <- tmpP[,tmpM$Distance < cutoffB]
  tmpM <- tmpM[tmpM$Distance < cutoffB,]
  varGenes <- apply(tmpP,1,function(x){var(x)/mean(x)})
  topVarGenes <- rownames(tmpP)[(log(varGenes) + apply(tmpP,1,mean))>1]
  print(paste("Number of genes used =", length(varGenes)))
  tmpPCA <- prcomp(t(tmpP[topVarGenes,]))
  tmpM$Distance <- apply(tmpPCA$x, 1, function(x){sqrt(sum(x*x))})
  tmpM <- cbind(tmpM, tmpPCA$x)
  return(list(tmpM,tmpPCA))
}
VPA_Global <- makeMyPCA2_Chim()




## Streamline a little bit, and just do what's necessary to get to the ANOVA stuff for each cell type:
## 
## 
##
library(ggbeeswarm)
#celltype <- "CPN"
for(celltype in c("aRG","oRG","IP","CFuPN","CPN","PN","ChPL","CR")){
  print(celltype)
  ETOH_Temp <- makeMyPCA_Chim(celltype=celltype, treatment=c("Control","VPA"))
  tmpdist <- as.matrix(dist(ETOH_Temp[[2]]$x))
  etoDistInfo <- data.frame(SampleA=c("No"), TreatmentA=c("No"), CTA=c("No"), GenotypeA=c("No"), BatchA=c("No"),MixA=c("No"),
                            SampleB=c("No"), TreatmentB=c("No"), CTB=c("No"), GenotypeB=c("No"), BatchB=c("No"),MixB=c("No"),
                            Distance=c(100))
  for(i in 1:(nrow(ETOH_Temp[[1]])-1)){
    for(j in (i+1):nrow(ETOH_Temp[[1]])){
      tmp <- data.frame(c(ETOH_Temp[[1]][i,1:6], ETOH_Temp[[1]][j,1:6], tmpdist[i,j]))
      colnames(tmp) <- colnames(etoDistInfo)
      etoDistInfo <- rbind(etoDistInfo, tmp)
    }
  }
  etoDistInfo <- etoDistInfo[-1,]
  etoDistInfo$TreatmentCompare <- paste(etoDistInfo$TreatmentA, etoDistInfo$TreatmentB, sep="-")
  etoDistInfo$GenoCompare <- paste(etoDistInfo$GenotypeA, etoDistInfo$GenotypeB, sep="-")
  etoDistInfo$MonoTreatment <- etoDistInfo$TreatmentA==etoDistInfo$TreatmentB
  etoDistInfo$MonoGenotype <- etoDistInfo$GenotypeA==etoDistInfo$GenotypeB
  etoDistInfo$WhichGen <- etoDistInfo$MonoGenotype
  etoDistInfo$WhichGen[etoDistInfo$WhichGen] <- etoDistInfo$GenotypeA[etoDistInfo$WhichGen]
  etoDistInfo$WhichGen[!etoDistInfo$MonoGenotype] <- "Between\nDonors"
  etoDistInfo$WhichGen <- factor(etoDistInfo$WhichGen, levels=c("CW","Mito210","H1","PGP1","11a","Between\nDonors"))
  etoDistInfo$BatchBatch <- paste(etoDistInfo$BatchA, etoDistInfo$BatchB, sep="-")
  etoDistInfo$MixMix <- paste(etoDistInfo$MixA, etoDistInfo$MixB, sep="-")
  ## Let's run some ANOVAs, I guess:
  Chim_Response <- etoDistInfo[etoDistInfo$TreatmentCompare %in% c("VPA-Control","Control-VPA") & etoDistInfo$MonoGenotype,]
  chimrespAnova <- aov(Distance ~ BatchBatch + WhichGen, dat=Chim_Response)
  summary(chimrespAnova)
  #Check normality:
  hist(chimrespAnova$residuals)
  car::qqPlot(chimrespAnova$residuals, id=FALSE)
  shapiro.test(chimrespAnova$residuals)
  #oneway.test(Distance~WhichGen, data=Chim_Response, var.equal=FALSE)
  ## Good;
  post_test <- glht(chimrespAnova, linfct = mcp(WhichGen = "Tukey"))
  summary(post_test)
  emmeans(chimrespAnova, list(pairwise ~ WhichGen), adjust = "tukey")
  TUKEY <- TukeyHSD(chimrespAnova, conf.level=.95)
  generate_label_df <- function(TUKEY, variable){
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    Tukey.labels$treatment=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
    return(Tukey.labels)
  }
  LABELS <- generate_label_df(TUKEY , "WhichGen")
  myBox <- boxplot(Chim_Response$Distance ~ Chim_Response$WhichGen)
  over <- 0.1*max( myBox$stats[nrow(myBox$stats),] )
  text( c(1:nlevels(Chim_Response$WhichGen)) , myBox$stats[nrow(myBox$stats),]+over , LABELS[,1] )
  LABELS$Positions <- sapply(LABELS$treatment, function(x){max(Chim_Response$Distance[Chim_Response$WhichGen==x]+4)})
  ggplot(Chim_Response, aes(x=WhichGen, y=Distance, fill=WhichGen)) +
    geom_violin(width=0.8) +
    scale_fill_manual(values=myCols2[c("CW","Mito210","H1","PGP1","11a")]) +
    theme(axis.text.x=element_text(angle=45, face="bold", size=12, hjust=1),
          axis.title=element_text(face="bold", size=16),
          panel.background = element_rect(fill="white", color="white"),
          plot.title=element_text(face="bold", size=14, hjust=0.5)) +
    ggtitle(paste0("Within-Donor Control-to-VPA Distances\n", celltype, " Only; ANOVA p-value = ", signif(summary(chimrespAnova)[[1]][["Pr(>F)"]][2],2), "\n(Chimeroids)")) +
    xlab("") + ylab("Control-to-VPA Distance") +
    NoLegend() +
    geom_text(data=LABELS,aes(x=treatment, y=Positions, label=Letters), inherit.aes = FALSE,
              size = 4, color = "black", fontface = "bold") +
    geom_boxplot(width=0.2) +
    geom_jitter(size=2, pch=21, width=0.06)
    #geom_point(position=position_jitterdodge(dodge.width=0.8), size=2, pch=21)
  ggsave(paste0("PCA_Plots/Distance_PCA_VlnAndPoint_ChimeroidsOnly_VPA_",celltype,"_Tukey_",Sys.Date(),".svg"), height=4, width=5)
  print(paste("Done with",celltype))
}

##Global
for(i in 1){
  ETOH_Temp <- makeMyPCA2_Chim(treatment=c("Control","VPA"))
  tmpdist <- as.matrix(dist(ETOH_Temp[[2]]$x))
  etoDistInfo <- data.frame(SampleA=c("No"), TreatmentA=c("No"), GenotypeA=c("No"), BatchA=c("No"),MixA=c("No"),
                            SampleB=c("No"), TreatmentB=c("No"), GenotypeB=c("No"), BatchB=c("No"),MixB=c("No"),
                            Distance=c(100))
  for(i in 1:(nrow(ETOH_Temp[[1]])-1)){
    for(j in (i+1):nrow(ETOH_Temp[[1]])){
      tmp <- data.frame(c(ETOH_Temp[[1]][i,1:5], ETOH_Temp[[1]][j,1:5], tmpdist[i,j]))
      colnames(tmp) <- colnames(etoDistInfo)
      etoDistInfo <- rbind(etoDistInfo, tmp)
    }
  }
  etoDistInfo <- etoDistInfo[-1,]
  etoDistInfo$TreatmentCompare <- paste(etoDistInfo$TreatmentA, etoDistInfo$TreatmentB, sep="-")
  etoDistInfo$GenoCompare <- paste(etoDistInfo$GenotypeA, etoDistInfo$GenotypeB, sep="-")
  etoDistInfo$MonoTreatment <- etoDistInfo$TreatmentA==etoDistInfo$TreatmentB
  etoDistInfo$MonoGenotype <- etoDistInfo$GenotypeA==etoDistInfo$GenotypeB
  etoDistInfo$WhichGen <- etoDistInfo$MonoGenotype
  etoDistInfo$WhichGen[etoDistInfo$WhichGen] <- etoDistInfo$GenotypeA[etoDistInfo$WhichGen]
  etoDistInfo$WhichGen[!etoDistInfo$MonoGenotype] <- "Between\nDonors"
  etoDistInfo$WhichGen <- factor(etoDistInfo$WhichGen, levels=c("CW","Mito210","H1","PGP1","11a","Between\nDonors"))
  etoDistInfo$BatchBatch <- paste(etoDistInfo$BatchA, etoDistInfo$BatchB, sep="-")
  etoDistInfo$MixMix <- paste(etoDistInfo$MixA, etoDistInfo$MixB, sep="-")
  ## Let's run some ANOVAs, I guess:
  Chim_Response <- etoDistInfo[etoDistInfo$TreatmentCompare %in% c("VPA-Control","Control-VPA") & etoDistInfo$MonoGenotype,]
  chimrespAnova <- aov(Distance ~ BatchBatch + WhichGen, dat=Chim_Response)
  summary(chimrespAnova)
  #Check normality:
  hist(chimrespAnova$residuals)
  car::qqPlot(chimrespAnova$residuals, id=FALSE)
  shapiro.test(chimrespAnova$residuals)
  oneway.test(Distance~BatchBatch + WhichGen, data=Chim_Response, var.equal=FALSE)
  ## Good;
  post_test <- glht(chimrespAnova, linfct = mcp(WhichGen = "Tukey"))
  summary(post_test)
  emmeans(chimrespAnova, list(pairwise ~ WhichGen), adjust = "tukey")
  TUKEY <- TukeyHSD(chimrespAnova, conf.level=.95)
  generate_label_df <- function(TUKEY, variable){
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    Tukey.labels$treatment=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
    return(Tukey.labels)
  }
  LABELS <- generate_label_df(TUKEY , "WhichGen")
  myBox <- boxplot(Chim_Response$Distance ~ Chim_Response$WhichGen)
  over <- 0.1*max( myBox$stats[nrow(myBox$stats),] )
  text( c(1:nlevels(Chim_Response$WhichGen)) , myBox$stats[nrow(myBox$stats),]+over , LABELS[,1] )
  LABELS$Positions <- sapply(LABELS$treatment, function(x){max(Chim_Response$Distance[Chim_Response$WhichGen==x]+4)})
  ggplot(Chim_Response, aes(x=WhichGen, y=Distance, fill=WhichGen)) +
    geom_violin(width=0.8) +
    scale_fill_manual(values=myCols2[c("CW","Mito210","H1","PGP1","11a")]) +
    theme(axis.text.x=element_text(angle=45, face="bold", size=12, hjust=1),
          axis.title=element_text(face="bold", size=16),
          panel.background = element_rect(fill="white", color="white"),
          plot.title=element_text(face="bold", size=14, hjust=0.5)) +
    ggtitle(paste0("Within-Donor Control-to-VPA Distances\nGlobal; ANOVA p-value = ", signif(summary(chimrespAnova)[[1]][["Pr(>F)"]][2],2), "\n(Chimeroids)")) +
    xlab("") + ylab("Control-to-VPA Distance") +
    NoLegend() +
    geom_text(data=LABELS,aes(x=treatment, y=Positions, label=Letters), inherit.aes = FALSE,
              size = 4, color = "black", fontface = "bold") +
    geom_boxplot(width=0.2) +
    geom_jitter(size=2, pch=21, width=0.06)
  #geom_point(position=position_jitterdodge(dodge.width=0.8), size=2, pch=21)
  ggsave(paste0("PCA_Plots/Distance_PCA_VlnAndPoint_ChimeroidsOnly_VPA_Global_Tukey_",Sys.Date(),".svg"), height=4, width=5)
}



## Make some figures for the PCAs themselves:
makePCAFigure <- function(inList, PCa="PC1", PCb="PC2"){
  inMat <- inList[[1]]
  inMat$PCa <- inMat[[PCa]]
  inMat$PCb <- inMat[[PCb]]
  ##Make centroids for each donor/treatment
  myCentroids <- data.frame(Sample="None", Treatment="None", Genotype="None", PCa=c(0), PCb=c(0), Group="None")
  for(i in unique(inMat$Treatment)){
    for(j in unique(inMat$Genotype)){
        tmp <- data.frame(Sample=paste(i,j, sep="-"), Treatment=i, Genotype=j,
                          PCa=c(mean(as.numeric(inMat[inMat$Treatment==i & inMat$Genotype==j,"PCa"]))),
                          PCb=c(mean(as.numeric(inMat[inMat$Treatment==i & inMat$Genotype==j,"PCb"]))),
                          Group=paste(j))
        myCentroids <- rbind(myCentroids, tmp)
    }
  }
  myCentroids <- myCentroids[-1,]
  myTitle <- paste(unique(inMat$Treatment)[2],inMat$CellType[1],sep="-treated ")
  percVar <- (inList[[2]]$sdev^2)/sum(inList[[2]]$sdev^2)*100
  names(percVar) <- colnames(inList[[2]]$x)
  ggplot(inMat, aes(x=PCa, y=PCb, color=Genotype, fill=Treatment, shape=Treatment)) +
    geom_point(size=3, stroke=1.75, alpha=0.4) +
    scale_shape_manual(values=c(1,17)) +
    scale_color_manual(values=myCols2[unique(inList[[1]][["Genotype"]])]) +
    scale_fill_manual(values=c("white", "black")) +
    ggtitle(myTitle) +
    theme(axis.title = element_text(face="bold"),
          strip.text = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust=0.5),
          panel.background = element_rect(fill="white", color="grey85"),
          panel.grid.major = element_line(color = "grey90", size = 0.3, linetype = 1)) +
    xlab(paste0(PCa," (", round(percVar[[PCa]], digits=1), "% of variance)")) +
    ylab(paste0(PCb," (", round(percVar[[PCb]], digits=1), "% of variance)")) +
    geom_point(data=myCentroids[myCentroids$Treatment=="Control",], inherit.aes = FALSE,
               aes(x=PCa, y=PCb, shape=Treatment), size=1, shape=16, color="black", stroke=1) +
    geom_point(data=myCentroids[myCentroids$Treatment!="Control",], inherit.aes = FALSE,
               aes(x=PCa, y=PCb, shape=Treatment), size=1, shape=17, color="black", stroke=1) +
    geom_line(data=myCentroids, aes(group=Group), lty=1, color="black")
}
makePCAFigure(VPA_CFuPN2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_CFuPN_",Sys.Date(),".svg"), width=5, height=4)

makePCAFigure(VPA_CPN2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_CPN_",Sys.Date(),".svg"), width=5, height=4)

makePCAFigure(VPA_oRG2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_oRG_",Sys.Date(),".svg"), width=5, height=4)

makePCAFigure(VPA_aRG2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_aRG_",Sys.Date(),".svg"), width=5, height=4)


makePCAFigure(VPA_IP2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_IP_",Sys.Date(),".svg"), width=5, height=4)


makePCAFigure(VPA_PN2)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_PN_",Sys.Date(),".svg"), width=5, height=4)


makeGlobalPCAFigure <- function(inList, PCa="PC1", PCb="PC2"){
  inMat <- inList[[1]]
  inMat$PCa <- inMat[[PCa]]
  inMat$PCb <- inMat[[PCb]]
  ##Make centroids for each donor/treatment
  myCentroids <- data.frame(Sample="None", Treatment="None", Genotype="None", PCa=c(0), PCb=c(0), Group="None")
  for(i in unique(inMat$Treatment)){
    for(j in unique(inMat$Genotype)){
      tmp <- data.frame(Sample=paste(i,j, sep="-"), Treatment=i, Genotype=j,
                        PCa=c(mean(as.numeric(inMat[inMat$Treatment==i & inMat$Genotype==j,"PCa"]))),
                        PCb=c(mean(as.numeric(inMat[inMat$Treatment==i & inMat$Genotype==j,"PCb"]))),
                        Group=paste(j))
      myCentroids <- rbind(myCentroids, tmp)
    }
  }
  myCentroids <- myCentroids[-1,]
  myTitle <- paste0("Global ",unique(inMat$Treatment)[2],"-treated")
  percVar <- (inList[[2]]$sdev^2)/sum(inList[[2]]$sdev^2)*100
  names(percVar) <- colnames(inList[[2]]$x)
  ggplot(inMat, aes(x=PCa, y=PCb, color=Genotype, fill=Treatment, shape=Treatment)) +
    geom_point(size=3, stroke=1.75, alpha=0.4) +
    scale_shape_manual(values=c(1,17)) +
    scale_color_manual(values=myCols2[unique(inList[[1]][["Genotype"]])]) +
    scale_fill_manual(values=c("white", "black")) +
    ggtitle(myTitle) +
    theme(axis.title = element_text(face="bold"),
          strip.text = element_text(face="bold"),
          plot.title = element_text(face="bold", hjust=0.5),
          panel.background = element_rect(fill="white", color="grey85"),
          panel.grid.major = element_line(color = "grey90", size = 0.3, linetype = 1)) +
    xlab(paste0(PCa," (", round(percVar[[PCa]], digits=1), "% of variance)")) +
    ylab(paste0(PCb," (", round(percVar[[PCb]], digits=1), "% of variance)")) +
    geom_point(data=myCentroids[myCentroids$Treatment=="Control",], inherit.aes = FALSE,
               aes(x=PCa, y=PCb, shape=Treatment), size=1, shape=16, color="black", stroke=1) +
    geom_point(data=myCentroids[myCentroids$Treatment!="Control",], inherit.aes = FALSE,
               aes(x=PCa, y=PCb, shape=Treatment), size=1, shape=17, color="black", stroke=1) +
    geom_line(data=myCentroids, aes(group=Group), lty=1, color="black")
}
makeGlobalPCAFigure(VPA_Global)
ggsave(paste0("/Users/tfaits/Documents/Chimera/OctoberSequencing2023/PCA_Plots/PCA_Plot_VPA_Chimeroids_Global_",Sys.Date(),".svg"), width=5, height=4)

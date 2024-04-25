## Merge_Seurats
## This assumes all objects are processed; it removes SCTransform assays (if they exist) to merge on raw UMI counts
## Inputs are a comma-separated list of input Seurat objects, followed by the name/location of the intended output file
## Rscript Merge_Seurats.R <inSeur1,inSeur2,...,inSeurN> <MergedSeurat.rds">
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
inSeurs <- as.vector(sapply(args[1], function(x){strsplit(x,",")[[1]]}))
outFile <- args[2]

myList <- c()
for(i in inSeurs){
  myList[[i]] <- readRDS(i)
  DefaultAssay(myList[[i]]) <- "RNA"
  myList[[i]][["SCT"]] <- NULL
}

seur <- merge(myList[[1]], myList[2:length(myList)])
rm(myList)

seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
print("Merge complete. Beginning SCTransform.")
seur <- SCTransform(seur, vars.to.regress="percent.mt")
seur <- RunPCA(seur, npcs=30)
seur <- RunUMAP(seur, reduction = "pca", dims = 1:30)
seur <- FindNeighbors(seur, dims=1:30)
seur <- FindClusters(seur, resolution=0.8)
saveRDS(seur, outFile)

## Make Seurat object from 10X outputs
## Make Seurat from an input folder name
library(Seurat)
args <- commandArgs(trailingOnly = T)
inFolder <- args[1]
sampleName <- strsplit(inFolder, "/")[[1]]
sampleName <- sampleName[length(sampleName)]
startData <- Read10X(data.dir = paste0(inFolder,"/outs/filtered_feature_bc_matrix"))
# Initialize the Seurat object with the raw (non-normalized data).
seur <- CreateSeuratObject(counts = startData, project = sampleName, min.cells = 3, min.features = 200)
seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
seur[["percent.rb"]] <- PercentageFeatureSet(seur, pattern = "^RP[SL]")
seur <- subset(seur, subset=nCount_RNA>500 & nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)
seur <- SCTransform(seur, vars.to.regress = "percent.mt")
seur <- RunPCA(seur, npcs = 30)
seur <- RunUMAP(seur, dims=1:30)
seur <- FindNeighbors(seur, dims=1:30)
seur <- FindClusters(seur, resolution = 0.3)
print(ncol(seur))
saveRDS(seur, paste0(sampleName, "_Seurat_BasicQC.RDS"))

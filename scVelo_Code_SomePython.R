## Convert the Seurat object into an Anndata object (in preparation for adding in velocyto .loom data)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

seur <- readRDS("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/Seurats/AllMerged.RDS")
seur@meta.data <- readRDS("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/Seurats/Metadata_AllMerged_Streamlined_November17_2023.RDS")

seur$barcode <- colnames(seur)
seur$Suffix <- sapply(colnames(seur), function(x){strsplit(x,"-")[[1]][2]})## This assumes that there's only one dash in the cell names
# seur$UMAP_1 <- seur@reductions$umap@cell.embeddings[,1]
# seur$UMAP_2 <- seur@reductions$umap@cell.embeddings[,2]
write.csv(seur@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seur, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seur@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)


## Prep in python
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['MergedAnnotation'], frameon=False, save=True)


# save dataset as anndata format
adata.write('my_data.h5ad')



#####
#####
## PYTHON STUFF ##
#####
#####
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('my_data.h5ad')
# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M1_CNT_1/velocyto/CH16_M1_CNT_1.loom", cache=True)
ldata2 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M1_CNT_2/velocyto/CH16_M1_CNT_2.loom", cache=True)
ldata3 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M1_CNT_3/velocyto/CH16_M1_CNT_3.loom", cache=True)
ldata4 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M1_VPA_1/velocyto/CH16_M1_VPA_1.loom", cache=True)
ldata5 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M1_VPA_2/velocyto/CH16_M1_VPA_2.loom", cache=True)
ldata6 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M2_CNT_1/velocyto/CH16_M2_CNT_1.loom", cache=True)
ldata7 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M2_CNT_2/velocyto/CH16_M2_CNT_2.loom", cache=True)
ldata8 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M2_VPA_1/velocyto/CH16_M2_VPA_1.loom", cache=True)
ldata9 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_M2_VPA_2/velocyto/CH16_M2_VPA_2.loom", cache=True)
ldata10 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_SingleDonor_R1_1/velocyto/CH16_SingleDonor_R1_1.loom", cache=True)
ldata11 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_SingleDonor_R1_2/velocyto/CH16_SingleDonor_R1_2.loom", cache=True)
ldata12 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_SingleDonor_R2_1/velocyto/CH16_SingleDonor_R2_1.loom", cache=True)
ldata13 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/October2023_ChimeroidRevisions/CellRangerOuts/CH16_SingleDonor_R2_2/velocyto/CH16_SingleDonor_R2_2.loom", cache=True)

ldata14 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S1_M1_CNT__2_SystemFailure/velocyto/R0_S1_M1_CNT__2_SystemFailure.loom", cache=True)
ldata15 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S2_M1_ETOH__2_SystemFailure/velocyto/R0_S2_M1_ETOH__2_SystemFailure.loom", cache=True)
ldata16 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S3_M1_VPA__2_SystemFailure/velocyto/R0_S3_M1_VPA__2_SystemFailure.loom", cache=True)
ldata17 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S4_MULTI_C1R2_SystemFailure/velocyto/R0_S4_MULTI_C1R2_SystemFailure.loom", cache=True)
ldata18 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S5_MULTI_C2R2_SystemFailure/velocyto/R0_S5_MULTI_C2R2_SystemFailure.loom", cache=True)
ldata19 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R0_S6_MULTI_C3R2_SystemFailure/velocyto/R0_S6_MULTI_C3R2_SystemFailure.loom", cache=True)
ldata20 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S1_MULTI_C1R3/velocyto/R2_S1_MULTI_C1R3.loom", cache=True)
ldata21 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S2_MULTI_C2R3_GoodRepeat/velocyto/R2_S2_MULTI_C2R3_GoodRepeat.loom", cache=True)
ldata22 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S3_MULTI_C3R3/velocyto/R2_S3_MULTI_C3R3.loom", cache=True)
ldata23 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S4_M2_CNT__2/velocyto/R2_S4_M2_CNT__2.loom", cache=True)
ldata24 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S5_M2_ETOH__2/velocyto/R2_S5_M2_ETOH__2.loom", cache=True)
ldata25 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S5_MULTI_C2R4_GoodRepeat/velocyto/R2_S5_MULTI_C2R4_GoodRepeat.loom", cache=True)
ldata26 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R2_S6_M2_VPA__2/velocyto/R2_S6_M2_VPA__2.loom", cache=True)
ldata27 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S1_M2_CNT__3/velocyto/R3_S1_M2_CNT__3.loom", cache=True)
ldata28 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S2_M2_ETOH__3/velocyto/R3_S2_M2_ETOH__3.loom", cache=True)
ldata29 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S3_M2_VPA__3/velocyto/R3_S3_M2_VPA__3.loom", cache=True)
ldata30 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S4_MULTI_C1R4/velocyto/R3_S4_MULTI_C1R4.loom", cache=True)
ldata31 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S6_M2_VPA__1_GoodRepeat/velocyto/R3_S6_M2_VPA__1_GoodRepeat.loom", cache=True)
ldata32 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R3_S6_MULTI_C3R4/velocyto/R3_S6_MULTI_C3R4.loom", cache=True)
ldata33 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R4_S1_M1_CNT__1/velocyto/R4_S1_M1_CNT__1.loom", cache=True)
ldata34 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R4_S2_M1_ETOH__1/velocyto/R4_S2_M1_ETOH__1.loom", cache=True)
ldata35 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R4_S3_M1_VPA__1/velocyto/R4_S3_M1_VPA__1.loom", cache=True)
ldata36 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R4_S4_M2_CNT__1/velocyto/R4_S4_M2_CNT__1.loom", cache=True)
ldata37 = scv.read("/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/Chimera/Data/CellrangerOuts/February_Runs/R4_S5_M2_ETOH__1/velocyto/R4_S5_M2_ETOH__1.loom", cache=True)

ldata38 = scv.read('/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/EtOH_VPA_July2022/CellRanger/VA_1/velocyto/VA_1.loom', cache=True)
ldata39 = scv.read('/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/EtOH_VPA_July2022/CellRanger/VA_2/velocyto/VA_2.loom', cache=True)
ldata40 = scv.read('/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/EtOH_VPA_July2022/CellRanger/ETOHre/velocyto/ETOHre.loom', cache=True)
ldata41 = scv.read('/net/bos-nfsisilon/ifs/rc_labs/arlotta_lab2/Users/Tyler/projects/EtOH_VPA_July2022/CellRanger/ETOH_2/velocyto/ETOH_2.loom', cache=True)


# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_1_2' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_2_2' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_3_2' for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_11_2' for bc in barcodes]
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_12_2' for bc in barcodes]
ldata5.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata6.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_4_2' for bc in barcodes]
ldata6.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata7.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_5_2' for bc in barcodes]
ldata7.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata8.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_13_2' for bc in barcodes]
ldata8.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata9.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_14_2' for bc in barcodes]
ldata9.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata10.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_2_1' for bc in barcodes]
ldata10.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_3_1' for bc in barcodes]
ldata11.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata12.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_4_1' for bc in barcodes]
ldata12.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata13.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_5_1' for bc in barcodes]
ldata13.obs.index = barcodes


barcodes = [bc.split(':')[1] for bc in ldata14.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_6_2' for bc in barcodes]
ldata14.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata15.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_24_2' for bc in barcodes]
ldata15.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata16.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_15_2' for bc in barcodes]
ldata16.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata17.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_3_1_1' for bc in barcodes]
ldata17.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata18.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_4_1_1' for bc in barcodes]
ldata18.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata19.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_5_1_1' for bc in barcodes]
ldata19.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata20.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_6_1_1' for bc in barcodes]
ldata20.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata21.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_7_1_1' for bc in barcodes]
ldata21.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata22.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_9_1_1' for bc in barcodes]
ldata22.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata23.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_7_2' for bc in barcodes]
ldata23.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata24.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_25_2' for bc in barcodes]
ldata24.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata25.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_2_1_1' for bc in barcodes]
ldata25.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata26.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_17_2' for bc in barcodes]
ldata26.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata27.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_8_2' for bc in barcodes]
ldata27.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata28.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_26_2' for bc in barcodes]
ldata28.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata29.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_18_2' for bc in barcodes]
ldata29.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata30.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_8_1_1' for bc in barcodes]
ldata30.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata31.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_19_2' for bc in barcodes]
ldata31.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata32.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_1_1_1' for bc in barcodes]
ldata32.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata33.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_9_2' for bc in barcodes]
ldata33.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata34.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_27_2' for bc in barcodes]
ldata34.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata35.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_16_2' for bc in barcodes]
ldata35.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata36.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_10_2' for bc in barcodes]
ldata36.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata37.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_28_2' for bc in barcodes]
ldata37.obs.index = barcodes


barcodes = [bc.split(':')[1] for bc in ldata38.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_20_2' for bc in barcodes]
ldata38.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata39.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_21_2' for bc in barcodes]
ldata39.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata40.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_22_2' for bc in barcodes]
ldata40.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata41.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_23_2' for bc in barcodes]
ldata41.obs.index = barcodes


# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
ldata7.var_names_make_unique()
ldata8.var_names_make_unique()
ldata9.var_names_make_unique()
ldata10.var_names_make_unique()
ldata11.var_names_make_unique()
ldata12.var_names_make_unique()
ldata13.var_names_make_unique()
ldata14.var_names_make_unique()
ldata15.var_names_make_unique()
ldata16.var_names_make_unique()
ldata17.var_names_make_unique()
ldata18.var_names_make_unique()
ldata19.var_names_make_unique()
ldata20.var_names_make_unique()
ldata21.var_names_make_unique()
ldata22.var_names_make_unique()
ldata23.var_names_make_unique()
ldata24.var_names_make_unique()
ldata25.var_names_make_unique()
ldata26.var_names_make_unique()
ldata27.var_names_make_unique()
ldata28.var_names_make_unique()
ldata29.var_names_make_unique()
ldata30.var_names_make_unique()
ldata31.var_names_make_unique()
ldata32.var_names_make_unique()
ldata33.var_names_make_unique()
ldata34.var_names_make_unique()
ldata35.var_names_make_unique()
ldata36.var_names_make_unique()
ldata37.var_names_make_unique()
ldata38.var_names_make_unique()
ldata39.var_names_make_unique()
ldata40.var_names_make_unique()
ldata41.var_names_make_unique()

# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3, ldata4, ldata5, ldata6, ldata7, ldata8, ldata9, ldata10, ldata11, ldata12, ldata13, ldata14, ldata15, ldata16, ldata17, ldata18, ldata19, ldata20, ldata21, ldata22, ldata23, ldata24, ldata25, ldata26, ldata27, ldata28, ldata29, ldata30, ldata31, ldata32, ldata33, ldata34, ldata35, ldata36, ldata37, ldata38, ldata39, ldata40, ldata41])
# merge matrices into the original adata object
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='MergedAnnotation', frameon=False, legend_loc='on data', title='Annotation', save='_celltypes.pdf')


sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='Clusters', save='_clusters.pdf')

## Save a checkpoint here
adata.write('Merged_data.h5ad')



##
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('Merged_data.h5ad')

## Collectively, this is both A) too large of an object to work with effectively, and B) maybe shouldn't be treated all together?
## Let's split it by treatment, calculating velocity and making the vector/stream maps separately (although still in the shared UMAP embedding)

## Start with Control:
adata = adata[adata.obs.Treatment == "Control"]

scv.pl.proportions(adata, groupby='MergedAnnotation')

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='Control_embedding.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='Control_embedding_grid.pdf', title='', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='Control_embedding_stream.pdf', title='')

adata.write("Controls_data.h5ad")

scv.tl.rank_velocity_genes(adata, groupby='MergedAnnotation', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.head()
#scv.pl.scatter(adata, df['SATB2'][:3], ylabel='SATB2', frameon=False, color='MergedAnnotation', size=10, linewidth=1.5, save="SATB2_MultiControl.pdf")
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="MultiControl_velLengthConfidence.pdf")
scv.pl.velocity_graph(adata, threshold=.1, color='MergedAnnotation', save="MultiControl_VelGraph.pdf")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="MultiControl_Pseudotime.pdf")

## Need to run recover_dynamics... will take a long time?
#scv.tl.latent_time(adata)
#scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="MultiControl_LatentTime.pdf")

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_subset, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100)


# adata.uns['neighbors']['distances'] = adata.obsp['distances']
# adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
# scv.tl.paga(adata, groups='MergedAnnotation')
# df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
# running PAGA using priors: ['velocity_pseudotime']
# 
# scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
#             min_edge_width=2, node_size_scale=1.5)
#scv.pl.velocity(adata, var_names=['SATB2'], color='MergedAnnotation')



## Moving on to VPA
## Collectively, this is both A) too large of an object to work with effectively, and B) maybe shouldn't be treated all together?
## Let's split it by treatment, calculating velocity and making the vector/stream maps separately (although still in the shared UMAP embedding)
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('Merged_data.h5ad')

adata = adata[adata.obs.Treatment == "VPA"]
adata.write("VPA_data.h5ad")

scv.pl.proportions(adata, groupby='MergedAnnotation', save="VPA_Proportions.pdf")

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='VPA_embedding.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='VPA_embedding_grid.pdf', title='', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='VPA_embedding_stream.pdf', title='')

adata.write("VPA_data.h5ad")















## Oops, this should probably have been done split into Single-donor vs Multi-donor, too...
## Well, let's get going.
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('SingleDonor_Basic.h5ad')

adata = adata[adata.obs.Treatment == "Control"]
#adata.write("SingleDonor_Control_data.h5ad")
scv.pl.proportions(adata, groupby='MergedAnnotation', save="SingleDonor_Control.pdf")
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='SingleDonor_Control_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='SingleDonor_Control_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='SingleDonor_Control_embedding_stream.pdf', title='')
adata.write("SingleDonor_Control_data.h5ad")






## Oops, this should probably have been done split into Single-donor vs Multi-donor, too...
## Well, let's get going.
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('SingleDonor_Basic.h5ad')

adata = adata[adata.obs.Treatment == "VPA"]
#adata.write("SingleDonor_Control_data.h5ad")
scv.pl.proportions(adata, groupby='MergedAnnotation', save="SingleDonor_VPA.pdf")
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='SingleDonor_VPA_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='SingleDonor_VPA_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='SingleDonor_VPA_embedding_stream.pdf', title='')
adata.write("SingleDonor_VPA_data.h5ad")


## Oops, this should probably have been done split into Single-donor vs Multi-donor, too...
## Well, let's get going.
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('MultiDonor_Basic.h5ad')

adata = adata[adata.obs.Treatment == "Control"]
#adata.write("MultiDonor_Control_data.h5ad")
scv.pl.proportions(adata, groupby='MergedAnnotation', save="MultiDonor_Control.pdf")
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='MultiDonor_Control_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='MultiDonor_Control_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='MultiDonor_Control_embedding_stream.pdf', title='')
adata.write("MultiDonor_Control_data.h5ad")






## Oops, this should probably have been done split into Single-donor vs Multi-donor, too...
## Well, let's get going.
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('MultiDonor_Basic.h5ad')

adata = adata[adata.obs.Treatment == "VPA"]
#adata.write("MultiDonor_Control_data.h5ad")
scv.pl.proportions(adata, groupby='MergedAnnotation', save="MultiDonor_VPA.pdf")
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='MultiDonor_VPA_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='MergedAnnotation', save='MultiDonor_VPA_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='MultiDonor_VPA_embedding_stream.pdf', title='')
adata.write("MultiDonor_VPA_data.h5ad")
scv.tl.rank_velocity_genes(adata, groupby='MergedAnnotation', min_corr=.3)
#df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
#df.head()
#scv.pl.scatter(adata, df['SATB2'][:3], ylabel='SATB2', frameon=False, color='MergedAnnotation', size=10, linewidth=1.5, save="SATB2_MultiControl.pdf")
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="MultiVPA_velLengthConfidence.pdf")
scv.pl.velocity_graph(adata, threshold=.1, color='MergedAnnotation', save="MultiVPA_VelGraph.pdf")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="MultiVPA_Pseudotime.pdf")
adata.write("MultiDonor_VPA_data_GeneRank_Pseudotime.h5ad")





import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('MultiDonor_Control_data.h5ad')
scv.tl.rank_velocity_genes(adata, groupby='MergedAnnotation', min_corr=.3)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="MultiControl_velLengthConfidence.pdf")
scv.pl.velocity_graph(adata, threshold=.1, color='MergedAnnotation', save="MultiControl_VelGraph.pdf")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="MultiControl_Pseudotime.pdf")
adata.write("MultiDonor_Control_data_GeneRank_Pseudotime.h5ad")


import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('SingleDonor_Control_data.h5ad')
scv.tl.rank_velocity_genes(adata, groupby='MergedAnnotation', min_corr=.3)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="SingleControl_velLengthConfidence.pdf")
scv.pl.velocity_graph(adata, threshold=.1, color='MergedAnnotation', save="SingleControl_VelGraph.pdf")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="SingleControl_Pseudotime.pdf")
adata.write("SingleDonor_Control_data_GeneRank_Pseudotime.h5ad")

adata = sc.read_h5ad('SingleDonor_VPA_data.h5ad')
scv.tl.rank_velocity_genes(adata, groupby='MergedAnnotation', min_corr=.3)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="SingleVPA_velLengthConfidence.pdf")
scv.pl.velocity_graph(adata, threshold=.1, color='MergedAnnotation', save="SingleVPA_VelGraph.pdf")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="SingleVPA_Pseudotime.pdf")
adata.write("SingleDonor_VPA_data_GeneRank_Pseudotime.h5ad")


## "UPSET" package! Give it many groups, it compares all of them, shows the intersections between groups!



## Begin CellRank analysis?
adata = sc.read_h5ad("SingleDonor_Control_data_GeneRank_Pseudotime.h5ad")

adataH1 = adata[adata.obs["Genotype"].isin(["H1"])]
adataPGP1 = adata[adata.obs["Genotype"].isin(["PGP1"])]
adataCW = adata[adata.obs["Genotype"].isin(["CW"])]
adata11a = adata[adata.obs["Genotype"].isin(["11a"])]
adataMito = adata[adata.obs["Genotype"].isin(["Mito210"])]

sc.pp.neighbors(adataH1, n_neighbors=15, use_rep='X_pca')
scv.pp.filter_and_normalize(adataH1)
scv.pp.moments(adataH1)




## Latent time (and trying out the dyanmical model instead of the stochastic?)
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('SingleDonor_Basic.h5ad')
adata = adata[adata.obs.Treatment == "Control"]
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
print("Done loading data. Begin recovering dynamics.")
scv.tl.recover_dynamics(adata)
print("Done recovering dynamics. Begin dynamic velocity algorithm.")
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='SD_Control_embedding_stream.pdf', title='')
print("Done with dynamic velocity. Begin latent time calculations.")
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="SD_Control_LatentTime.pdf")
print("Done with latent time. Saving object.")
adata.write("SD_Control_LatentTime_data.h5ad")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes[:50], sortby='latent_time', col_color='MergedAnnotation', n_convolve=500, save="SD_Control_GeneHeatmapLT.pdf")


## Mem tea, tea infusor, etc


## Look into the DEGs between Control/VPA for each donor line in each protocol, and go for that correlation stuff.
## 

## Latent Time Multidonor Control
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('MultiDonor_Basic.h5ad')
adata = adata[adata.obs.Treatment == "Control"]
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
print("Done loading data. Begin recovering dynamics.")
scv.tl.recover_dynamics(adata)
print("Done recovering dynamics. Begin dynamic velocity algorithm.")
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='MD_Control_embedding_stream.pdf', title='')
print("Done with dynamic velocity. Begin latent time calculations.")
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="MD_Control_LatentTime.pdf")
print("Done with latent time. Saving object.")
adata.write("MD_Control_LatentTime_data.h5ad")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes[:50], sortby='latent_time', col_color='MergedAnnotation', n_convolve=500, save="MD_Control_GeneHeatmapLT.pdf")

##VPA Latent Time
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('SingleDonor_Basic.h5ad')
adata = adata[adata.obs.Treatment == "VPA"]
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
print("Done loading data. Begin recovering dynamics.")
scv.tl.recover_dynamics(adata)
print("Done recovering dynamics. Begin dynamic velocity algorithm.")
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='SD_VPA_embedding_stream.pdf', title='')
print("Done with dynamic velocity. Begin latent time calculations.")
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="SD_VPA_LatentTime.pdf")
print("Done with latent time. Saving object.")
adata.write("SD_VPA_LatentTime_data.h5ad")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes[:50], sortby='latent_time', col_color='MergedAnnotation', n_convolve=500, save="SD_VPA_GeneHeatmapLT.pdf")



## Latent Time Multidonor VPA
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('MultiDonor_Basic.h5ad')
adata = adata[adata.obs.Treatment == "VPA"]
scv.pp.filter_and_normalize(adata)
scv.pp.neighbors(adata)
scv.pp.moments(adata)
print("Done loading data. Begin recovering dynamics.")
scv.tl.recover_dynamics(adata)
print("Done recovering dynamics. Begin dynamic velocity algorithm.")
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['MergedAnnotation', 'Batch'], save='MD_VPA_embedding_stream.pdf', title='')
print("Done with dynamic velocity. Begin latent time calculations.")
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="MD_VPA_LatentTime.pdf")
print("Done with latent time. Saving object.")
adata.write("MD_VPA_LatentTime_data.h5ad")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes[:50], sortby='latent_time', col_color='MergedAnnotation', n_convolve=500, save="MD_VPA_GeneHeatmapLT.pdf")







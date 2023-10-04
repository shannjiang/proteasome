library(Seurat)
library(SeuratDisk)
library(plyr)
#library(pheatmap)
library(data.table)

work_dir = '/home/shann/Documents/synapse/MIT_ROSMAP_Multiomics/'
out_dir = '/home/shann/Documents/Natura/proteosome/MIT_ROSMAP_multiomics/'
rosmap_meta = read.csv(file = '/home/shann/Documents/synapse/ROSMAP/meta/ROSMAP_clinical_metadata.csv', header = T)

celltype = 'Immune_cells'
seu = readRDS(paste0(work_dir,celltype,'.rds'))
meta_df = seu@meta.data
meta_df2 = join(meta_df,rosmap_meta,by = 'projid')
rownames(meta_df2) = rownames(meta_df)
cts = as.data.frame(t(as.matrix(seu@assays$RNA@counts)))

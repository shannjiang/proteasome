library(Seurat)
library(SeuratDisk)
library(plyr)
#library(pheatmap)
library(data.table)

work_dir = '/home/shann/Documents/synapse/MIT_ROSMAP_Multiomics/'
out_dir = '/home/shann/Documents/Natura/proteosome/MIT_ROSMAP_multiomics/'
rosmap_meta = read.csv(file = '/home/shann/Documents/synapse/ROSMAP/meta/ROSMAP_clinical_metadata.csv', header = T)

celltype = 'Immune_cells'
celltype = 'Vasculature_cells'
seu = readRDS(paste0(work_dir,celltype,'.rds'))
meta_df = seu@meta.data
meta_df2 = join(meta_df,rosmap_meta,by = 'projid')
rownames(meta_df2) = rownames(meta_df)
cts = as.data.frame(t(as.matrix(seu@assays$RNA@counts)))

if(celltype == 'Immune_cells'){
  celltype_df = data.frame(broadClass = c('CAMs',rep('Microglia',3),'T_cells'), subClass = c('CAMs','Mic MKI67','Mic P2RY12','Mic TPT1','T cells'))
  celltypes = c('CAMs','Microglia','T_cells')
  for (celltype in celltypes){
    target_celltype_df = celltype_df[celltype_df$broadClass %in% celltype,]
    target_meta_df = meta_df2[meta_df2$cell_type_high_resolution %in% target_celltype_df$subClass,]
    target_cts = cts[rownames(cts) %in% rownames(target_meta_df),]
    target_cts$projid = target_meta_df$projid
    pseudoBulk_cts = aggregate(.~projid,target_cts,sum)
    rownames(pseudoBulk_cts) = pseudoBulk_cts$projid
    pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'projid']
    pseudoBulk_meta = rosmap_meta[rosmap_meta$projid %in% rownames(pseudoBulk_cts),]
    rownames(pseudoBulk_meta) = pseudoBulk_meta$projid
    pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]
    save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(out_dir,celltype,'_pseudoBulk.RData'))
    remove(target_cts,cts,target_meta_df)
  }
}

if(celltype == 'Vasculature_cells'){
  celltype_df = data.frame(broadClass = c('End',rep('Fib',2),'Per','SMC'), subClass = c('End','Fib FLRT2','Fib SLC4A4','Per','SMC'))
  celltypes = c('End')
  for (celltype in celltypes){
    target_celltype_df = celltype_df[celltype_df$broadClass %in% celltype,]
    target_meta_df = meta_df2[meta_df2$cell_type_high_resolution %in% target_celltype_df$subClass,]
    target_cts = cts[rownames(cts) %in% rownames(target_meta_df),]
    target_cts$projid = target_meta_df$projid
    pseudoBulk_cts = aggregate(.~projid,target_cts,sum)
    rownames(pseudoBulk_cts) = pseudoBulk_cts$projid
    pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'projid']
    pseudoBulk_meta = rosmap_meta[rosmap_meta$projid %in% rownames(pseudoBulk_cts),]
    rownames(pseudoBulk_meta) = pseudoBulk_meta$projid
    pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]
    save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(out_dir,celltype,'_pseudoBulk.RData'))
    remove(target_cts,cts,target_meta_df)
  }
}

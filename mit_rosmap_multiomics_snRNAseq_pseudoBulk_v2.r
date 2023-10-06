library(Seurat)
#library(SeuratDisk)
library(plyr)
#library(pheatmap)
library(data.table)

work_dir = '/home/shann/Documents/synapse/MIT_ROSMAP_Multiomics/'
out_dir = '/home/shann/Documents/synapse/MIT_ROSMAP_Multiomics/pseudoBulk/'
if(!dir.exists(out_dir)){dir.create(out_dir,recursive = T)}
rosmap_meta = read.csv(file = '/home/shann/Documents/synapse/ROSMAP/meta/ROSMAP_clinical_metadata.csv', header = T)
subset_size = 2500

#helper
#seu_obj = seu_list[[1]]
seu2pseudoBulk = function(seu_obj){
  meta_df = seu_obj@meta.data
  meta_df2 = join(meta_df,rosmap_meta,by = 'projid')
  rownames(meta_df2) = rownames(meta_df)
  cts = as.data.frame(t(as.matrix(seu_obj@assays$RNA@counts)))
  cts$projid = meta_df2$projid
  aggCts = aggregate(.~projid,cts,sum)
  return(aggCts)
}


#celltype = 'Immune_cells'
#Immune_cells sub celltypes CAMs Mic MKI67 Mic P2RY12 Mic TPT1 T cells
celltype = 'Vasculature_cells'
#Vasculature cells sub celltypes End Fib FLRT2 Fib SLC4A4 Per SMC
for(celltype in celltypes){
seu = readRDS(paste0(work_dir,celltype,'.rds'))
#screen for subset of seu based on cell type high resolution
#tmp_seu = seu[,rownames(seu@meta.data)[seu@meta.data$cell_type_high_resolution %in% c('Fib FLRT2','Fib SLC4A4')]]
#meta_df = seu@meta.data
#meta_df2 = join(meta_df,rosmap_meta,by = 'projid')
#rownames(meta_df2) = rownames(meta_df)
#cts = as.data.frame(t(as.matrix(seu@assays$RNA@counts)))
if(celltype == 'Immune_cells'){
  celltype_df = data.frame(broadClass = c('CAMs',rep('Microglia',3),'T_cells'), subClass = c('CAMs','Mic MKI67','Mic P2RY12','Mic TPT1','T cells'))
  celltypes = c('CAMs','Microglia','T_cells')
  for (target_celltype in celltypes){
    target_celltype_df = celltype_df[celltype_df$broadClass %in% target_celltype,]
    target_seu = seu[,rownames(seu@meta.data)[seu@meta.data$cell_type_high_resolution %in% target_celltype_df$subClass]]
    #subset target_seu each with subset_size cells
    sec_num = ceiling(dim(target_seu)[2]/subset_size)
    seu_list = list()
    for(i in 1:(sec_num-1)){
      #print(paste0('process ',i))
      seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(i-1)):(subset_size*i)]
    }
    seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(sec_num-1)):dim(target_seu)[2]]
    aggCts_list = lapply(c(1:length(seu_list)),function(x){seu2pseudoBulk(seu_list[[x]])})
    aggCts = do.call(rbind,aggCts_list)
    pseudoBulk_cts = aggregate(.~projid,aggCts,sum)
    rownames(pseudoBulk_cts) = pseudoBulk_cts$projid
    pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'projid']
    pseudoBulk_meta = rosmap_meta[rosmap_meta$projid %in% rownames(pseudoBulk_cts),]
    rownames(pseudoBulk_meta) = pseudoBulk_meta$projid
    pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]
    save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(out_dir,target_celltype,'_pseudoBulk.RData'))
    remove(seu_list,target_seu,aggCts_list)
  }
} else if(celltype == 'Vasculature_cells'){
  celltype_df = data.frame(broadClass = c('End',rep('Fib',2),'Per','SMC'), subClass = c('End','Fib FLRT2','Fib SLC4A4','Per','SMC'))
  celltypes = c('End')
  for (target_celltype in celltypes){
    target_celltype_df = celltype_df[celltype_df$broadClass %in% target_celltype,]
    target_seu = seu[,rownames(seu@meta.data)[seu@meta.data$cell_type_high_resolution %in% target_celltype_df$subClass]]
    #subset target_seu each with subset_size cells
    sec_num = ceiling(dim(target_seu)[2]/subset_size)
    seu_list = list()
    for(i in 1:(sec_num-1)){
      #print(paste0('process ',i))
      seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(i-1)):(subset_size*i)]
    }
    seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(sec_num-1)):dim(target_seu)[2]]
    aggCts_list = lapply(c(1:length(seu_list)),function(x){seu2pseudoBulk(seu_list[[x]])})
    aggCts = do.call(rbind,aggCts_list)
    pseudoBulk_cts = aggregate(.~projid,aggCts,sum)
    rownames(pseudoBulk_cts) = pseudoBulk_cts$projid
    pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'projid']
    pseudoBulk_meta = rosmap_meta[rosmap_meta$projid %in% rownames(pseudoBulk_cts),]
    rownames(pseudoBulk_meta) = pseudoBulk_meta$projid
    pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]
    save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(out_dir,target_celltype,'_pseudoBulk.RData'))
    remove(seu_list,target_seu,aggCts_list)
  }
} else {
  target_seu = seu
  #subset target_seu each with subset_size cells
  sec_num = ceiling(dim(target_seu)[2]/subset_size)
  seu_list = list()
  for(i in 1:(sec_num-1)){
    #print(paste0('process ',i))
    seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(i-1)):(subset_size*i)]
  }
  seu_list[length(seu_list)+1] = target_seu[,(1+subset_size*(sec_num-1)):dim(target_seu)[2]]
  aggCts_list = lapply(c(1:length(seu_list)),function(x){seu2pseudoBulk(seu_list[[x]])})
  aggCts = do.call(rbind,aggCts_list)
  pseudoBulk_cts = aggregate(.~projid,aggCts,sum)
  rownames(pseudoBulk_cts) = pseudoBulk_cts$projid
  pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'projid']
  pseudoBulk_meta = rosmap_meta[rosmap_meta$projid %in% rownames(pseudoBulk_cts),]
  rownames(pseudoBulk_meta) = pseudoBulk_meta$projid
  pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]
  save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(out_dir,celltype,'_pseudoBulk.RData'))
  remove(seu_list,target_seu,aggCts_list)
}
}

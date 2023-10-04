library(Seurat)
library(SeuratDisk)
library(plyr)
library(pheatmap)
library(data.table)
library(ggplot2)

norm_expr_dir = '/home/shann/Documents/Natura/proteosome/SEAAD/normalized_expr/'
meta_dir = '/home/shann/Documents/Natura/proteosome/SEAAD/normalized_expr/'
graph_dir = '/home/shann/Documents/Natura/proteosome/SEAAD/graph/'
if(!dir.exists(graph_dir)){dir.create(graph_dir, recursive = T)}
prot_gene_id_file = '/home/shann/Documents/Natura/proteosome/gene_id.txt'

#helper function
exp_pct_cal = function(vec){
  return(sum(vec != 0)/length(vec))
}

sem_cal = function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

Zscore_normalization = function(df){
  impute_df = as.data.frame(do.call(cbind, lapply(c(1:length(colnames(df))),function(x){df[,x][is.na(df[,x])] = mean(df[,x],na.rm=T);return(df[,x])})))
  norm_df = as.data.frame(do.call(cbind, lapply(c(1:length(colnames(impute_df))),function(x){(impute_df[,x]-mean(impute_df[,x]))/sd(impute_df[,x])})))
  colnames(norm_df) = colnames(df)
  rownames(norm_df) = rownames(df)
  return(norm_df)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

prot_gene_by_braak_heatmap_df_gen = function(braak_aggregated_df,gene_substructure_df){
  rownames(braak_aggregated_df) = braak_aggregated_df$braak
  braak_aggregated_df = braak_aggregated_df[,!colnames(braak_aggregated_df) %in% 'braak']
  braak_aggregated_df_norm = Zscore_normalization(braak_aggregated_df)
  braak_aggregated_df_norm_t = as.data.frame(t(braak_aggregated_df_norm))
  gene_substructure_df = gene_substructure_df[gene_substructure_df$gene %in% intersect(gene_substructure_df$gene,rownames(braak_aggregated_df_norm_t)),]
  rownames(gene_substructure_df) = gene_substructure_df$gene
  braak_aggregated_df_norm_t = braak_aggregated_df_norm_t[rownames(braak_aggregated_df_norm_t) %in% gene_substructure_df$gene,]
  braak_aggregated_df_norm_t = braak_aggregated_df_norm_t[rownames(gene_substructure_df),]
  colnames(braak_aggregated_df_norm_t) = paste0('braak',c(0,2:6))
  return(braak_aggregated_df_norm_t)
}

log10_by_braak_df_gen = function(celltype = celltype, brainRegion = brainRegion) {
  norm_expr_df = read.csv(file = paste0(norm_expr_dir,'prot_',celltype,'_',brainRegion,'_normDat.csv'), header = T, row.names = 1)
  log10_expr_df = log1p(norm_expr_df)
  gene_annot_df = read.csv(file = paste0(meta_dir,celltype,'_',brainRegion,'_gene_annot.csv'), header = T, row.names = 1)
  gene_annot_df = gene_annot_df[rownames(gene_annot_df) %in% colnames(log10_expr_df),]
  gene_annot_df = gene_annot_df[colnames(log10_expr_df),]
  colnames(log10_expr_df) = gene_annot_df$feature_name
  meta_df = read.csv(file = paste0(meta_dir,celltype,'_',brainRegion,'_meta.csv'), header = T, row.names = 1)
  meta_df = meta_df[rownames(log10_expr_df),]
  meta_df$braak = with(meta_df,ifelse(Braak.stage == 'Braak 0','braak0',ifelse(Braak.stage == 'Braak II','braak2',ifelse(Braak.stage == 'Braak III','braak3',ifelse(Braak.stage == 'Braak IV','braak4', ifelse(Braak.stage == 'Braak V','braak5',ifelse(Braak.stage == 'Braak VI','braak6','reference')))))))
  log10_expr_df$braak = meta_df$braak
  log10_expr_df = log10_expr_df[!log10_expr_df$braak %in% 'reference',]
  return(log10_expr_df)
}

#proteasome_gene_annotation
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1','MAFF','MAFG','MAFK')
proteasome_tf_genes = c('NFE2L1')
proteasome_receptor_genes = c('ADRM1')

substructure = rep(c('19S','20S','Immunoproteasome','Chaperon','Activator','TF','Receptor'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes),length(proteasome_receptor_genes)))
gene = c(proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes,proteasome_tf_genes,proteasome_receptor_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene
gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% 'Receptor',]

proteasome_gene_id_df = read.table(file = prot_gene_id_file, header = T, sep = '\t')
rownames(proteasome_gene_id_df) = proteasome_gene_id_df$gene
proteasome_gene_id_df$gene = proteasome_gene_id_df$OGS
gene_substructure_df = join(gene_substructure_df,proteasome_gene_id_df,by = 'gene')
rownames(gene_substructure_df) = gene_substructure_df$OGS

subClasses = c('L23_IT','Astrocyte','Chandelier','Endothelial','L4_IT','L5_ET','L5_IT','L56_NP','L6b','L6_CT','L6_IT_Car3','L6_IT','Lamp5','Lamp5_Lhx6','Microglia_PVM','Oligodendrocyte','OPC','Pax6','Pvalb','Sncg','Sst_Chodl','Sst','Vip','VLMC')
broadClasses = c('Glutamatergic','Astrocyte','GABAergic','Endothelial','Glutamatergic','Glutamatergic','Glutamatergic','Glutamatergic','Glutamatergic','Glutamatergic','Glutamatergic','Glutamatergic','GABAergic','GABAergic','Microglia','Oligodendrocyte','OPC','GABAergic','GABAergic','GABAergic','GABAergic','GABAergic','GABAergic','VLMC')
brainRegions = c('MTG','DLPFC')

celltype_df = data.frame(broadClass = broadClasses, subClass = subClasses)
celltypes = unique(broadClasses)

celltype = 'Microglia'
brainRegion = 'DLPFC'
for (brainRegion in brainRegions){
for (celltype in celltypes){
targetCelltype_df = celltype_df[celltype_df$broadClass %in% celltype,]
log10_by_braak_list = lapply(targetCelltype_df$subClass, function(x){log10_by_braak_df_gen(x,brainRegion)})
log10_by_braak_df = do.call(rbind,log10_by_braak_list)
log10_by_braak_mean_df = aggregate(.~braak,log10_by_braak_df,mean)

prot_expr_pct_by_braak_norm_t = prot_gene_by_braak_heatmap_df_gen(braak_aggregated_df = log10_by_braak_mean_df, gene_substructure_df = gene_substructure_df)

col_annot_df = data.frame(braak = paste0('braak',c(0,2:6)))
rownames(col_annot_df) = paste0('braak',c(0,2:6))
row_annot_df = data.frame(substructure = gene_substructure_df$substructure)
rownames(row_annot_df) = rownames(gene_substructure_df)

fig = pheatmap(prot_expr_pct_by_braak_norm_t, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, fontsize_row = 6, annotation_col = col_annot_df, annotation_row = row_annot_df)
save_pheatmap_pdf(fig, paste0(graph_dir,'SEAAD_',celltype,'_',brainRegion,'_gene_expr_pct_by_braak_stage_heatmap.pdf'), width = 3.7, height = 7)
}}

library(Seurat)
library(SeuratDisk)
library(plyr)
library(pheatmap)
library(data.table)
library(ggplot2)

work_dir = '/home/shann/Documents/Natura/proteosome/vilas_1m/'
out_dir = work_dir
proteasome_seu_dir = paste0(out_dir,'proteasome_seurat/')
if(!dir.exists(proteasome_seu_dir)){dir.create(proteasome_seu_dir, recursive = T)}
graph_dir = paste0(out_dir,'graphs/')
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
  colnames(braak_aggregated_df_norm_t) = paste0('braak',c(0:6))
  return(braak_aggregated_df_norm_t)
}

prot_subunit_expr_pct_by_braak_barplot = function(prot_cts, substructure_annot, braak_annot){
  prot_cts = prot_cts[gene_substructure_df$gene,]
  prot_cts$substructure = substructure_annot
  subunit_cts = aggregate(.~substructure,prot_cts,sum)
  rownames(subunit_cts) = subunit_cts$substructure
  subunit_cts = subunit_cts[,-1]
  subunit_cts_t = as.data.frame(t(subunit_cts))
  subunit_cts_t$braak = braak_annot
  subunit_cts_t = subunit_cts_t[!subunit_cts_t$braak %in% 'braakNA',]
  subunit_cts_expr_pct_by_braak = aggregate(.~braak,subunit_cts_t,exp_pct_cal)
  subunit_cts_expr_pct_by_braak_long_dt = melt(data = as.data.table(subunit_cts_expr_pct_by_braak), id.vars = 'braak', variable.name = 'subunit', value.name = 'expression_pct')

  p<- ggplot(subunit_cts_expr_pct_by_braak_long_dt, aes(x=subunit, y=expression_pct, fill=braak)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge())
  return(p)
}

prot_subunit_log10_by_braak_barplot = function(prot_cts, substructure_annot, braak_annot){
  prot_cts = prot_cts[gene_substructure_df$gene,]
  prot_cts$substructure = substructure_annot
  subunit_cts = aggregate(.~substructure,prot_cts,sum)
  rownames(subunit_cts) = subunit_cts$substructure
  subunit_cts = subunit_cts[,-1]
  subunit_log10 = log10(subunit_cts + 1)
  subunit_log10_t = as.data.frame(t(subunit_log10))
  subunit_log10_t$braak = braak_annot
  subunit_log10_t = subunit_log10_t[!subunit_log10_t$braak %in% 'braakNA',]
  subunit_log10_mean_by_braak = aggregate(.~braak,subunit_log10_t,mean)
  subunit_log10_sem_by_braak = aggregate(.~braak,subunit_log10_t,sem_cal)
  subunit_log10_mean_by_braak_long_dt = melt(data = as.data.table(subunit_log10_mean_by_braak), id.vars = 'braak', variable.name = 'subunit', value.name = 'mean')
  subunit_log10_sem_by_braak_long_dt = melt(data = as.data.table(subunit_log10_sem_by_braak), id.vars = 'braak', variable.name = 'subunit', value.name = 'sem')
  subunit_log10_by_braak_long_dt = subunit_log10_mean_by_braak_long_dt
  subunit_log10_by_braak_long_dt$sem = subunit_log10_sem_by_braak_long_dt$sem
  p<- ggplot(subunit_log10_by_braak_long_dt, aes(x=subunit, y=mean, fill=braak)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                   position=position_dodge(.9)) +
    ylab('log10(cts + 1)')
  return(p)
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

###########celltype-specific heatmap and barplot
celltype = 'endo'

seu <- LoadH5Seurat(paste0(work_dir,celltype,'/',celltype,'.h5Seurat'))
prot_seu = seu[rownames(seu) %in% gene_substructure_df$gene,]
SaveH5Seurat(prot_seu, filename = paste0(proteasome_seu_dir,'prot_',celltype,'_.h5Seurat'))

#gene heatmap expr_pct
prot_cts = as.data.frame(prot_seu@assays$RNA@counts)
meta_df = prot_seu@meta.data
prot_cts_t = as.data.frame(t(prot_cts))
prot_cts_t$braak = paste0('braak',meta_df$braaksc)
prot_cts_t = prot_cts_t[!prot_cts_t$braak %in% 'braakNA',]
prot_expr_pct_by_braak_df = aggregate(.~braak, prot_cts_t,exp_pct_cal)

prot_expr_pct_by_braak_norm_t = prot_gene_by_braak_heatmap_df_gen(braak_aggregated_df = prot_expr_pct_by_braak_df, gene_substructure_df = gene_substructure_df)

col_annot_df = data.frame(braak = paste0('braak',c(0:6)))
rownames(col_annot_df) = paste0('braak',c(0:6))
row_annot_df = data.frame(substructure = gene_substructure_df$substructure)
rownames(row_annot_df) = rownames(gene_substructure_df)

fig = pheatmap(prot_expr_pct_by_braak_norm_t, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, fontsize_row = 6, annotation_col = col_annot_df, annotation_row = row_annot_df)
save_pheatmap_pdf(fig, paste0(graph_dir,'vilas_1m_',celltype,'_gene_expr_pct_by_braak_stage_heatmap.pdf'), width = 3.7, height = 7)

#gene heatmap long10(cts + 1)
prot_cts = as.data.frame(prot_seu@assays$RNA@counts)
meta_df = prot_seu@meta.data
prot_cts_t = as.data.frame(t(prot_cts))
prot_log10_t = log10(prot_cts_t + 1)
prot_log10_t$braak = paste0('braak',meta_df$braaksc)
prot_log10_t = prot_log10_t[!prot_log10_t$braak %in% 'braakNA',]
prot_log10_by_braak_df = aggregate(.~braak, prot_log10_t, mean)

prot_log10_by_braak_norm_t = prot_gene_by_braak_heatmap_df_gen(braak_aggregated_df = prot_log10_by_braak_df, gene_substructure_df = gene_substructure_df)

col_annot_df = data.frame(braak = paste0('braak',c(0:6)))
rownames(col_annot_df) = paste0('braak',c(0:6))
row_annot_df = data.frame(substructure = gene_substructure_df$substructure)
rownames(row_annot_df) = rownames(gene_substructure_df)

fig = pheatmap(prot_log10_by_braak_norm_t, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, fontsize_row = 6, annotation_col = col_annot_df, annotation_row = row_annot_df)
save_pheatmap_pdf(fig, paste0(graph_dir,'vilas_1m_',celltype,'_gene_log10_by_braak_stage_heatmap.pdf'), width = 3.7, height = 7)

#subunit barplot expr_pct
prot_cts = as.data.frame(prot_seu@assays$RNA@counts)
p = prot_subunit_expr_pct_by_braak_barplot(prot_cts = prot_cts, substructure_annot = gene_substructure_df$substructure, braak_annot = paste0('braak',meta_df$braaksc))
pdf(file = paste0(graph_dir,'vilas_1m_',celltype,'_subunit_expr_pct_by_braak_stage_barplot.pdf'), width = 8, height = 7)
p
dev.off()

#subunit barplot log10(cts + 1)
prot_cts = as.data.frame(prot_seu@assays$RNA@counts)
p = prot_subunit_log10_by_braak_barplot(prot_cts = prot_cts, substructure_annot = gene_substructure_df$substructure, braak_annot = paste0('braak',meta_df$braaksc))
pdf(file = paste0(graph_dir,'vilas_1m_',celltype,'_subunit_log10_by_braak_stage_barplot.pdf'), width = 8, height = 7)
p
dev.off()

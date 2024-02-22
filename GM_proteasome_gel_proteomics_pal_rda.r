library(DEP)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

work_dir = 'C:/Users/shann/Documents/Natura/proteasome_project/proteomics/gel_proteomics/files/'
z = load(file = paste0(work_dir,'data_intensity_norm_imp_thr9_jan2022.rda'))

#helper function
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

proteasome_gene_df = read.table(file = paste0('C:/Users/shann/Documents/Natura/', 'gene_id.txt'), header = T, sep = '\t')

############# proteasome and related gene heatmap
#proteasome_gene_annotation
e_ligase_genes = c('UFC1','UBE2O','TRIM2','HECTD3','HUWE1','RBX1')
dubs_deubiquitylating_enzyme_genes = c('USP39','USP14','UCHL1','EIF3F','EIF3H','COPS5','COPS6','VCP')
hsp_protein_chaperon_genes = c('HSPB1','HSPA4','HSPA4L','HSPA2','HSPA8','HSP90AA1','HSP90AB1','CDC37','HSPH1')
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1')
#proteasome_receptor_genes = c('ADRM1')
substrLevel = c('E_ligase','DUBs_deubiquitylating_enzyme','HSP_protein_chaperon','P19S','P20S','IP','Assembly_chaperon','Activator')

substructure = rep(substrLevel,c(length(e_ligase_genes),length(dubs_deubiquitylating_enzyme_genes),length(hsp_protein_chaperon_genes),length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes)))
gene = c(e_ligase_genes,dubs_deubiquitylating_enzyme_genes,hsp_protein_chaperon_genes,proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene

datExpr = data_intensity_norm_imp
colnames(datExpr) = gsub('Control','CTR',colnames(datExpr))
sharedGenes = intersect(rownames(datExpr),rownames(gene_substructure_df))
row_annot_df = gene_substructure_df[rownames(gene_substructure_df) %in% sharedGenes,]
row_annot_df = row_annot_df[c('substructure')]
row_annot_df$substructure = factor(row_annot_df$substructure, levels = substrLevel)
datExpr2 = datExpr[rownames(datExpr) %in% sharedGenes,]
datExpr2 = datExpr2[rownames(row_annot_df),]
col_annot_df = data.frame(diagnosis = c(rep('CTR',12),rep('AD',14)))
rownames(col_annot_df) = c(paste0('CTR_',c(1:12)),paste0('AD_',c(1:14)))
datExpr2 = datExpr2[,rownames(col_annot_df)]

annot_colors = list(condition = c(AD = "#1B9E77", CTR = "#D95F02"),substructure = c(E_ligase = '#8DD3C7',DUBs_deubiquitylating_enzyme = '#FFFFB3',HSP_protein_chaperon = '#BEBADA',P19S = '#FB8072',P20S = '#80B1D3',IP = '#FDB462',Assembly_chaperon = '#B3DE69',Activator = '#FCCDE5'))
fig = pheatmap(datExpr2,color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), cellheight = 8, cellwidth = 8, scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, fontsize = 10,annotation_colors = annot_colors, annotation_col = col_annot_df, annotation_row = row_annot_df)
save_pheatmap_pdf(fig, paste0(work_dir,'gel_proteomics_GM_proteasome_gene_heatmap_v2.pdf'), width = 8, height = 10)


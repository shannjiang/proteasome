library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(pheatmap)

work_dir = '/home/shann/Documents/Natura/proteosome/vilas_1m/proteasome_seurat/'
out_dir = work_dir

#proteasome_gene_annotation
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1','MAFF','MAFG','MAFK')
proteasome_tf_genes = c('NFE2L1','IRF1','STAT1','EGR1')
proteasome_receptor_genes = c('ADRM1')
housekeeper_genes = c('ACTB','GAPDH')

substructure = rep(c('19S','20S','Immunoproteasome','Chaperon','Activator','TF','Receptor','housekeeper'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes),length(proteasome_receptor_genes),length(housekeeper_genes)))
gene = c(proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes,proteasome_tf_genes,proteasome_receptor_genes,housekeeper_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene
#gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% 'Receptor',]
#rownames(gene_substructure_df) = gene_substructure_df$OGS


celltype = 'endo'
seu = LoadH5Seurat(paste0(work_dir,'prot_',celltype,'.h5Seurat'))
meta_df = seu@meta.data
meta_df2 = meta_df[,colnames(meta_df) %in% c('orig.ident','projid','sex','pmAD','Cdx','age_death','apoe_genotype','cogdx_stroke','stroke_bi','braaksc')]
cts = as.data.frame(t(as.matrix(seu@assays$SCT@counts)))
cts$orig.ident = meta_df2$orig.ident
pseudoBulk_cts = aggregate(.~orig.ident, cts, sum)
pseudoBulk_meta = unique(meta_df2)
rownames(pseudoBulk_cts) = pseudoBulk_cts$orig.ident
pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'orig.ident']
rownames(pseudoBulk_meta) = pseudoBulk_meta$orig.ident
pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]

save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(work_dir,celltype,'_pseudoBulk.RData'))

##vst normalization
pseudoBulk_log2 = log2(pseudoBulk_cts+1)
pseudoBulk_log2$braak = 

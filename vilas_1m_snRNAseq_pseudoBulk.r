library(Seurat)
library(SeuratDisk)
library(plyr)
#library(pheatmap)
library(data.table)

work_dir = '/mnt/mfs/ctcn/datasets/rosmap/rnaseq/dlpfcTissue/snrnaseq/annotation.2022-05-10/Seurat/'
out_dir = '/home/sj3195/proteasome/out/snrnaseq/vilas_1m/'
pseudoBulk_dir = paste0(out_dir,'pseudoBulk/')
if(!dir.exists(pseudoBulk_dir)){dir.create(pseudoBulk_dir, recursive = T)}
#out_dir = proteasome_seu_dir

#celltype = 'endo'
#celltypes = c('endo','astrocytes','cux2+','cux2-','inhibitory','microglia','oligodendrocytes','opcs')
celltypes = c('cux2+','cux2-','inhibitory','microglia','oligodendrocytes','opcs')
for(celltype in celltypes) {
seu <- LoadH5Seurat(paste0(work_dir,celltype,'/',celltype,'.h5Seurat'))

#celltype = 'endo'
#seu = LoadH5Seurat(paste0(work_dir,'prot_',celltype,'.h5Seurat'))
meta_df = seu@meta.data
meta_df2 = meta_df[,colnames(meta_df) %in% c('orig.ident','projid','sex','pmAD','Cdx','age_death','apoe_genotype','cogdx_stroke','stroke_bi','braaksc','ceradsc','plaq_d_mf','plaq_n_mf','nft_mf','plaq_d','plaq_n','nft','amyloid','tangles','tdp_stage4','cvda_4gp2','caa_4gp','pmi','cogng_demog_slope','cognep_demog_slope','dlbdx','niareagansc','pAD','sqrt.amyloid','sqrt.tangles','AD_Cdx','sqrt.amy.tan')]
cts = as.data.frame(t(as.matrix(seu@assays$SCT@counts)))
cts$orig.ident = meta_df2$orig.ident
pseudoBulk_cts = aggregate(.~orig.ident, cts, sum)
pseudoBulk_meta = unique(meta_df2)
rownames(pseudoBulk_cts) = pseudoBulk_cts$orig.ident
pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'orig.ident']
rownames(pseudoBulk_meta) = pseudoBulk_meta$orig.ident
pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]

save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(pseudoBulk_dir,celltype,'_pseudoBulk.RData'))
remove(seu,cts,meta_df,meta_df2)
}

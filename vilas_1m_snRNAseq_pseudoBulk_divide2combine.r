library(Seurat)
library(SeuratDisk)
library(plyr)
#library(pheatmap)
library(data.table)

work_dir = '/home/shann/Documents/Natura/proteosome/vilas_1m/endo/'
out_dir = work_dir
pseudoBulk_dir = paste0(out_dir,'pseudoBulk/')
if(!dir.exists(pseudoBulk_dir)){dir.create(pseudoBulk_dir, recursive = T)}
#out_dir = proteasome_seu_dir

celltype = 'endo'
#celltypes = c('endo','astrocytes','cux2+','cux2-','inhibitory','microglia','oligodendrocytes','opcs')
celltypes = c('cux2+','cux2-','inhibitory','microglia','oligodendrocytes','opcs')
for(celltype in celltypes) {
seu <- LoadH5Seurat(paste0(work_dir,celltype,'.h5Seurat'))

#celltype = 'endo'
#seu = LoadH5Seurat(paste0(work_dir,'prot_',celltype,'.h5Seurat'))
#seu = total_seu[,c(1:dim(total_seu)[2]/2)]
meta_df = seu@meta.data
meta_df2 = meta_df[,colnames(meta_df) %in% c('orig.ident','projid','sex','pmAD','Cdx','age_death','apoe_genotype','cogdx_stroke','stroke_bi','braaksc','ceradsc','plaq_d_mf','plaq_n_mf','nft_mf','plaq_d','plaq_n','nft','amyloid','tangles','tdp_stage4','cvda_4gp2','caa_4gp','pmi','cogng_demog_slope','cognep_demog_slope','dlbdx','niareagansc','pAD','sqrt.amyloid','sqrt.tangles','AD_Cdx','sqrt.amy.tan')]
remove(meta_df)
cts = as.data.frame(t(as.matrix(seu@assays$SCT@counts)))
remove(seu)
cts$orig.ident = meta_df2$orig.ident
#divide cts
cts_p1 = cts[c(1:round(dim(cts)[1]/2)),]
cts_p2 = cts[c((round(dim(cts)[1]/2)+1):dim(cts)[1]),]
remove(cts)
pseudoBulk_cts_p1 = aggregate(.~orig.ident, cts_p1, sum)
remove(cts_p1)
pseudoBulk_cts_p2 = aggregate(.~orig.ident, cts_p2, sum)
remove(cts_p2)
pseudoBulk_cts = rbind(pseudoBulk_cts_p1,pseudoBulk_cts_p2)
pseudoBulk_cts = aggregate(.~orig.ident, pseudoBulk_cts, sum)
pseudoBulk_meta = unique(meta_df2)
rownames(pseudoBulk_cts) = pseudoBulk_cts$orig.ident
pseudoBulk_cts = pseudoBulk_cts[,!colnames(pseudoBulk_cts) %in% 'orig.ident']
rownames(pseudoBulk_meta) = pseudoBulk_meta$orig.ident
pseudoBulk_meta = pseudoBulk_meta[rownames(pseudoBulk_cts),]

save(pseudoBulk_cts, pseudoBulk_meta, file = paste0(pseudoBulk_dir,celltype,'_pseudoBulk.RData'))
remove(meta_df2)
}

library(VennDiagram)
library(WebGestaltR)
library(ggplot2)


work_dir = 'C:/Users/shann/Documents/Natura/proteasome_project/proteomics/purified_proteomics/files/'
proteasome_gene_df = read.table(file = paste0('C:/Users/shann/Documents/Natura/', 'gene_id.txt'), header = T, sep = '\t')
gm_res = read.csv(file = paste0(work_dir,'GM_proteinGroups_data_imp_res.csv'), header = T, row.names = 1)
wm_res = read.csv(file = paste0(work_dir,'WM_proteinGroups_data_imp_res.csv'), header = T, row.names = 1)

##Venn diagram
#remove rows with NA for AD_vs_Control_p.val column
#gm_df2 = gm_df[!is.na(gm_df$AD_vs_Control_p.val),]
#wm_df2 = wm_df[!is.na(wm_df$AD_vs_Control_p.val),]
gm_df2 = gm_res
wm_df2 = wm_res
gm_df2$name = toupper(gm_df2$name)
wm_df2$name = toupper(wm_df2$name)
gm_df2$name = gsub('\\..*','',gm_df2$name)
wm_df2$name = gsub('\\..*','',wm_df2$name)
gm_df2 = gm_df2[!gm_df2$name %in% 0,]
gm_df2$name = gsub('_','',gm_df2$name)
wm_df2$name = gsub('_','',wm_df2$name)

gm_genes = gm_df2$name
#gm_genes = gm_genes[-c(1:7)]
wm_genes = wm_df2$name
gene_list = list(gm_genes,wm_genes)
names(gene_list) = c('grey','white')
venn.diagram(gene_list, height = 300, width = 300, fill = c("#FF0000","#00FF00"), lwd = c(1,1), filename = paste0(work_dir,'GMvsWM_venn.tiff'), cex = 0.2, cat.cex = 0.2)

##GM vs WM shared genes p-value Q-Q plot
sharedGenes = intersect(gm_res$name,wm_res$name)
gm_shared_res = gm_res[gm_res$name %in% sharedGenes,]
wm_shared_res = wm_res[wm_res$name %in% sharedGenes,]
shared_res = gm_shared_res[c('name','AD_vs_CTR_p.val')]
colnames(shared_res)[2] = paste0('gm_',colnames(shared_res)[2])
rownames(shared_res) = shared_res$name
rownames(wm_shared_res) = wm_shared_res$name
wm_shared_res = wm_shared_res[rownames(shared_res),]
shared_res$wm_AD_vs_CTR_p.val = wm_shared_res$AD_vs_CTR_p.val
shared_res$gm_log10p = -log10(shared_res$gm_AD_vs_CTR_p.val)
shared_res$wm_log10p = -log10(shared_res$wm_AD_vs_CTR_p.val)
shared_res$color = ifelse(shared_res$name %in% proteasome_gene_df$OGS, 'black', 'grey')
prot_shared_res = 

p = ggplot(shared_res,aes(x=wm_log10p,y=gm_log10p)) + geom_point(color = 'grey',size = 0.5) + geom_point(data = subset(shared_res, shared_res$color %in% 'black'), aes(x = wm_log10p,y = gm_log10p), color = 'black', size = 0.5) + theme(axis.text = element_text(size = 7),axis.title=element_text(size=7),legend.position = "none") + xlab('White Matter -Log10(P value)') + ylab('Grey Matter -Log10(P value)') + geom_abline(intercept = 0, slope = 1, color = 'red')
pdf(file = paste0(work_dir,'GM_vs_WM_QQplot.pdf'), width = 2.4, height = 2.4)
print(p)
dev.off()
p


###########enrichment analysis
#WebGestaltR helper function
WebGestaltR_GO_KEGG = function(projectName = projectName, enrichmentDir = enrichmentDir, refFile = refFile, targetFile = targetFile){
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase="geneontology_Biological_Process_noRedundant", interestGeneFile=targetFile,
  interestGeneType="genesymbol", referenceGeneFile=refFile,
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=enrichmentDir, projectName=paste0(projectName,'_BP'), fdrThr=1)

enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase="geneontology_Cellular_Component_noRedundant", interestGeneFile=targetFile,
  interestGeneType="genesymbol", referenceGeneFile=refFile,
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=enrichmentDir, projectName=paste0(projectName,'_CC'), fdrThr=1)
  
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase="geneontology_Molecular_Function_noRedundant", interestGeneFile=targetFile,
  interestGeneType="genesymbol", referenceGeneFile=refFile,
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=enrichmentDir, projectName=paste0(projectName,'_MF'), fdrThr=1)

enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase="pathway_KEGG", interestGeneFile=targetFile,
  interestGeneType="genesymbol", referenceGeneFile=refFile,
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=enrichmentDir, projectName=paste0(projectName,'_KEGG'), fdrThr=1)
}

maxPower = function(vector){
	vector = vector[vector != 0]
	return(ceiling(-log10(vector[1])))
}

enrichmentTop5Df = function(projectName = projectName, enrichmentDir = enrichmentDir){
bp_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_BP/enrichment_results_',projectName,'_BP.txt'), header = T, sep = '\t')
cc_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_CC/enrichment_results_',projectName,'_CC.txt'), header = T, sep = '\t')
mf_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_MF/enrichment_results_',projectName,'_MF.txt'), header = T, sep = '\t')
kegg_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_KEGG/enrichment_results_',projectName,'_KEGG.txt'), header = T, sep = '\t')

bp_df$FDR = gsub('^0$',paste0('1e-',maxPower(bp_df$FDR)),bp_df$FDR)
bp_df$FDR = as.numeric(bp_df$FDR)
cc_df$FDR = gsub('^0$',paste0('1e-',maxPower(cc_df$FDR)),cc_df$FDR)
cc_df$FDR = as.numeric(cc_df$FDR)
mf_df$FDR = gsub('^0$',paste0('1e-',maxPower(mf_df$FDR)),mf_df$FDR)
mf_df$FDR = as.numeric(mf_df$FDR)
kegg_df$FDR = gsub('^0$',paste0('1e-',maxPower(kegg_df$FDR)),kegg_df$FDR)
kegg_df$FDR = as.numeric(kegg_df$FDR)
bp_df$log10fdr = -log10(bp_df$FDR)
cc_df$log10fdr = -log10(cc_df$FDR)
mf_df$log10fdr = -log10(mf_df$FDR)
kegg_df$log10fdr = -log10(kegg_df$FDR)

#add category
bp_df$cat = 'BP'
cc_df$cat = 'CC'
mf_df$cat = 'MF'
kegg_df$cat = 'KEGG'

#combine the top five terms from each category
top5_df = rbind(bp_df[1:5,],cc_df[1:5,],mf_df[1:5,],kegg_df[1:5,])
top5_df$description = factor(top5_df$description, levels = rev(top5_df$description))
return(top5_df)
}

enrichmentTop5Df_GOonly = function(projectName = projectName, enrichmentDir = enrichmentDir){
bp_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_BP/enrichment_results_',projectName,'_BP.txt'), header = T, sep = '\t')
cc_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_CC/enrichment_results_',projectName,'_CC.txt'), header = T, sep = '\t')
mf_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_MF/enrichment_results_',projectName,'_MF.txt'), header = T, sep = '\t')

bp_df$FDR = gsub('^0$',paste0('1e-',maxPower(bp_df$FDR)),bp_df$FDR)
bp_df$FDR = as.numeric(bp_df$FDR)
cc_df$FDR = gsub('^0$',paste0('1e-',maxPower(cc_df$FDR)),cc_df$FDR)
cc_df$FDR = as.numeric(cc_df$FDR)
mf_df$FDR = gsub('^0$',paste0('1e-',maxPower(mf_df$FDR)),mf_df$FDR)
mf_df$FDR = as.numeric(mf_df$FDR)
bp_df$log10fdr = -log10(bp_df$FDR)
cc_df$log10fdr = -log10(cc_df$FDR)
mf_df$log10fdr = -log10(mf_df$FDR)

#add category
bp_df$cat = 'BP'
cc_df$cat = 'CC'
mf_df$cat = 'MF'

#combine the top five terms from each category
top5_df = rbind(bp_df[1:5,],cc_df[1:5,],mf_df[1:5,])
top5_df$description = factor(top5_df$description, levels = rev(top5_df$description))
return(top5_df)
}

#pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
#dot_df3 %>% filter(cell_exp > 0, cell_exp_pct > 1) %>%
enrichmentDotPlot = function(top5_df) {
p = ggplot(top5_df, aes(x=enrichmentRatio, y = description, color = log10fdr, size = size)) +
  geom_point() +
  scale_colour_gradient2(name = '-log10(FDR)', low = 'blue', mid = 'grey', high = 'red', midpoint = midpoint) +
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 8, angle = 0, vjust = 1, hjust=1), axis.title=element_text(size=8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  ylab('') +
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = c(5.5,10.5,15.5),linetype = 'dashed')
#dev.off()
return(p)
}

enrichmentDotPlot_GOonly = function(top5_df) {
p = ggplot(top5_df, aes(x=enrichmentRatio, y = description, color = log10fdr, size = size)) +
  geom_point() +
  scale_colour_gradient2(name = '-log10(FDR)', low = 'blue', mid = 'grey', high = 'red', midpoint = midpoint) +
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 8, angle = 0, vjust = 1, hjust=1), axis.title=element_text(size=8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  ylab('') +
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = c(5.5,10.5),linetype = 'dashed')
#dev.off()
return(p)
}

###enrichment
#library(WebGestaltR)
library(cowplot)
enrichmentDir = paste0(work_dir,'enrichmentDir/')
if(!dir.exists(enrichmentDir)){dir.create(enrichmentDir,recursive = T)}
refFile = 'C:/Users/shann/Documents/files4enrichment/protein-coding_gene_v2.txt'

#GM and WM shared genes
targetGenes = intersect(gm_genes, wm_genes)
write.table(targetGenes,file = paste0(work_dir,'target_gene.txt'), quote = F, row.names = F, col.names = F)
targetFile = paste0(work_dir,'target_gene.txt')
projectName = 'GM_and_WM_shared_genes'
WebGestaltR_GO_KEGG(projectName,enrichmentDir,refFile,targetFile)
top5_df = enrichmentTop5Df(projectName,enrichmentDir)
midpoint = mean(top5_df$log10fdr)
p = enrichmentDotPlot(top5_df)
pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
print(p)
dev.off()
p

#GM specific genes
targetGenes = gm_genes[!gm_genes %in% wm_genes]
write.table(targetGenes,file = paste0(work_dir,'target_gene.txt'), quote = F, row.names = F, col.names = F)
targetFile = paste0(work_dir,'target_gene.txt')
projectName = 'GM_specific_genes'
WebGestaltR_GO_KEGG(projectName,enrichmentDir,refFile,targetFile)
top5_df = enrichmentTop5Df(projectName,enrichmentDir)
midpoint = mean(top5_df$log10fdr)
p = enrichmentDotPlot(top5_df)
pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
print(p)
dev.off()
p

#WM specific genes
targetGenes = wm_genes[!wm_genes %in% gm_genes]
write.table(targetGenes,file = paste0(work_dir,'target_gene.txt'), quote = F, row.names = F, col.names = F)
targetFile = paste0(work_dir,'target_gene.txt')
projectName = 'WM_specific_genes'
WebGestaltR_GO_KEGG(projectName,enrichmentDir,refFile,targetFile)
top5_df = enrichmentTop5Df_GOonly(projectName,enrichmentDir)
midpoint = mean(top5_df$log10fdr)
p = enrichmentDotPlot_GOonly(top5_df)
pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
print(p)
dev.off()
p

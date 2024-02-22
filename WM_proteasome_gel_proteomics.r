# Loading the package for proteasome analysis
library(DEP)
# Loading a package required for data handling
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#helper function
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

work_dir = 'C:/Users/shann/Documents/Natura/proteasome_project/proteomics/gel_proteomics/files/'
proteasome_gene_df = read.table(file = paste0('C:/Users/shann/Documents/Natura/', 'gene_id.txt'), header = T, sep = '\t')
prot_dat = read.csv(file = paste0(work_dir,'WM_proteinGroups.csv'), header = T)
prot_dat[is.na(prot_dat)] = 0
colnames(prot_dat) = gsub('LFQ.intensity.Control','LFQ.intensity.CTR',colnames(prot_dat))
colnames(prot_dat) = gsub('LFQ.intensity.Ctr','LFQ.intensity.CTR',colnames(prot_dat))
colnames(prot_dat) = gsub('LFQ.intensity.CTR_','LFQ.intensity.CTR',colnames(prot_dat))
colnames(prot_dat) = gsub('LFQ.intensity.AD_','LFQ.intensity.AD',colnames(prot_dat))

# Are there any duplicated gene names?
prot_dat$Gene.names %>% duplicated() %>% any()

# Make a table of duplicated gene names
prot_dat %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
prot_unique <- make_unique(prot_dat, "Gene.names", "Protein.IDs", delim = ";")

# Are there any duplicated names?
prot_unique$name %>% duplicated() %>% any()

# Generate a SummarizedExperiment object using an experimental design
prot_LFQ_columns <- grep("LFQ.", colnames(prot_unique)) # get LFQ column numbers
#experimental_design <- UbiLength_ExpDesign
experimental_design = data.frame(label = colnames(prot_dat)[prot_LFQ_columns])
experimental_design$label = gsub('LFQ.intensity.','',experimental_design$label)
experimental_design$condition = rep(c('AD','CTR'),c(8,6))
experimental_design$replicate = c(c(1:8),c(1:6))
#experimental_design = data.frame(label = c('AD_5650',), condition, replicate)
data_se <- make_se(prot_unique, prot_LFQ_columns, experimental_design)

data_se_parsed <- make_se_parse(prot_unique, prot_LFQ_columns)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

##############filter out proteins of low quality
# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

################imputation
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation data_imp
plot_imputation(data_norm, data_imp)

# Plot intensity distributions before and after imputation data_imp_man
plot_imputation(data_norm, data_imp_man)

# Plot intensity distributions before and after imputation data_imp_knn
plot_imputation(data_norm, data_imp_knn)
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "CTR")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1.5))

#################visualization of the results
# Plot the first and second principal components (PCA)
plot_pca(dep, x = 1, y = 2, indicate = c("condition"), n = 50, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 2, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a volcano plot for the contrast "AD vs CTR""
plot_volcano(dep, contrast = "AD_vs_CTR", label_size = 2, add_names = TRUE)
# output result
data_imp_res = get_results(dep)
write.csv(data_imp_res, file = paste0(work_dir,'WM_protein_data_imp_res.csv'))

############## all gene heatmap
datExpr = as.data.frame(data_imp@assays@data)
datExpr = datExpr[,!colnames(datExpr) %in% c('group','group_name')]
col_annot_df = data.frame(condition = rep(c('CTR','AD'),c(6,8)))
rownames(col_annot_df) = c(paste0('CTR_',c(1:6)),paste0('AD_',c(1:8)))
datExpr = datExpr[,rownames(col_annot_df)]

fig = pheatmap(datExpr, scale = 'row', cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, fontsize_row = 6, annotation_col = col_annot_df, annotation_row = NA)
save_pheatmap_pdf(fig, paste0(work_dir,'WM_gel_proteomics_proteasome_gene_heatmap.pdf'), width = 3.7, height = 7)

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

sharedGenes = intersect(rownames(datExpr),rownames(gene_substructure_df))
row_annot_df = gene_substructure_df[rownames(gene_substructure_df) %in% sharedGenes,]
row_annot_df = row_annot_df[c('substructure')]
row_annot_df$substructure = factor(row_annot_df$substructure, levels = substrLevel)
datExpr2 = datExpr[rownames(datExpr) %in% sharedGenes,]
datExpr2 = datExpr2[rownames(row_annot_df),]

annot_colors = list(condition = c(AD = "#1B9E77", CTR = "#D95F02"),substructure = c(E_ligase = '#8DD3C7',DUBs_deubiquitylating_enzyme = '#FFFFB3',HSP_protein_chaperon = '#BEBADA',P19S = '#FB8072',P20S = '#80B1D3',IP = '#FDB462',Assembly_chaperon = '#B3DE69',Activator = '#FCCDE5'))
fig = pheatmap(datExpr2, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdBu")))(100), cellheight = 12, cellwidth = 12, scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = F, fontsize = 14, annotation_colors = annot_colors, annotation_col = col_annot_df, annotation_row = row_annot_df)
save_pheatmap_pdf(fig, paste0(work_dir,'gel_proteomics_WM_proteasome_gene_heatmap.pdf'), width = 8, height = 10)


############## Plot a volcano plot with genes passing nominal P = 0.05 highlighted
library(ggplot2)
#library(ggrepel)
# read in data_imp_result file and proteasome gene id list
data_imp_res$log10P = -log10(data_imp_res$AD_vs_CTR_p.val)
data_imp_res$volcano_color = with(data_imp_res, ifelse(AD_vs_CTR_p.val < 0.05, 'black', 'grey'))
volcano_colr = as.character(unique(data_imp_res$volcano_color))
#mark a new threshold passing -log10(P) > 3
data_imp_res$annot_thre = with(data_imp_res,ifelse(-log10(AD_vs_CTR_p.val) > -log10(0.05), 'TRUE', 'FALSE'))
#options(ggrepel.max.overlaps = Inf)
pdf(paste0(work_dir,'gel_proteomics_WM_all_gene_volcano.pdf'),width = 3.3, height = 3.3)
ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P)) + geom_point(aes(color=volcano_color),size=0.5,alpha=0.5) + scale_color_manual(breaks = unique(data_imp_res$volcano_color), values = volcano_colr) + theme(axis.text = element_text(size = 12),axis.title=element_text(size=12),legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_text(data=subset(data_imp_res,annot_thre %in% c('TRUE')),aes(label=name),size=3,nudge_x = 0,nudge_y = 0.2,check_overlap = TRUE)
dev.off()


############## Plot a volcano plot with proteasome genes
#library(ggplot2)
#library(ggrepel)
# read in data_imp_result file and proteasome gene id list
#data_imp_res = read.csv(file = paste0(work_dir, 'protein_data_imp_result.csv'), header = T, row.names = 1)
data_imp_res$log10P = -log10(data_imp_res$AD_vs_CTR_p.val)
data_imp_res$volcano_color = with(data_imp_res, ifelse(name %in% proteasome_gene_df$OGS & AD_vs_CTR_ratio < 0, 'blue', ifelse(name %in% proteasome_gene_df$OGS & AD_vs_CTR_ratio > 0, 'red', 'grey')))
data_imp_res$volcano_color = factor(data_imp_res$volcano_color, levels = c('grey', 'blue', 'red'))
volcano_colr = as.character(unique(data_imp_res$volcano_color))
#options(ggrepel.max.overlaps = 30)
pdf(paste0(work_dir,'gel_proteomics_WM_proteasome_gene_volcano.pdf'),width = 3.3, height = 3.3)
ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P)) + geom_point(aes(color=volcano_color),size=0.5,alpha=0.5) + geom_point(data=subset(data_imp_res,volcano_color %in% c('blue','red')),aes(color=volcano_color),size=0.5,alpha=0.5) + scale_color_manual(breaks = unique(data_imp_res$volcano_color), values = volcano_colr) + theme(axis.text = element_text(size = 12),axis.title=element_text(size=12),legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_vline(xintercept=-0.5,linetype='dashed',color='black') + geom_vline(xintercept=0.5,linetype='dashed',color='black') + geom_text(data=subset(data_imp_res,volcano_color %in% c('blue','red')),aes(label=name),size=3,nudge_x = 0,nudge_y = 0.2,check_overlap = TRUE)
dev.off()

############## Plot a volcano plot with all proteasome genes
#data_imp_res = read.csv(file = paste0(work_dir, 'protein_data_imp_result.csv'), header = T, row.names = 1)
#proteasome_gene_df = read.table(file = paste0(work_dir, 'gene_id.txt'), header = T, sep = '\t')
data_imp_res$log10P = -log10(data_imp_res$AD_vs_CTR_p.val)
data_imp_res$volcano_color = with(data_imp_res, ifelse(name %in% proteasome_gene_df$OGS & AD_vs_CTR_ratio < 0, 'blue', ifelse(name %in% proteasome_gene_df$OGS & AD_vs_CTR_ratio > 0, 'red', 'grey')))
data_imp_res$volcano_color = factor(data_imp_res$volcano_color, levels = c('grey', 'blue', 'red'))
volcano_colr = as.character(unique(data_imp_res$volcano_color))
options(ggrepel.max.overlaps = Inf)
pdf(paste0(work_dir,'gel_proteomics_WM_proteasome_gene_volcano.pdf'),width = 7, height = 7)
ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P)) + geom_point(aes(color=volcano_color),size=0.75,alpha=0.5) + geom_point(data=subset(data_imp_res,volcano_color %in% c('blue','red')),aes(color=volcano_color),size=0.75,alpha=1) + scale_color_manual(breaks = unique(data_imp_res$volcano_color), values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_vline(xintercept=-0.5,linetype='dashed',color='black') + geom_vline(xintercept=0.5,linetype='dashed',color='black') + geom_text_repel(data=subset(data_imp_res,volcano_color %in% c('blue','red')),aes(label=name),size=2)
dev.off()

############# Plot a volcano plot proteasome subunit color-coded
library(plyr)

#data_imp_res = read.csv(file = paste0(work_dir, 'protein_data_imp_result.csv'), header = T, row.names = 1)
#proteasome_gene_df = read.table(file = paste0(work_dir, 'gene_id.txt'), header = T, sep = '\t')
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
rownames(proteasome_gene_df) = proteasome_gene_df$gene
proteasome_gene_df$gene = proteasome_gene_df$OGS
gene_substructure_df = join(gene_substructure_df,proteasome_gene_df,by = 'gene')
rownames(gene_substructure_df) = gene_substructure_df$ENSG_ID

data_imp_res$substructure = with(data_imp_res, ifelse(name %in% proteasome_19s_genes, '19S', ifelse(name %in% proteasome_20s_genes, '20S', ifelse(name %in% immunoproteasome_genes, 'Immunoproteasome', ifelse(name %in% proteasome_chaperon_genes, 'Chaperon', ifelse(name %in% proteasome_activator_genes, 'Activator', ifelse(name %in% proteasome_tf_genes, 'TF', ifelse(name %in% proteasome_receptor_genes, 'Receptor', 'non-proteasome'))))))))

data_imp_res$log10P = -log10(data_imp_res$AD_vs_CTR_p.val)
data_imp_res$volcano_color = with(data_imp_res, ifelse(substructure == 'non-proteasome', 'grey', ifelse(substructure == '19S', 'red', ifelse(substructure == '20S', 'blue', ifelse(substructure == 'Activator', 'yellow', ifelse(substructure == 'Chaperon', 'green', ifelse(substructure == 'Immunoproteasome', 'purple', 'brown')))))))
data_imp_res$volcano_color = factor(data_imp_res$volcano_color, levels = c('grey', 'red', 'blue', 'yellow', 'green', 'purple', 'brown'))
volcano_colr = as.character(unique(data_imp_res$volcano_color))
substructure = as.character(unique(data_imp_res$substructure))
data_imp_res$shape = 16

options(ggrepel.max.overlaps = Inf)
ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P)) + geom_point(aes(color=volcano_color),size=0.75,alpha=0.5) + geom_point(data=subset(data_imp_res,volcano_color %in% c('red', 'blue', 'yellow', 'green', 'purple', 'brown')),aes(color=volcano_color),size=0.75,alpha=1) + scale_color_manual(breaks = unique(data_imp_res$volcano_color), values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_vline(xintercept=-0.5,linetype='dashed',color='black') + geom_vline(xintercept=0.5,linetype='dashed',color='black') + geom_text_repel(data=subset(data_imp_res,volcano_color %in% c('red', 'blue', 'yellow', 'green', 'purple', 'brown')),aes(label=name),size=2)

ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P,color=volcano_color)) + geom_point(aes(color=volcano_color),size=0.75,alpha=0.5) + geom_point(data=subset(data_imp_res,volcano_color %in% c('red', 'blue', 'yellow', 'green', 'purple', 'brown')),aes(color=volcano_color),size=0.75,alpha=1) + scale_color_manual(name = 'substructure', breaks = unique(data_imp_res$volcano_color), labels = substructure, values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_vline(xintercept=-0.5,linetype='dashed',color='black') + geom_vline(xintercept=0.5,linetype='dashed',color='black') + geom_text_repel(data=subset(data_imp_res,volcano_color %in% c('red', 'blue', 'yellow', 'green', 'purple', 'brown')),aes(label=name),size=2) + theme(legend.position = "bottom")

############# Plot a volcano plot proteasome subunit color-coded without receptor
library(plyr)

#data_imp_res = read.csv(file = paste0(work_dir, 'protein_data_imp_result.csv'), header = T, row.names = 1)
#proteasome_gene_df = read.table(file = paste0(work_dir, 'gene_id.txt'), header = T, sep = '\t')
#proteasome_gene_annotation
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1','MAFF','MAFG','MAFK')
proteasome_tf_genes = c('NFE2L1')
#proteasome_receptor_genes = c('ADRM1')

substructure = rep(c('19S','20S','Immunoproteasome','Chaperon','Activator','TF'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes)))
gene = c(proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes,proteasome_tf_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene
rownames(proteasome_gene_df) = proteasome_gene_df$gene
proteasome_gene_df$gene = proteasome_gene_df$OGS
gene_substructure_df = join(gene_substructure_df,proteasome_gene_df,by = 'gene')
rownames(gene_substructure_df) = gene_substructure_df$ENSG_ID

data_imp_res$substructure = with(data_imp_res, ifelse(name %in% proteasome_19s_genes, '19S', ifelse(name %in% proteasome_20s_genes, '20S', ifelse(name %in% immunoproteasome_genes, 'Immunoproteasome', ifelse(name %in% proteasome_chaperon_genes, 'Chaperon', ifelse(name %in% proteasome_activator_genes, 'Activator', ifelse(name %in% proteasome_tf_genes, 'TF', 'non-proteasome')))))))

data_imp_res$log10P = -log10(data_imp_res$AD_vs_CTR_p.val)
data_imp_res$volcano_color = with(data_imp_res, ifelse(substructure == 'non-proteasome', 'grey', ifelse(substructure == '19S', 'red', ifelse(substructure == '20S', 'blue', ifelse(substructure == 'Activator', 'brown', ifelse(substructure == 'Chaperon', 'green', 'purple'))))))
data_imp_res$volcano_color = factor(data_imp_res$volcano_color, levels = c('grey', 'red', 'blue', 'brown', 'green', 'purple'))
volcano_colr = as.character(unique(data_imp_res$volcano_color))
substructure = as.character(unique(data_imp_res$substructure))
data_imp_res$shape = 16

options(ggrepel.max.overlaps = Inf)
ggplot(data_imp_res,aes(x=AD_vs_CTR_ratio,y=log10P,color=volcano_color)) + geom_point(aes(color=volcano_color),size=0.75,alpha=0.5) + geom_point(data=subset(data_imp_res,volcano_color %in% c('red', 'blue', 'green', 'purple', 'brown')),aes(color=volcano_color),size=0.75,alpha=1) + scale_color_manual(name = 'substructure', breaks = unique(data_imp_res$volcano_color), labels = substructure, values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + geom_vline(xintercept=-0.5,linetype='dashed',color='black') + geom_vline(xintercept=0.5,linetype='dashed',color='black') + geom_text_repel(data=subset(data_imp_res,name %in% c('PSMG4')),aes(label=name),size=2) + theme(legend.position = "bottom")


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

#pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
#dot_df3 %>% filter(cell_exp > 0, cell_exp_pct > 1) %>%
enrichmentDotPlot = function(top5_df) {
p = ggplot(top5_df, aes(x=enrichmentRatio, y = description, color = log10fdr, size = size)) +
  geom_point() +
  scale_colour_gradient2(name = '-log10(FDR)', low = 'blue', mid = 'grey', high = 'red', midpoint = midpoint) +
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 12, angle = 0, vjust = 1, hjust=1), axis.title=element_text(size=12), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  ylab('') +
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = c(5.5,10.5,15.5),linetype = 'dashed')
#dev.off()
return(p)
}

###enrichment
library(WebGestaltR)
library(cowplot)
enrichmentDir = paste0(work_dir,'enrichmentDir/')
if(!dir.exists(enrichmentDir)){dir.create(enrichmentDir,recursive = T)}
refFile = 'C:/Users/shann/Documents/files4enrichment/protein-coding_gene_v2.txt'

#all WM genes
targetGenes = data_imp_res$name
write.table(targetGenes,file = paste0(work_dir,'target_gene.txt'), quote = F, row.names = F, col.names = F)
targetFile = paste0(work_dir,'target_gene.txt')
projectName = 'all_WM_genes'
WebGestaltR_GO_KEGG(projectName,enrichmentDir,refFile,targetFile)
top5_df = enrichmentTop5Df(projectName,enrichmentDir)
midpoint = mean(top5_df$log10fdr)
p = enrichmentDotPlot(top5_df)
pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
print(p)
dev.off()
p


#WM nominal sig genes
targetGenes = data_imp_res[data_imp_res$AD_vs_CTR_p.val<0.05,]$name
write.table(targetGenes,file = paste0(work_dir,'target_gene.txt'), quote = F, row.names = F, col.names = F)
targetFile = paste0(work_dir,'target_gene.txt')
projectName = 'nominal_sig_WM_genes'
WebGestaltR_GO_KEGG(projectName,enrichmentDir,refFile,targetFile)
top5_df = enrichmentTop5Df(projectName,enrichmentDir)
midpoint = mean(top5_df$log10fdr)
p = enrichmentDotPlot(top5_df)
pdf(file = paste0(enrichmentDir, projectName,'_top5_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top5_df$description))
print(p)
dev.off()
p

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
#print(ad.__version__)

#set environment variable
work_dir = '/home/sj3195/proteasome/data/snrnaseq/SEAAD/'

#read in h5ad file
#L2n3_IT = sc.read_h5ad(filename = work_dir + 'L2n3_IT.h5ad')
celltypes = ['L2n3_IT','Astro','Chandelier','Endo','L4_IT','L5_ET','L5_IT','L5n6_NP','L6b','L6_CT','L6_IT_Car3','L6_IT','Lamp5','Lamp5_Lhx6','Micro-PVM','Oligo','OPC','Pax6','Pvalb','Sncg','Sst_Chodl','Sst','Vip','VLMC']
for celltype in celltypes:
    adata = sc.read_h5ad(filename = work_dir + celltype + '.h5ad')
    #adata.var['ENSG_ID'] = adata.var_names

    ####keep only proteasome genes
    gene_id_df = pd.read_csv(work_dir + 'htr4_proj_gene_id.txt', sep = '\t',header = 0)
    #remove ADRM1 MAFF MAFG and MAFK
    #gene_id_df = gene_id_df[gene_id_df['OGS'].str.contains('ADRM1|MAFF|MAFG|MAFK')==False]
    #prot_adata = adata[:,gene_id_df['ENSG_ID'].to_list()]
    ####remove reference samples
    adata = adata[~adata.obs.ADNC.isin(['Reference'])]
    meta_df = adata.obs
    adata_var_df = adata.var
    cts = pd.DataFrame(adata.raw.X.todense())
    cts.index = meta_df.index
    cts.columns = adata_var_df.index
    #screen prot cts
    prot_cts = cts[gene_id_df['ENSG_ID'].to_list()]
    prot_cts[gene_id_df['ENSG_ID'].to_list()] = prot_cts[gene_id_df['ENSG_ID'].to_list()].astype(int)
    prot_cts.to_csv(work_dir + celltype + '_htr4_proj_gene_cts.csv')
    #meta_df.to_csv(work_dir + celltype + '_meta.csv')

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix


# In[9]:


work_dir = '/home/shann/Documents/cziscience/SEAAD/'
out_dir = work_dir
gene_id_df = pd.read_csv('/home/shann/Documents/Natura/proteosome/gene_id.txt', sep = '\t', header = 0)
prefix = 'prot'


# In[3]:


celltypes = ['L23_IT','Astrocyte','Chandelier','Endothelial','L4_IT','L5_ET','L5_IT','L56_NP','L6b','L6_CT','L6_IT_Car3','L6_IT','Lamp5','Lamp5_Lhx6','Microglia_PVM','Oligo','OPC','Pax6','Pvalb','Sncg','Sst_Chodl','Sst','Vip','VLMC']
brainRegions = ['MTG','DLPFC']


# In[4]:


celltype = 'Endothelial'
brainRegion = 'MTG'


# In[5]:


for brainRegion in brainRegions:
    for celltype in celltypes:
        adata = sc.read_h5ad(filename = work_dir + celltype + '_' + brainRegion + '.h5ad')
        meta_df = adata.obs
        adata_var_df = adata.var
        #normalized expr
        normDat = pd.DataFrame(adata.X.todense())
        normDat.index = meta_df.index
        normDat.columns = adata_var_df.index
        targetGeneNormDat = normDat[gene_id_df['ENSG_ID'].to_list()]
        targetGeneNormDat.to_csv(out_dir + prefix + '_' + celltype + '_' + brainRegion + '_normDat.csv')


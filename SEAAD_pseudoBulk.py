#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import pandas as pd
#pd.options.display.max_columns = None
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix


# In[43]:


#set environment variable
work_dir = '/home/shann/Documents/cziscience/SEAAD/pseudoBulk/'
out_dir = work_dir


# In[3]:


#celltype = 'Sncg'
#brainRegion = 'MTG'
celltypes = ['L23_IT','Astrocyte','Chandelier','Endothelial','L4_IT','L5_ET','L5_IT','L56_NP','L6b','L6_CT','L6_IT_Car3','L6_IT','Lamp5','Lamp5_Lhx6','Microglia_PVM','Oligo','OPC','Pax6','Pvalb','Sncg','Sst_Chodl','Sst','Vip','VLMC']
brainRegions = ['MTG','DLPFC']

for brainRegion in brainRegions:
    for celltype in celltypes:
        adata = sc.read_h5ad(filename = work_dir + celltype + '_' + brainRegion + '.h5ad')


        # In[4]:


        meta_df = adata.obs
        adata_var_df = adata.var
        cts = pd.DataFrame(adata.raw.X.todense())


        # In[31]:


        cts.index = meta_df.index
        cts.columns = adata_var_df.index


        # In[37]:


        ENSG_ID = cts.columns.tolist()
        cts['donor_id'] = meta_df['donor_id']


        # In[54]:


        indAggCts = cts.groupby('donor_id')[ENSG_ID].sum()


        # In[67]:


        for c in indAggCts.columns:
            if indAggCts[c].dtype == np.float32:
                indAggCts[c] = indAggCts[c].astype(int)


        # In[69]:


        meta_df.to_csv(out_dir + celltype + '_' + brainRegion + '_meta.csv')
        adata_var_df.to_csv(out_dir + celltype + '_' + brainRegion + '_gene_annot.csv')
        indAggCts.to_csv(out_dir + celltype + '_' + brainRegion + '_indAggCts.csv')

#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix


def gene_len_adj(cts_list):
    adj_cts_list = [m / n for m, n in zip(cts_list, CDS_length)]
    return adj_cts_list


# set environment variable
work_dir = "/home/sj3195/cziscience/SEAAD/"
exon_len_df = pd.read_csv(
    "/home/sj3195/references/ensembl/ensembl.GRCh38.108.exon_lengths_nodupOGS.txt",
    sep="\t",
    header=0,
    index_col="ENSG_ID"
)
# celltype = "L6b"
celltypes = [
    "L2n3_IT",
    "Astro",
    "Chandelier",
    "Endo",
    "L4_IT",
    "L5_ET",
    "L5_IT",
    "L5n6_NP",
    "L6b",
    "L6_CT",
    "L6_IT_Car3",
    "L6_IT",
    "Lamp5",
    "Lamp5_Lhx6",
    "Micro-PVM",
    "Oligo",
    "OPC",
    "Pax6",
    "Pvalb",
    "Sncg",
    "Sst_Chodl",
    "Sst",
    "Vip",
    "VLMC",
]
for celltype in celltypes:
    adata = sc.read_h5ad(filename=work_dir + celltype + ".h5ad")
    meta_df = adata.obs
    adata_var_df = adata.var
    cts = pd.DataFrame(adata.raw.X.todense())
    cts.index = meta_df.index
    cts.columns = adata_var_df.index
    cts_T = cts.T
    share_ENSG_list = cts_T.index.intersection(exon_len_df.index).tolist()
    share_cts = cts_T[cts_T.index.isin(share_ENSG_list)]
    share_exon_len_df = exon_len_df[exon_len_df.index.isin(share_ENSG_list)]
    share_exon_len_df = share_exon_len_df.reindex(share_cts.index)
    CDS_length = share_exon_len_df["CDS_length"].tolist()
    adj_cts = share_cts.divide(CDS_length, axis=0)
    adj_cts_1e4 = adj_cts * 1e4
    tp10k = adj_cts_1e4.divide(adj_cts.sum(axis=0), axis=1)
    tp10k_spa_mtx = csr_matrix(tp10k.T.values)
    tp10k_adata = ad.AnnData(tp10k_spa_mtx, dtype="float32")
    tp10k_adata.obs_names = tp10k.columns.tolist()
    tp10k_adata.var_names = tp10k.index.tolist()
    tp10k_adata.write(work_dir + celltype + "_tp10k.h5ad", compression="gzip")

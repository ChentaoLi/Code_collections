import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import scanpy as sc
import os
#载入 loom文件
N1= anndata.read_loom("./N1/velocyto/N2.loom")
N2= anndata.read_loom("./N2/velocyto/N2.loom")
N3= anndata.read_loom("./N3/velocyto/N3.loom")
R1= anndata.read_loom("./R1/velocyto/R1.loom")
R2= anndata.read_loom("./R2/velocyto/R2.loom")
R3= anndata.read_loom("./R3/velocyto/R3.loom")
#修改列名
barcodes = [bc.split(':')[1] for bc in N1.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['N1_'+ bc[0:len(bc)] for bc in barcodes]
N1.obs.index = barcodes
N1.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in N2.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['N2_'+ bc[0:len(bc)] for bc in barcodes]
N2.obs.index = barcodes
N2.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in N3.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['N3_'+ bc[0:len(bc)] for bc in barcodes]
N3.obs.index = barcodes
N3.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in R1.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['R1_'+ bc[0:len(bc)] for bc in barcodes]
R1.obs.index = barcodes
R1.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in R2.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['R2_'+ bc[0:len(bc)] for bc in barcodes]
R2.obs.index = barcodes
R2.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in R3.obs.index.tolist()] 
barcodes=[bc[0:len(bc)-1]+'-1' for bc in barcodes]
barcodes=['R3_'+ bc[0:len(bc)] for bc in barcodes]
R3.obs.index = barcodes
R3.var_names_make_unique()
ldata= N1.concatenate([N2,N3,R1,R2,R3])
adata=sc.read_h5ad('./data/project1/scRNA.T.h5ad')
adata= scv.utils.merge(adata,ldata)
scv.pp.filter_and_normalize(adata,min_shared_counts=20,n_top_genes=2000)
scv.pp.moments(adata,n_pcs=30,n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata,mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata,basis='umap',color='celltype', title='', smooth=.8, min_mass=4)
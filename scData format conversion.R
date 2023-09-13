# Transform single-cell object
# 1. SeuratDisk ha5d to Seurat
library(Seurat)
library(SeuratData)
library(SeuratDisk)
# load data
url <- "https://seurat.nygenome.org/pbmc3k_final.h5ad"
curl::curl_download(url, basename(url))
Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")

# 2. Transform files with sceasy
# https://github.com/cellgeni/sceasy
library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
loompy <- reticulate::import('loompy')

sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                      outFile='filename.h5ad')

sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='filename.rds')

sceasy::convertFormat(seurat_object, from="seurat", to="sce",
                      outFile='filename.rds')

sceasy::convertFormat(sce_object, from="sce", to="anndata",
                      outFile='filename.h5ad')

sceasy::convertFormat(sce_object, from="sce", to="loom",
                      outFile='filename.loom')

sceasy::convertFormat('filename.loom', from="loom", to="anndata",
                      outFile='filename.h5ad')

sceasy::convertFormat('filename.loom', from="loom", to="sce",
                      outFile='filename.rds')

# 3. Transform files with scDIOR
# https://github.com/JiekaiLab/scDIOR
# devtools::install_github('JiekaiLab/dior')
# for python; pip install diopy

# Single-cell data from R to Python
dior::write_h5(data, file='scdata.h5' object.type = 'singlecellexperiment')
# in Python; adata = diopy.input.read_h5(file = 'scdata.h5')

# Single-cell data from Python to R
# in Python; diopy.output.write_h5(data_py, file = 'scdata.h5')
adata = dior::read_h5(file='scdata.h5', target.object = 'seurat')

# spatial omics data
dior::write_h5(data, file='scdata.h5', object.type = 'singlecellexperiment')
# in Python; adata = diopy.input.read_h5(file = 'scdata.h5')


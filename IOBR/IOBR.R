# Immuno-Oncology Biological Research
# https://iobr.github.io/IOBR/IOBR-VIGNETTE.html
###########################################
# 1. 安装
## 安装依赖包
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
          "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
          "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr')

for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}
# BiocManager::install("clusterProfiler")
# BiocManager::install("maftools")
options(timeout=9999999)
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR",force=TRUE)

###########################################
# 2. 数据准备
library(dplyr)
library(IOBR)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
# BiocManager::install("UCSCXenaTools")
library(UCSCXenaTools) # for TCGA

eset_stad<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Stomach Cancer (STAD)") %>% 
  XenaFilter(filterDatasets    = "TCGA-STAD.htseq_counts.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()
#TCGA数据清洗
# Remove the version numbers in Ensembl ID.
eset_stad$Ensembl_ID<-substring(eset_stad$Ensembl_ID, 1, 15)
eset_stad<-column_to_rownames(eset_stad, var = "Ensembl_ID")
# Revert back to original format because the data from UCSC was log2(x+1)transformed.
eset_stad<-(2^eset_stad)+1
# In this process, *biomaRt* R package is utilized to acquire the gene length of each Ensembl ID and calculate the TPM of each sample. If identical gene symbols exists, these genes would be ordered by the mean expression levels. The gene symbol with highest mean expression level is selected and remove others.
eset_stad<-count2tpm(countMat = eset_stad, idType = "Ensembl", org="hsa", source = "local" )

#鉴定超出值
res <- find_outlier_samples(eset = eset_stad, project = "ACRG", show_plot = TRUE)
#剔除超出样本
eset1 <- eset_stad[, !colnames(eset_stad)%in%res]

library("GEOquery")
eset_geo<-getGEO(GEO     = "GSE100935",
                 getGPL  = F,
                 destdir = "./")
eset    <-eset_geo[[1]]
eset    <-exprs(eset)
eset<-anno_eset(eset       = eset,
                annotation = anno_hug133plus2,
                symbol     = "symbol",
                probe      = "probe_id",
                method     = "mean")


####----PCA分析----####
data("pdata_acrg")
res<- iobr_pca(data       = eset1,
               is.matrix   = TRUE,
               scale       = TRUE,
               is.log      = FALSE,
               pdata       = pdata_acrg, 
               id_pdata    = "ID", 
               group       = "Subtype",
               geom.ind    = "point", 
               cols        = "normal",
               palette     = "jama", 
               repel       = FALSE,
               ncp         = 5,
               axes        = c(1, 2),
               addEllipses = TRUE)
res + geom_point()

####----批次校正----####
# NOTE: This process may take a few minutes which depends on the internet connection speed. Please wait for its completion.
eset_geo <- getGEO(GEO     = "GSE57303", getGPL  = F, destdir = "./")
eset2    <- eset_geo[[1]]
eset2    <- exprs(eset2)
eset2[1:5,1:5]

#注释信息
eset2<-anno_eset(eset       = eset2,
                 annotation = anno_hug133plus2,
                 symbol     = "symbol",
                 probe      = "probe_id",
                 method     = "mean")
eset2[1:5, 1:5]

#批次效应
eset_com <- remove_batcheffect( eset1       = eset1,  
                                eset2       = eset2,   
                                eset3       = NULL,
                                id_type     = "symbol",
                                data_type   = "array", 
                                cols        = "normal", 
                                palette     = "jama", 
                                log2        = TRUE, 
                                check_eset  = TRUE,
                                adjust_eset = TRUE,
                                repel       = FALSE,
                                path        = "result")

dim(eset_com)

####----肿瘤生态分析----####
library("GEOquery")
# NOTE: This process may take a few minutes which depends on the internet connection speed. Please wait for its completion.
eset_geo <- getGEO(GEO     = "GSE62254", getGPL  = F, destdir = "./")
eset    <- eset_geo[[1]]
eset    <- exprs(eset)
eset[1:5,1:5]

#基因注释
# Conduct gene annotation using `anno_hug133plus2` file; If identical gene symbols exists, these genes would be ordered by the mean expression levels. The gene symbol with highest mean expression level is selected and remove others. 
eset<-anno_eset(eset       = eset,
                annotation = anno_hug133plus2,
                symbol     = "symbol",
                probe      = "probe_id",
                method     = "mean")
eset[1:5, 1:3]

##胃癌TME分型
#安装TME分型包
if (!requireNamespace("TMEclassifier", quietly = TRUE)) devtools::install_github("IOBR/TMEclassifier")
library(TMEclassifier)
tme <- tme_classifier(eset = eset, scale = TRUE)
head(tme)

##差异分析
pdata <- tme[!tme$TMEcluster=="IS", ]
deg  <-   iobr_deg(eset         = eset,
                   annoation    = NULL,
                   pdata        = pdata,
                   group_id     = "TMEcluster",
                   pdata_id     = "ID",
                   array        = TRUE,
                   method       = "limma",
                   contrast     = c("deg_group","IA","IE"),
                   path         = NULL,
                   padj_cutoff  = 0.01,
                   logfc_cutoff = 0.5)

##GSEA分析
sig_list <- signature_collection[c("TMEscoreB_CIR", "TMEscoreA_CIR", "DNA_replication", "Base_excision_repair",
                                   "Pan_F_TBRs", "TGFb.myCAF", "Ferroptosis", "TLS_Nature", "Glycolysis")]
sig_list

gsea <-  sig_gsea(deg,
                  genesets          = sig_list,
                  path              = "GSEA",
                  gene_symbol       = "symbol",
                  logfc             = "log2FoldChange",
                  org               = "hsa",
                  show_plot         = TRUE,
                  msigdb            = TRUE,
                  category          = "H",
                  subcategory       = NULL,
                  palette_bar       = "set2")
gsea

####----差异分析，方法2----####
library(Seurat)
res <- find_markers_in_bulk(pdata      = tme, 
                            eset       = eset, 
                            group      = "TMEcluster", 
                            nfeatures  = 2000, 
                            top_n      = 20, 
                            thresh.use = 0.15, 
                            only.pos   = TRUE, 
                            min.pct    = 0.10)

top15 <-  res$top_markers %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(15, avg_log2FC)
top15$gene

#定义分型对应的颜色
cols <- c('#2692a4','#fc0d3a','#ffbe0b')
p1 <- DoHeatmap(res$sce, top15$gene, group.colors = cols )+
  scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)))

#提取表达矩阵
input <- combine_pd_eset(eset = eset, pdata = tme, feas = top15$gene, scale = T)
p2 <- sig_box(input, variable = "TMEcluster", signature = "IFNG", jitter = TRUE,
              cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4)

p3 <- sig_box(input, variable = "TMEcluster", signature = "IL1A",  
              jitter = TRUE, cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4)

#合并以上数据
if (!requireNamespace("patchwork", quietly = TRUE))   install.packages("patchwork")
library(patchwork)
p <- (p1|p2/p3) + plot_layout(widths = c(2.3,1))
p + plot_annotation(tag_levels = 'A')

####----Signature鉴定----####
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)
sig_tme <- t(column_to_rownames(sig_tme, var = "ID"))
sig_tme[1:5, 1:3]

#寻找变量特征
res <- find_markers_in_bulk(pdata = tme, eset = sig_tme, group = "TMEcluster", nfeatures = 1000, top_n = 20, min.pct = 0.10)

top15 <-  res$top_markers %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(15, avg_log2FC)

p1 <- DoHeatmap(res$sce, top15$gene, group.colors = cols)+
  scale_fill_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)))

#可视化结果
top15$gene  <- gsub(top15$gene, pattern = "\\-", replacement = "\\_")
input <- combine_pd_eset(eset = sig_tme, pdata = tme, feas = top15$gene, scale = T)

p2 <- sig_box(input, variable = "TMEcluster", signature = "IFNG_signature_Ayers_et_al", jitter = TRUE,
              cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4, size_of_font = 6)

p3 <- sig_box(input, variable = "TMEcluster", signature = "Neutrophils_Bindea_et_al",  
              jitter = TRUE, cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4, size_of_font = 6)

p <- (p1|p2/p3) + plot_layout(widths = c(2.3,1))
p + plot_annotation(tag_levels = 'A')

#生存分析
library(survminer)
data(pdata_acrg, package = "IOBR")
input <- merge(pdata_acrg, input, by = "ID")
p1<-surv_group(input_pdata       = input,
               target_group      = "TMEcluster",
               ID                = "ID",
               reference_group   = "High",
               project           = "ACRG",
               cols              = cols, 
               time              = "OS_time",
               status            = "OS_status",
               time_type         = "month",
               save_path         = "result")

p1

p2 <- percent_bar_plot(input, x = "TMEcluster" , y = "Subtype", palette = "jama")

p3 <- percent_bar_plot(input, x = "TMEcluster" , y = "Lauren", palette = "jama")

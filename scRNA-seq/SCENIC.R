# 不要轻易用自己的电脑尝试

#################
# 分析
# 安装和加载R包
# 下载基因注释配套数据库
library(Seurat)
library(tidyverse)
library(foreach)
library(RcisTarget)
library(doParallel)
library(SCopeLoomR)
library(AUCell)
BiocManager::install(c("doMC", "doRNG"))
library(doRNG)
BiocManager::install("GENIE3")
library(GENIE3)
#if (!requireNamespace("devtools", quietly = TRUE)) 
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
library(SCENIC)
#这里下载人的
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) 
}

# 构建分析数据
exprMat <- as.matrix(immune@assays$RNA@data)#表达矩阵
exprMat[1:4,1:4]#查看数据
cellInfo <- immune@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
#构建scenicOptions对象，接下来的SCENIC分析都是基于这个对象的信息生成的
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "F:/cisTarget_databases", nCores = 1)

# Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat.filtered <- exprMat[genesKept, ]
exprMat.filtered[1:4,1:4]
runCorrelation(exprMat.filtered, scenicOptions)
exprMat.filtered.log <- log2(exprMat.filtered + 1)
runGenie3(exprMat.filtered.log, scenicOptions)
# 非常费时间

# AUCell打分 也费时间
# Build and score the GRN
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
exprMat_log <- log2(exprMat + 1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions,exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

#################
# 可视化结果
# 可视化转录因子与seurat细胞分群联动
#regulons AUC
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(immune, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'immuneauc.rds')

#二进制regulo AUC
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(immune, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'immunebin.rds')

# FeaturePlot 或者小提琴图等方式可视化，参考GSVA分析

# 热图结束
library(pheatmap)
celltype = subset(cellInfo,select = 'CellType')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
regulons <- c('CEBPB_extended_1035g','RUNX1_extended_200g',
              'FOXC1_extended_100g','MYBL1_extended_112g',
              'IRF1_extended_1785g',
              'ELF1_1760g','ELF1_extended_2165g',
              'IRF1_extended_1765g','ETS1_extended_2906g',
              'YY1_extended_1453g','ATF3_extended_1022g',
              'E2F4_extended_637g',
              'KLF2_12g','HES1_13g',
              'GATA3_20g','HOXB2_extended_362g',
              'SOX4_extended_10g',
              'RUNX3_extended_170g','EGR3_extended_23g',
              'MXD4_extended_182g','HOXD9_extended_25g')
AUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%regulons,]
BINmatrix <- BINmatrix[rownames(BINmatrix)%in%regulons,]
pheatmap(AUCmatrix, show_colnames=F, annotation_col=celltype,
         width = 6, height = 5)
pheatmap(BINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 6, height = 5)
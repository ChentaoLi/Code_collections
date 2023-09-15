# Get the target subclusters
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
## 提取只含T细胞的子集
scRNA.T <- subset(scRNA,idents="Tc")
## 重新降维聚类
scRNA.T <- NormalizeData(scRNA.T, normalization.method = "LogNormalize", scale.factor = 1e4) 
scRNA.T <- FindVariableFeatures(scRNA.T, selection.method = 'vst', nfeatures = 3000)
scRNA.T <- ScaleData(scRNA.T,rownames(scRNA.T))
scRNA.T <- RunPCA(scRNA.T, features = VariableFeatures(object = scRNA.T)) 
library(harmony)
scRNA.T <- RunHarmony(scRNA.T, "orig.ident")
scRNA.T <- RunUMAP(scRNA.T,dims = 1:30, 
                   reduction = "harmony")
scRNA.T <- FindNeighbors(scRNA.T, reduction = "harmony",
                         dims = 1:30) 
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  scRNA.T=FindClusters(scRNA.T, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
apply(scRNA.T@meta.data[,grep("RNA_snn",colnames(scRNA.T@meta.data))],2,table)
Idents(scRNA.T)<- "RNA_snn_res.0.1"
p1 <- DimPlot(scRNA.T, reduction = "umap",
              label = T,label.box = T)+xlim(-15,15) 

# Monocle3
# ##创建CDS对象并预处理数据
data <- GetAssayData(scRNA.T, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA.T@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA.T, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype1") + ggtitle('int.umap')
## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE,color_cells_by = "celltype", group_label_size=4,cell_size=1.5)
# 确定root
cds = order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

# velocyto.R
#loom 文件生成
cellranger_gtf=$HOME/reference/scgenome.gtf
rmsk_gtf=$HOME/reference/Sus_repeat_rmsk.gtf  ## UCSC 版本确定
cellranger_outDir=/ifs/home/yusenwei/data/project/N1
velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf
R
library(velocyto.R)
library(Seurat)
library(SeuratWrappers)
scRNA.T <- readRDS("./scRNA.T.rds")
#载入 loom文件
N1= read.loom.matrices("./N1/velocyto/N2.loom")
N2= read.loom.matrices("./N2/velocyto/N2.loom")
N3= read.loom.matrices("./N3/velocyto/N3.loom")
R1= read.loom.matrices("./R1/velocyto/R1.loom")
R2= read.loom.matrices("./R2/velocyto/R2.loom")
R3= read.loom.matrices("./R3/velocyto/R3.loom")
#修改列名
N1$spliced[1:5, 1:5]
N1 <- lapply(N1, function(x) {
  colnames(x) <- sub("x", "-1", sub("N1", "N1_", colnames(x)))
  x
})
N2<- lapply(N2, function(x) {
  colnames(x) <- sub("x", "-1", sub("N2", "N2_", colnames(x)))
  x
})
N3<- lapply(N3, function(x) {
  colnames(x) <- sub("x", "-1", sub("N3", "N3_", colnames(x)))
  x
})
R1 <- lapply(R1, function(x) {
  colnames(x) <- sub("x", "-1", sub("R1", "R1_", colnames(x)))
  x
})
R2 <- lapply(IUGR2, function(x) {
  colnames(x) <- sub("x", "-1", sub("R2", "R2_", colnames(x)))
  x
})
R3 <- lapply(R3, function(x) {
  colnames(x) <- sub("x", "-1", sub("R3", "R3_", colnames(x)))
  x
})
ldat <- list(
  N1=N1,
  N2=N2,
  N3=N3,
  R1=R1,
  R2=R2,
  R3=R3
)
matrix.name <- names(ldat$NBW1)
ldat <- lapply(matrix.name, function(x){
  dat.list <- lapply(ldat, function(y){
    y[[x]]
  })
  dat.merged <- do.call(cbind, dat.list)
  dat.merged
})
names(ldat) <- matrix.name
cells.id <- colnames(scRNA.T)
cells.id[1:5]
ldat <- lapply(ldat, function(x) {
  x[, cells.id]
})
lapply(ldat, dim)
emat <- ldat$spliced
nmat <- ldat$unspliced
cell_cluster <- as.character(scRNA.T@meta.data$celltype1)
names(cell_cluster) <- cells.id
emat <- filter.genes.by.cluster.expression(emat, cell_cluster, min.max.cluster.average = .1)
nmat <- filter.genes.by.cluster.expression(nmat, cell_cluster, min.max.cluster.average = .1)
length(intersect(rownames(emat), rownames(nmat)))

# RNA Velocity analysis
### 参数设置
fit.quantile = 0.05 # 官方教程设定为 0.05
deltaT = 1 # default: 1
kCells = 10 # default: 10
### RNA velocity分析
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT = deltaT, kCells = kCells, fit.quantile = fit.quantile)


### 参数设定
emb = t(reducedDimS(monocle_cds)) # DDRTree坐标 monocle2
emb = Embeddings(scRNA.T, reduction = "umap") # monocle3 seurat
vel = rvel.qf # velocity estimates
n = 100 # 最邻近细胞的数量
scale = "sqrt" # scale方法
# 散点的颜色（用以区分不同的细胞状态）
cell.colors = plyr::mapvalues(cell_cluster, names(colors), colors)
cell.alpha = 0.2 # 散点颜色的透明度
cell.cex = 1 # 散点的尺寸
arrow.scale = 1 # 箭头的长度
arrow.lwd = 1.5 # 箭头的粗细
grid.n = 50 # grids的数量

### plot
pdf(file="velo.pdf")
show.velocity.on.embedding.cor(emb, vel, n, scale=scale, 
                               cell.colors = ac(cell.colors, alpha = cell.alpha),
                               cex = cell.cex, arrow.scale = arrow.scale, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 1,
                               grid.n = grid.n, arrow.lwd = arrow.lwd,cell.border.alpha = 0.1)
dev.off()

# scvelo-python
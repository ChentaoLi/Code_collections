# Monocle2
# http://cole-trapnell-lab.github.io/monocle-release/docs/
#  BiocManager::install("monocle")
library(monocle)#monocle构建CDS需要3个矩阵：expr.matrix、pd、featuredata
# Macro - seruat 对象
sample_ann <-  Macro@meta.data 
#构建featuredata，一般featuredata需要两个col，一个是gene_id,一个是gene_short_name,row对应counts的rownames
gene_ann <- data.frame(gene_short_name = rownames(Macro@assays$RNA),
                       row.names =  rownames(Macro@assays$RNA))
#head(gene_ann)
pd <- new("AnnotatedDataFrame",data=sample_ann)
fd <- new("AnnotatedDataFrame",data=gene_ann)
#构建matrix
ct=as.data.frame(Macro@assays$RNA@counts)#单细胞counts矩阵

#构建monocle对象
sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 筛选基因
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F)

# 数据降维
cds <- reduceDimension(cds, max_components = 2, num_dim = 20,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds))

# 将拟时与seurat分群对应，并挑选显著性基因可视化。
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$celltype)
pData(cds)$Cluster=pData(cds)$celltype
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Cluster")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
#  挑选差异最显著的基因可视化
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL)

# 前面差异基因筛选后，开始拟时推测。
# 第一步: 挑选合适的基因
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
#第二步降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
save(cds, "Monocle2.rds") # 保存
#可视化细胞分化轨迹
plot_cell_trajectory(cds, color_by = "Cluster")

# 可视化基因时序图。
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster")

# 先做一个细胞群的谱系分化图。从这个图可以看出我们关注的细胞分化轨迹。
library(ggsci)
plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm()
plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1) # 安装state分图

# 处理细胞谱系拟时可视化，我们还关注分化轨迹过程中基因的情况。选定关注的基因，看看其在拟时中的表达。
pData(cds)$TGFBR2 = log2( exprs(cds)['TGFBR2',]+1)

# 通过拟时基因表达模式聚类。
cds$id <- rownames(cds)
library(dplyr)
cds %>% arrange(qval) %>% head(10) %>% select(id) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$id
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds[gene_to_cluster,],
                                                 num_clusters = 3,
                                                 cores = 8,
                                                 show_rownames = TRUE)

# BEAM进行统计分析。
BEAM_res <- BEAM(my_cds_subset, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
table(BEAM_res$qval < 1e-4)
plot_genes_branched_heatmap(my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 8,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)

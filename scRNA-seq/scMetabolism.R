library(scMetabolism)
countexp.Seurat<-sc.metabolism.Seurat(obj = TSNE, method = "AUCell", imputation =F, ncores = 2, metabolism.type = "KEGG")
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score

row_means <- apply(metabolism.matrix, 1, mean)
sorted_means <- order(row_means, decreasing = TRUE)
top_50_means_rows <- sorted_means[1:40]
common_rows <- intersect(top_20_rows, top_50_means_rows)
common_rows_data <- metabolism.matrix[common_rows, ]

# 计算每行的标准差
row_sd <- apply(common_rows_data, 1, sd)
# 根据标准差对行进行排序
sorted_rows <- order(row_sd, decreasing = TRUE)
# 获取前20个变化最显著的行
top_20_rows <- sorted_rows[1:8]
# 提取原始矩阵中相应的行
top_20_rows_data <- common_rows_data[top_20_rows, ]

input.pathway<- rownames(top_20_rows_data)

DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "method", norm = "y")

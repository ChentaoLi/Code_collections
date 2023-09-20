# devtools::install_github("arc85/celltalker")
library(celltalker)
library(Seurat)
library(dplyr)
library(magrittr)

# 数据输入对象：一个Seurat对象（已经做好细胞分类）
## Run celltalker
#创建一个细胞类型向量，这里的pbmc.combined是一个已经定义好cell type的Seurat对象，其cell type name存在active.ident中，
cell_types <- pbmc.combined@active.ident
unique(cell_types)
# 将细胞类型向量添加到meta.data中
pbmc.combined <- AddMetaData(pbmc.combined, metadata = cell_types, col.name = "Cell_Type")
head(pbmc.combined@meta.data)

##运行celltalk函数---很快
pbmccombin_inter <- celltalk(input_object=pbmc.combined,
                             metadata_grouping="Cell_Type",
                             ligand_receptor_pairs=ramilowski_pairs,
                             number_cells_required=100,
                             min_expression=1000,
                             max_expression=20000,
                             scramble_times=10)    #函数参数介绍见下图，ramilowski_pairs这一配体-受体对只适用于人，要进行鼠的分析需进行小写转换
pbmccombin_inter  


## Identify top statistically significant interactions—识别具有统计显著性的顶级互作，并进行可视化。
top_stats <- pbmccombin_inter  %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05) %>%
  group_by(cell_type1) %>%
  top_n(3,interact_ratio) %>%
  ungroup()


#配色后绘制圈图
library(RColorBrewer)   #加载RColorBrewer包，该包提供了一系列优美的颜色调色板
n_colors <- max(13, length(unique(pbmc.combined$Cell_Type)))  # 依据细胞类型最少需要 13 种颜色
palette <- rep(brewer.pal(n = 8, name = "Set2"), length.out = n_colors)  # 重复使用 Set2 调色板，使用brewer.pal函数从Set2调色板中选择8种颜色，并将其重复使用以满足所需的颜色数量
colors_use <- palette[1:n_colors]  # 选择需要的颜色数量
circos_plot(ligand_receptor_frame=top_stats,   #顶级互作框
            cell_group_colors=colors_use,
            ligand_color="blue",
            receptor_color="red",
            cex_outer=0.5,
            cex_inner=0.4)

# References: https://zhuanlan.zhihu.com/p/628699666?utm_id=0
######################################################################
DimPlot(hca_bm,group.by="cell_types")
## Run celltalker
hca_bm_interactions <- celltalk(input_object=hca_bm,
                                metadata_grouping="cell_types",
                                ligand_receptor_pairs=ramilowski_pairs,
                                number_cells_required=100,
                                min_expression=1000,
                                max_expression=20000,
                                scramble_times=10)

## Identify top statistically significant interactions
top_stats <- hca_bm_interactions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05) %>%
  group_by(cell_type1) %>%
  top_n(3,interact_ratio) %>%
  ungroup()

## Generate a circos plot
colors_use <- RColorBrewer::brewer.pal(n=length(unique(hca_bm$cell_types)),"Set2")

circos_plot(ligand_receptor_frame=top_stats,
            cell_group_colors=colors_use,
            ligand_color="blue",
            receptor_color="red",
            cex_outer=0.5,
            cex_inner=0.4)
# https://github.com/arc85/celltalker
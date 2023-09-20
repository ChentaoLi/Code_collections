# 使用CellChat 对多个数据集细胞通讯进行比较分析
#######################
# 准备工作
library(CellChat)
library(patchwork)
# 创建目录以保存图片
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

# 加载每个数据集的cellchat对象，然后合并在一起
# 用户需要在每个数据集上单独运行 CellChat，然后将不同的 CellChat 对象合并在一起。
cellchat.NL <- readRDS()
cellchat.LS <- readRDS()
object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#######################
# 第一部分：预测细胞通信的一般原理
# 细胞-细胞通信是否增强
# 细胞类型显著变化之间的相互作用
# 主要来源和目标如何从一个条件到为另一个条件变化的

# 比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# 比较不同细胞群之间的相互作用数量和强度
# 两个数据集之间细胞通信网络中交互或交互强度的差异数可以使用圆图可视化， 与第一个数据集相比，[红色]（或[蓝色]边表示信号在第二个数据集中增加或[减少]）。
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# 热图
# 顶部彩色条形图表示热图（传入信号）中显示的列值的总和。右边的彩色条形图表示一行值（传出信号）的总和。在色条中红色或蓝色表示第二个数据集中与第一个数据集相比增加或[减少]信号。
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# 差异网络分析仅适用于配对数据集。如果有更多的数据集进行比较，我们可以直接显示每个数据集中任意两个细胞群之间的交互次数或交互强度。
# 为了更好地控制不同数据集中推断网络的节点大小和边缘权重，我们计算每个细胞组的最大细胞数量以及所有数据集中交互（或交互权重）的最大数量。
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# 不同细胞类型之间相互作用或交互强度的差异
# 为了简化复杂的网络，并深入了解细胞类型级别的细胞通信，我们可以根据定义的细胞群聚合细胞-细胞通信。在这里，我们将细胞群分为三种细胞类型，然后重新合并CellChat对象列表。
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
# 然后，我们可以显示每个数据集中任意两个细胞类型之间的交互次数或交互强度。
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# 此外，我们还可以使用圆图显示任意两种细胞类型之间的交互或交互强度的差异。与第一个数据集相比，红色（或蓝色）色边缘表示第二个数据集中的信号增加（或减少）。
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

# 比较 2D 空间中的主要来源和目标
# 比较二D空间中的传出和传入交互强度，可以识别不同数据集之间显著变化的发送或接收信号的细胞群。
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#######################
# 第二部分：识别保守和环境特异的信号通路
# CellChat 可以根据其在多种生物条件下的细胞通信网络，识别差异较大（或更少）的信号网络、信号组以及基于其细胞通信网络的保守和环境特异的信号通路。
# 根据信号/结构的相似性识别差异较大（或更少）的信号网络以及信号组
# CellChat 根据推断的通信网络的功能和拓扑相似性，对其进行联合多重学习和分类。NB：此类分析适用于两个以上的数据集。
# 功能相似性：功能相似度高表示主要发射器和接收器相似，可解释为两个信号通路或两个配体受体对具有相似的作用。NB：功能相似性分析不适用于具有不同细胞类型成分的多个数据集。
# 结构相似性：结构相似性用于比较其信号网络结构，而不考虑发送器和接收器的相似性。NB：结构相似性分析适用于具有相同细胞类型组成或截然不同的细胞类型组成多个数据集。
# 在这里，我们可以根据功能相似性运行多重和分类学习分析，因为两个数据集具有相同的单元类型组成。

# 根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# 基于结构相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2

# 计算和可视化通路距离
# 我们可以根据信号网络在共享双维空间中的欧几里德距离来识别差异较大（或更少）的信号网络。更大的距离意味着两个数据集之间的通信网络在功能或结构相似性方面存在更大的差异。NB：我们只计算两个数据集之间重叠信号通路的距离。此处未考虑仅在一个数据集中标识的信号通路。如果有三个以上的数据集，可以通过在函数rankSimilarity中定义comparison进行对比。
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

# 识别并可视化保守和环境特异的信号通路
# 通过比较每个信号通路的信息流/交互，我们可以识别信号通路，（i） 关闭，（ii） 减少，（iii） 打开或（iv） 增加。
# 比较每个信号通路的整体信息流
# 我们可以通过简单地比较每个信号通路的信息流来识别保守和环境特异的信号通路，该信息流由推断网络中所有一对细胞群之间的通信概率之和（即网络中的总权重）定义。
# 此条形图可在堆叠模式下绘制。根据 NL 和 LS 皮肤之间推断的网络中整体信息流的差异对重要信号通路进行排名。红色的顶部信号通路富含 NL 皮肤，绿色在 LS 皮肤中得到了富集。
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# 比较与每个细胞群相关的传出（或传入）信号
# 上述分析将传出和传入信号的信息汇总在一起。我们还可以比较两个数据集之间的传出（或传入）信号模式，从而识别显示不同信号模式的信号通路/受配体。
# 我们可以将来自不同数据集的所有已识别的信号通路进行组合，从而并排比较它们，包括传出信号、传入信号和整体信号，方法是将传出和传入信号聚合在一起。NB：rankNet还显示了整体信号的比较，但它没有显示特定细胞群中的信号强度。
library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##############################3
# 第三部分：识别上调和下调的信号配体对
# 我们可以比较由某些细胞群到其他细胞组的配体受体对调节的通信概率。这可以通过设置comparison在函数netVisual_bubble中来完成。
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

# 此外，我们可以在一个数据集中识别与另一个数据集相比，上升（增加）和下降调节（减少）信号配体受体对。这可以通过指定max.dataset和min.dataset在函数netVisual_bubble中完成。信号增加意味着这些信号在一个数据集中与其他数据集相比具有更高的通信概率（强度）。
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

# NB：气泡图中显示的配体受体对可以通过signaling.LSIncreased = gg1$data 访问。
# 通过比较每个 L-R 对和每对细胞组的两个数据集之间的通信概率，可以采用上述方法来识别上调和下调的信号。另外，我们可以根据微分基因表达分析来识别上调和下调的信号配体对。具体来说，我们对每个细胞组执行两种生物条件（即NL和LS）之间的微分表达分析，然后根据发送者细胞中配体和接收器细胞中受体的折叠变化获得上调和下调的信号。此类分析可如下所示。
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
# 由于信号基因在多亚单位中可能很复杂，我们可以使用net.upnet.down进一步的来获得单个信号基因。

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# 然后，我们使用气泡图或和弦图可视化上调和向下调的信号配体对。
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# 使用和弦图可视化上调和下调的信号配体对
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

##############################################################
# 第四部分：使用层次结构图、圆图或和弦图可视比较细胞-细胞通信
# 与单个数据集的 CellChat 分析类似，我们可以使用层次结构图、圆图或和弦图可视化细胞通信网络。
# 边缘颜色/重量、节点颜色/大小/形状：在所有可视化图中，边缘颜色与发送者源一致，边缘权重与交互强度成正比。较厚的边缘线表示信号更强。在层次结构图和圆图中，圆的大小与每个细胞组中的细胞数量成正比。在层次图中，实心和开放的圆分别代表源和目标。在和弦图中，内条颜色表示从相应的外条接收信号的目标。内条大小与目标接收的信号强度成正比。这种内条有助于解释复杂的和弦图。请注意，有一些内条没有任何和弦的一些细胞组，请忽略他们，因为这是一个本包尚未解决的问题。

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# netVisual_chord_cell对于和弦图，CellChat 具有独立函数，通过调整circlize包中的不同参数来灵活可视化信号网络。例如，我们可以定义一个group命名的字符矢量，以创建多组和弦图，将细胞群集分组到不同的细胞类型。

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# 使用和弦图，CellChat 提供两个函数netVisual_chord_cell和netVisual_chord_gene，可视化具有不同目的和不同级别的细胞-细胞通信。 netVisual_chord_cell用于可视化不同细胞群之间的细胞-细胞通信（和弦图中的每个部分是细胞组），netVisual_chord_gene用于可视化由多个配体受体或信号通路调解的细胞-细胞通信（和弦图中的每个部分都是配体、受体或信号通路）。
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'MIF'.
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'CXCL '.
# NB：在生成绘图时，请忽略注释，例如"Note: The first link end is drawn out of sector ‘MIF’"。如果基因名称重叠，您可以通过降低small.gap值来调整参数。

#####################################################
# 第五部分：比较不同数据集之间的信号基因表达分布
# 我们可以利用seurat包装的函数plotGeneExpression绘制与L-R对或信号通路相关的信号基因的基因表达分布图。
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.

# 保存合并的cellchat对象
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_humanSkin_NL_vs_LS.rds")

# https://www.jianshu.com/p/49a0a0b50987







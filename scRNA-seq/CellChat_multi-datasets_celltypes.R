# 使用CellChat 对多个数据集细胞通讯进行比较分析
#######################
# 准备工作
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)

#######################
# 第一部分：比较分析具有略有不同细胞类型成分的多个数据集
# 对于具有稍微不同的细胞类型（组）组成的数据集，CellChat 可以使用函数liftCellChat将细胞组提升到所有数据集的相同细胞标记，然后执行比较分析，作为对具有相同细胞类型成分的数据集的联合分析。
# 加载每个数据集的CellChat对象
# 用户需要在每个数据集上单独运行 CellChat，然后将不同的 CellChat 对象合并在一起。在这里，我们也使用updateCellChat这样做，因为这两个对象是使用早期版本（<0.5.0）的CellChat获得。
cellchat.E13 <- readRDS(url("https://ndownloader.figshare.com/files/25957094"))
cellchat.E13 <- updateCellChat(cellchat.E13)
#> Update slot 'var.features' from a vector to a list
cellchat.E14 <- readRDS(url("https://ndownloader.figshare.com/files/25957634"))
cellchat.E14 <- updateCellChat(cellchat.E14)
#> Update slot 'var.features' from a vector to a list

# 升级cellchat对象并合并在一起
# Define the cell labels to lift up
group.new = levels(cellchat.E14@idents)
cellchat.E13 <- liftCellChat(cellchat.E13, group.new)
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
object.list <- list(E13 = cellchat.E13, E14 = cellchat.E14)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
#> Warning in mergeCellChat(object.list, add.names = names(object.list),
#> cell.prefix = TRUE): Prefix cell names!
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

# 使用提升的对象可视化推断的信号网络
# Hierarchy plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Circle plot
pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#> Note: The second link end is drawn out of sector 'Basal-P'.
#> Note: The second link end is drawn out of sector 'FIB-A'.

# 对于和弦图，CellChat 具有独立函数netVisual_chord_cell，通过调整circlize包中的不同参数来灵活可视化信号网络。例如，我们可以定义一个group命名的字符矢量，以创建多组和弦图，例如，将细胞群集分组到不同的细胞类型。
# Chord diagram
group.merged <- c(rep("Dermal", 10), rep("Epidermal", 3)) # grouping cell clusters into dermal and epidermal cells to study the cell-cell communication between dermal and epidermal
names(group.merged) <- levels(object.list[[1]]@idents)
pathways.show <- c("WNT") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.merged, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The second link end is drawn out of sector 'Basal-P'.
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The second link end is drawn out of sector 'FIB-A'.

#######################
# 第二部分：对具有截然不同的细胞类型成分的多个数据集的比较分析
# CellChat 可用于比较来自截然不同的生物背景的两个 scRNA-seq 数据集之间的细胞-细胞通信模式。如胚胎形态形成细胞与引起伤口修复的细胞。对于具有截然不同的细胞类型（组）组成的数据集，除了以下两个方面外，大多数 CellChat 的功能都可以应用：
# 不能用于比较不同细胞群之间相互作用的差异数和相互作用强度。
# 但是，用户仍然可以使用函数 netVisual_diffInteraction和netVisual_circle来显示交互次数和交互强度。``
# 不能使用computeNetSimilarityPairwise(cellchat, type = "functional")的功能相似性识别信号组。


# References: https://www.jianshu.com/p/5743e29ffd03
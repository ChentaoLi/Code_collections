##########################
# 安装包
devtools::install_github("sqjin/CellChat")

##########################
# 载入包
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

##########################
# 第一部分：CellChat对象的数据输入、处理及初始化
## CellChat 需要两个输入：一个是细胞的基因表达数据，另一个是用户分配的细胞标签（即基于标签的模式）或单细胞数据的低维表示（即无标签模式）。对于后者，CellChat 通过根据低维空间或伪时间轨迹空间中的细胞距离构建共享的邻近图自动对细胞进行分组。

# 加载数据
## 对于基因表达数据矩阵，要求基因为行名，细胞为列名。需要将标准化数据作为 CellChat 分析的输入。如果用户提供count数据，我们提供一个函数normalizeData来计算文库大小，然后进行log转换。对于分组信息，需要使用带有行名的数据作为CellChat 的输入。
# Here we load a scRNA-seq data matrix and its associated cell meta data
load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
#>  [1] Inflam. FIB  FBN1+ FIB    APOE+ FIB    COL11A1+ FIB cDC2        
#>  [6] LC           Inflam. DC   cDC1         CD40LG+ TC   Inflam. TC  
#> [11] TC           NKT         
#> 12 Levels: APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 ... NKT

# 创建CellChat 对象
## 用户可以从数据矩阵、Seurat 或SingleCellExperiment对象创建新的 CellChat 对象。如果输入是 Seurat 或SingleCellExperiment对象，则默认情况下将使用对象中的meta data，用户必须提供该数据来定义细胞分组。例如，group.by=Seurat 对象中默认的细胞标识。
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC


# 将细胞信息添加到对象的meta slot中
## 如果在创建cellchat对象时未添加细胞meta信息，用户也可以稍后添加该信息，并使用setIdent设置该对象默认的细胞标识。
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# 设置配体受体交互数据库
## 我们的数据库 CellChatDB 是一个手动整理的文献支持的配体受体在人和小鼠中的交互数据库。小鼠中的CellChatDB包含2，021个经验证的分子相互作用，包括60%的自分泌/旁分泌信号相互作用、21%的细胞外基质（ECM）受体相互作用和19%的细胞-细胞接触相互作用。人的CellChatDB包含1，939个经验证的分子相互作用，包括61.8%的自分泌/旁分泌信号相互作用、21.7%的细胞外基质（ECM）受体相互作用和16.5%的细胞-细胞接触相互作用。
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# 预处理用于细胞通信分析的表达数据
## 为了推断细胞状态特异的通信，我们识别一个细胞组中过度表达的配体或受体，然后识别过度表达的配体受体相互作用，是否过表达。
## 我们还提供将基因表达数据投影到蛋白质-蛋白质相互作用 （PPI） 网络上的功能。具体来说，投影过程根据高度可信的实验验证的蛋白质-蛋白质网络中定义的基因表达值来平滑基因的表达值。此功能在分析具有浅测序深度的单细胞数据时很有用，因为投影可减少信号基因的dropput效应，特别是对于配体/受体的可能的零表达。人们可能担心这种投影过程可能引入的人为因素，但是，这几乎可以忽略不计。用户还可以通过在函数computeCommunProb()中国设置raw.use = TRUE跳过此步骤。

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

##########################
# 第二部分：细胞通信网络的推断
## CellChat 通过分配具有概率值的每个相互作用并进行排列检验，来推断具有生物学意义的细胞-细胞通信。CellChat通过将基因表达与先前已知的信号配体、受体及其同因子之间的相互作用知识相结合，利用大量作用规律，对细胞-细胞通信的概率进行模拟。
## 推断的配体受体对的数量显然取决于计算每个细胞组平均基因表达的方法。默认情况下，CellChat 使用一种统计学上强大的均值方法，称为"trimean"，与其他方法相比，它产生的相互作用更少。然而，我们发现 CellChat 在预测更强的交互方面表现良好，这非常有助于缩小交互范围，以便进一步进行实验验证。在computeCommunProb中，我们提供了一个选项，用于使用其他方法，如5%和10%截断均值，来计算平均基因表达。值得注意的是，"trimean"大约是25%的截断平均值，这意味着如果一组表达细胞的百分比低于25%，则平均基因表达为零。要使用 10% 截断的平均值，用户可以设置type = "truncatedMean"和对trim = 0.1。
## 在分析未分类的单细胞转录组时，假设丰富的细胞群倾向于发送比稀有细胞群更强的信号，CellChat 还可以在概率计算中考虑每个细胞组中细胞比例的影响。用户可以设置population.size = TRUE

# 计算通信概率并推断cellchat网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 提取推断的cellchat网络作为数据框架
# 我们提供一个函数subsetCommunication，以轻松访问推断感兴趣的细胞-细胞通信。例如
# df.net <- subsetCommunication(cellchat)返回一个数据框架，该数据框架由配体/受体级别的所有推断细胞通信组成。设置slot.name = "netP"可以在信号通路级别访问推断的通信
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))将推断的细胞-细胞通信从细胞组1和2发送到细胞组4和5。
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))通过向WNT和TGFb发出信号来调节推断的细胞通信。
# 在信号通路级别推断细胞-细胞通信
## CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率。
## NB：每个配体受体对和每个信号通路的推断细胞间通信网络分别存储在插槽"net"和"netP"中。
cellchat <- computeCommunProbPathway(cellchat)

# 计算整合的细胞通信网络
## 我们可以通过计算链接数或汇总通信概率来计算整合的细胞通信网络。用户还可以通过设置sources.use和targets.use`
cellchat <- aggregateNet(cellchat)
# 我们还可以可视化整合的细胞通信网络。例如，使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度（比重）。
groupSize <- as.numeric(table(cellchat@idents))par(mfrow = c(1,2), xpd=TRUE)netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# 由于细胞通信网络复杂，我们可以检查每个细胞组发送的信号。在这里，我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边缘权重。
mat <- cellchat@net$weightpar(mfrow = c(3,4), xpd=TRUE)for (i in 1:nrow(mat)) {  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  mat2[i, ] <- mat[i, ]  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}

##########################
# 第三部分：细胞通信网络的可视化
# 在推断细胞通信网络后，CellChat 为进一步的数据探索、分析和可视化提供了各种功能。
# 它提供了几种可视化细胞通信网络的方法，包括分层图、圆图、和弦图和气泡图。
# 它提供了一个易于使用的工具，用于提取和可视化推断网络的高阶信息。例如，它允许对细胞群的主要信号输入和输出以及这些群和信号如何协调功能进行现成预测。
# 它可以通过结合通讯网络分析、模式识别和多重学习方法，使用综合方法对推断出的细胞-细胞通信网络进行定量表征和比较。

# 使用层次结构图、圆图或和弦图可视化每个信号通路
pathways.show <- c("CXCL") # Hierarchy plot# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells vertex.receiver = seq(1,4) # a numeric vector. netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle 
plotpar(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord 
diagrampar(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
# Heatmap
Heatmappar(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")#> Do heatmap based on a single object
# Chord 
diagramgroup.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cellsnames(group.cellType) <- levels(cellchat@idents)netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))#> Plot the aggregated cell-cell communication network at the signaling pathway level#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# 计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
netAnalysis_contribution(cellchat, signaling = pathways.show)
# 我们还可以可视化由单个配体受体对调节的细胞-细胞通信。我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair# Hierarchy plotvertex.receiver = seq(1,4) # a numeric vectornetVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle 
plotnetVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord 
diagramnetVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# 自动保存所有推断网络的模块以进行快速探索
## 在实际使用中，用户可以使用‘for … loop’自动保存所有推断网络快速探索使用。 netVisual，netVisual支持svg、png和pdf格式的输出。
# Access all the signaling pathways showing significant communicationspathways.show.all <- cellchat@netP$pathways# check the order of cell identity to set suitable vertex.receiverlevels(cellchat@idents)vertex.receiver = seq(1,4)for (i in 1:length(pathways.show.all)) {  # Visualize communication network associated with both signaling pathway and individual L-R pairs  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)}
# 可视化由多个配体受体或信号通路调节的细胞通信
# 气泡图
# 我们还可以使用netVisual_bubble显示从某些细胞组到其他细胞组的所有重要相互作用（L-R对）。
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) associated with certain signaling 
pathwaysnetVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)#> Comparing communications on a single object

# 和弦图
# 类似于气泡图，cellchat提供了绘制和弦图的功能netVisual_chord_gene
# 显示从某些细胞组到其他细胞组的所有相互作用（L-R对或信号通路）。两个特殊情况：一个显示从一个细胞组发送的所有交互，另一个显示一个细胞组接收的所有交互。
# 显示用户输入的交互或用户定义的某些信号通路
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')# show all the interactions sending from Inflam.FIBnetVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)#> Note: The first link end is drawn out of sector 'MIF'
# show all the interactions received by Inflam.DCnetVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathwaysnetVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)#> Note: The second link end is drawn out of sector 'CXCR4 '.#> Note: The first link end is drawn out of sector 'CXCL12 '.
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'MIF'.#> Note: The second link end is drawn out of sector ' '.#> Note: The first link end is drawn out of sector 'CXCL '.
# NB：在生成图时，请忽略注释，例如 “Note: The first link end is drawn out of sector ‘MIF’。如果基因名称重叠，您可以通过降低small.gap值来调整参数。

# 使用小提琴/点图绘制信号基因表达分布
# 我们可以利用Seurat 包装的函数plotGeneExpression绘制与L-R对或信号通路相关的信号基因的基因表达分布图。
# plotGeneExpression(cellchat, signaling = "CXCL")#> Registered S3 method overwritten by 'spatstat':#>   method     from#>   print.boxx cli#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.#> Scale for 'y' is already present. Adding another scale for 'y', which will#> replace the existing scale.
# 默认情况下，用户可以通过plotGeneExpression只显示与推断的重要通信相关的信号基因的表达。
# plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
# 或者，用户可以使用extractEnrichedLR提取与推断的L-R对或信号通路相关的信号基因，然后使用Seurat包绘制基因表达图。


##########################
# 第四部分：细胞通信网络系统分析
## 为了便于对复杂的细胞间通信网络进行解释，CellChat 通过从图形理论、模式识别和多重学习中抽象的方法对网络进行量化。
## 它可以使用网络分析的集中度措施确定给定信号网络中的主要信号源和目标以及调节者和影响者
## 它可以通过利用模式识别方法预测特定细胞类型的关键传入和传出信号，以及不同细胞类型之间的协调响应。
## 它可以通过定义相似度测量方法和从功能和拓扑角度进行多重学习来分组信号通路。
## 它可以通过对多个网络的联合多重学习来描绘保存上下文特定的信号通路。
## 识别细胞组的信号角色（例如，占主导地位的发送器、接收器）以及主要贡献信号
## CellChat 允许通过计算每个细胞组的多个网络中心测量，随时识别细胞间通信网络中占主导地位的发送者、接收者、调解者和影响者。具体来说，我们在加权导向网络中采用了措施，包括度外、度内、介于两者之间流动和信息集中度，分别识别细胞间通信的主要发送者、接收者、调解者和影响者。在以权重为计算通信概率的加权定向网络中，将外向度计算为来自细胞组的传出信号的通信概率之和，并计算为传入信号对单元组通信概率的总和的度内，可用于分别识别信号网络的主要单元件发送器和接收器。有关信息中心之间流动的定义，请查看本文及相关参考文献。

# 计算和可视化网络中心分数
# Compute the network centrality 
scorescellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groupsnetAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# 在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
## 我们还提供了另一种直观方法，使用散点图在 2D 空间中可视化占主导地位的发射器（源）和接收器（目标）。
# Signaling role analysis on the aggregated cell-cell communication network from all signaling 
pathwaysgg1 <- netAnalysis_signalingRole_scatter(cellchat)#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways# Signaling role analysis on the cell-cell communication networks of interestgg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))#> Signaling role analysis on the cell-cell communication network from user's inputgg1 + gg2

# 识别对某些细胞组的传出或传入信号贡献最大的信号
## 我们还可以回答以下问题：哪些信号对某些细胞组的传出或传入信号贡献最大。
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathwaysht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interestht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

# 确定全局通信模式，探索多个细胞类型和信号通路如何协调在一起
## 除了探索单个通路的详细通信外，一个重要问题是多个细胞组和信号通路如何协调功能。CellChat 采用模式识别方法识别全局通信模式。
## 随着模式数量的增加，可能会有多余的模式，因此很难解释通信模式。我们选择了五种模式作为默认模式。一般来说，它具有生物学意义，模式数量要大于2。此外，我们还提供了一个函数selectK来推断模式的数量，该数基于 NMF R 包中已实施的两个指标Cophenetic和Silhouette。这两个指标都根据共识矩阵的分层聚类来衡量特定数量模式的稳定性。对于一系列模式，适当的模式数量是Cophenetic 和 Silhouette值开始突然下降的模式。
#识别和可视化分泌细胞的传出通信模式
## 传出模式揭示了发送者细胞（即作为信号源的细胞）如何相互协调，以及它们如何与某些信号通路协调以驱动通信。
## 为了直观地显示潜在模式与细胞群和配体受体对或信号通路的关联，我们使用了河流（冲积）图。我们首先将每行 W 和 H 的每列标准化为 [0，1]，然后在 W 和 H 中设置为零，如果它们小于 0.5。这种阈值允许发现与每个推断模式相关的最丰富的细胞组和信号通路，即每个细胞组或信号通路仅与一个推断模式相关联。这些阈值矩阵 W 和 H 用作创建冲积图的输入。
## 为了将细胞群与其丰富的信号通路直接联系起来，如果 W 和 H 中的元素少于 1/R（R 是潜在模式数），则我们将它们中的元素设置为零。通过使用不太严格的阈值，可以获得与每个细胞组相关的更丰富的信号通路。我们使用每个细胞组对通过乘以 W 乘以 H 计算的每个信号通路的贡献分数，构建了一个点图，其中点大小与贡献分数成正比，以显示细胞组与其丰富信号通路之间的关联。用户还可以降低参数cutoff，以显示每个细胞组关联的更丰富的信号通路。

# 通信模式分析所需的包
library(NMF)#> Loading required package: pkgmaker#> Loading required package: registry#> Loading required package: rngtools#> Loading required package: cluster#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16#>   To enable shared memory capabilities, try: install.extras('#> NMF#> ')#> #> Attaching package: 'NMF'#> The following objects are masked from 'package:igraph':#> #>     algorithm, comparelibrary(ggalluvial)
# 在这里，我们运行selectK推断模式的数量。
selectK(cellchat, pattern = "outgoing")
# 当传出模式数为 3 时，Cophenetic 和Silhouette值都开始突然下降。
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plotnetAnalysis_river(cellchat, pattern = "outgoing")#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plotnetAnalysis_dot(cellchat, pattern = "outgoing")

# 识别和可视化目标细胞的传入通信模式
# 传入模式显示目标细胞（即信号接收器中的细胞）如何相互协调，以及它们如何与某些信号通路协调以响应传入的信号。
selectK(cellchat, pattern = "incoming")
# 当传入模式的数量为 4 时，Cophenetic 值开始下降。
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plotnetAnalysis_river(cellchat, pattern = "incoming")#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plotnetAnalysis_dot(cellchat, pattern = "incoming")

# 信号网络的多重和分类学习分析
## 此外，CellChat 能够量化所有重要信号通路之间的相似性，然后根据其CellChat 网络的相似性对其进行分组。分组可以基于功能或结构相似性进行。
## 功能相似性：功能相似度高表示主要发送器和接收器相似，可解释为两个信号通路或两个配体受体对具有相似的作用。功能相似性分析要求两个数据集之间的细胞群组成相同。
## 结构相似性：结构相似性用于比较其信号网络结构，而不考虑发送器和接收器的相似性。

# 根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "functional")cellchat <- netEmbedding(cellchat, type = "functional")#> Manifold learning of the signaling networks for a single datasetcellchat <- netClustering(cellchat, type = "functional")#> Classification learning of the signaling networks for a single dataset# Visualization in 2D-spacenetVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
基于结构相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "structural")cellchat <- netEmbedding(cellchat, type = "structural")#> Manifold learning of the signaling networks for a single datasetcellchat <- netClustering(cellchat, type = "structural")#> Classification learning of the signaling networks for a single dataset# Visualization in 2D-spacenetVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

##########################
# 第五部分：保存cellchat对象
saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")



























































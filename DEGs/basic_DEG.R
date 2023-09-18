################## 差异分析三巨头 —— DESeq2, edgeR 和 limma 包 #################

# 加载要用到的包，大家如果没有可以先自己安装一下
library(data.table)  # 用于高效处理大数据集
library(dplyr)       # 用于数据操作和转换
library(ggplot2)     # 画图图
library(pheatmap)    # 绘制热图
library(DESeq2)      # 差异分析一号选手
library(edgeR)       # 差异分析二号选手
library(limma)       # 差异分析三号选手
library(tinyarray)   # 用于绘制各种图表，今天用它绘制韦恩图

############################ 表达矩阵和分组信息整理 ##############################
# 读取表达矩阵
exp_brca <- fread("./data/TCGA-BRCA.htseq_counts.tsv.gz", header = T, sep = '\t', data.table = F)

# ensembl id 转换为 gene symbol
gene_id <- fread("./data/gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)
gene_id <- gene_id[ , c(1, 2)]
exp_brca <- merge(gene_id, exp_brca, by.y  = "Ensembl_ID", by.x = "id" )

# 去重
exp_brca <- distinct(exp_brca, gene, .keep_all = T)

# 上面的操作是如果存在重复行，直接保留第一行
# 我们也可以对重复基因名取平均表达量再去重

# 可以用limma包中的avereps函数，或者也可以使用aggregate函数取平均
# library(limma)
# exp_brca <- avereps(exp_brca, exp_brca$gene)

# 把基因名转换为行名
rownames(exp_brca) <- exp_brca$gene
exp_brca <- exp_brca[ , -c(1,2)]
dim(exp_brca) # 58387  1217

# 分组信息
# TCGA数据通过样本ID就可以进行分组
# https://mp.weixin.qq.com/s/8Q7Rb7K-Mrs_fA5c5OvyqA

# 选取样本名14和15位置元素，因为它们可以代表样本类型
gp <- substring(colnames(exp_brca), 14, 15)

# 分组
brca_tumor <- exp_brca[, as.numeric(gp) < 10]
brca_normal <- exp_brca[, as.numeric(gp) >= 10]

# 按顺序存储在一个矩阵中（可以方便后面画热图）
exp_brca <- cbind(brca_tumor, brca_normal)

group <- c(rep('tumor', ncol(brca_tumor)), rep('normal', ncol(brca_normal)))
group <- factor(group, levels = c("normal", "tumor"))

########################### 使用 DESeq2 进行差异分析 ###########################
library(DESeq2)

# 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
colData <- data.frame(row.names = colnames(exp_brca),
                      group = group)
colData$group <- factor(colData$group, levels = c("normal", "tumor"))
# DESeq2 要求输入数据是由整数组成的矩阵，且没有经过标准化
# 我们在Xena下载的数据是 log2(count+1)，所以需要进行处理
exp_brca_int <- 2^(exp_brca) - 1
exp_brca_int <- apply(exp_brca_int, 2, as.integer)
rownames(exp_brca_int) <- rownames(exp_brca)

# 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
dds <- DESeqDataSetFromMatrix(countData = exp_brca_int, # 表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息

# 进行差异表达分析
dds <- DESeq(dds)

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
# res <- results(dds, contrast = c("group", levels(group)[2], levels(group)[1]))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
DEG_deseq2 <- na.omit(DEG)

# 将处理后的差异表达结果保存为R数据文件
save(DEG_deseq2, file = './data/DEG_deseq2.Rdata')

# 画图
load("./data/DEG_deseq2.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))

# 火山图
p <- ggplot(data = DEG_deseq2, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_deseq2.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
deg_opt <- DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_deseq2.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

########################### 使用 edgeR 进行差异分析 ############################

library(edgeR)

# 创建 DGEList 对象，用于存储基因表达数据和组信息，还是使用原始计数矩阵
d <- DGEList(counts = exp_brca_int, group = group)

# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2

# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)

# edgeR包中的 filterByExpr() 函数，它提供了自动过滤基因的方法，可保留尽可能多的有足够表达计数的基因。
# 此函数默认选取最小的组内的样本数量为最小样本数，保留至少在这个数量的样本中有10个或更多序列片段计数的基因。 
# 过滤标准是，以最小对组内样本数为标准，（此例最小组内样本为3），如果有基因在所有样本中表达数（count）
# 小于10的个数超过最小组内样本数，就剔除该基因。

# 这里补充解释一下 CPM 的含义。
# 常用的标度转换有 CPM（counts per million）、log-CPM、FPKM、RPKM 等。
# CPM 是将 Counts 转变为 counts per million，消除测序深度影响。
# log-CPM 是将 CPM 进行 log2 计算。cpm函数会在进行 log2 转换前给 CPM 值加上一个弥补值。
# 默认的弥补值是 2/L，其中 2 是“预先计数”，而 L 是样本总序列数（以百万计）的平均值，
# 所以 log-CPM 值是根据 CPM 值通过 log2(CPM + 2/L) 计算得到的。 
# 对于一个基因，CPM 值为 1 相当于在序列数量约 2 千万的样品中，有 20 个计数，
# 或者在序列数量约 7.6 千万有 76 个计数。

# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]

# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)

# 归一化，TMM 方法
d <- calcNormFactors(d)

# 注意：归一化并不会直接在counts数值上修改，而是归一化系数会被自动存在 d$samples$norm.factors

# 顺便介绍一下归一化的意义
# 归一化不是绝对必要的，但是推荐进行归一化。
# 有重复的样本中，应该不具备生物学意义的外部因素会影响单个样品的表达
# 例如中第一批制备的样品会总体上表达高于第二批制备的样品，
# 假设所有样品表达值的范围和分布都应当相似，
# 需要进行归一化来确保整个实验中每个样本的表达分布都相似。

# 将归一化后的数据赋值给 dge 变量
dge = d

# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)

# edgeR 涉及到差异表达分析的函数有很多：exactTest、glmFit、glmLRT、glmQLFit、glmQLFTest。 
# qCML估计离散度需要搭配 exact test 进行差异表达分析，对应 exactTest 函数。 
# 而其他四个glm都是与GLM模型搭配使用的函数。其中，glmFit 和 glmLRT 函数是配对使用的，
# 用于 likelihood ratio test (似然比检验)，而 glmQLFit和 glmQLFTest则配对用于 quasi-likelihood F test (拟极大似然F检验)。

# 使用 LRT（Likelihood Ratio Test）计算差异表达
# 注意这里的 contrast 和 DESeq2 不一样，这里我们只需要输入 c(-1, 1) 即可
# -1 对应 normal，1 对应 tumor
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG <- topTags(lrt, n = nrow(dge))

# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_edgeR, file = './data/DEG_edgeR.Rdata')

# 画图图 —— 火山图和热图

load("./data/DEG_edgeR.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)
# 火山图
p <- ggplot(data = DEG_edgeR, 
            aes(x = logFC, 
                y = -log10(PValue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_edgeR.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
deg_opt <- DEG_edgeR %>% filter(DEG_edgeR$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_edgeR.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

########################### 使用 limma 进行差异分析 ############################

library(limma)

# limma 包建议使用对数转换后的数据，其实我们从 Xena 下载的直接就是log2(count+1)后的数据，
# 但是这里为了展示 limma 包内部的标准化方法的使用，我们这里还是使用原始计数矩阵作为输入数据。
exprSet <- exp_brca_int

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)

# 创建 DGEList 对象
dge <- DGEList(counts = exprSet)

# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析

# 使用线性模型进行拟合
fit <- lmFit(v, design)

# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")

# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_limma_voom, file = './data/DEG_limma_voom.Rdata')

# 画图图 —— 火山图和热图

load("./data/DEG_limma_voom.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))

# 火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_limma_voom.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
deg_opt <- DEG_limma_voom %>% filter(DEG_limma_voom$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_limma_voom.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

# References: https://mp.weixin.qq.com/s/0DVzPXrNCQgG7kTtMw74jQ
# https://github.com/mikelove/DESeq2
# https://github.com/jianjinxu/edgeR
# https://github.com/gangwug/limma
# 麻了 忘了保存了 我又重写了一遍 哭死了

#########################3
# 先分析
## Bulk-RNAseq 基因表达矩阵
library(readxl)
library(dplyr)
dat <- read_excel("data_test.xlsx") 
dat <- dat %>% data.frame()
row.names(dat) <- dat$gene
dat <- dat[,-1]
head(dat)

#BiocManager::install("GSEABase")
BiocManager::install("GSVA", version = "3.14") # R 4.1.2 注意版本w
library('GSEABase')
library(GSVA)
geneSets <- getGmt('h.all.v7.5.1.symbols.gmt')    ###下载的基因集
GSVA_hall <- gsva(expr=as.matrix(dat), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)

## 差异基因分析
## limma
#BiocManager::install('limma')
library(limma)
# 设置或导入分组
group <- factor(c(rep("Tumor", 3), rep("Normal", 4)), levels = c('Tumor', 'Normal'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design
# Tunor VS Normal
compare <- makeContrasts(Tumor - Normal, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

##############################33
# 发散条形图绘制
## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggtheme)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签
# 此处参考了：https://mp.weixin.qq.com/s/eCMwWCnjTyQvNX2wNaDYXg
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
ggsave("gsva_bar.pdf",p,width = 8,height  = 8)

# References：https://blog.csdn.net/weixin_45822007/article/details/123013713
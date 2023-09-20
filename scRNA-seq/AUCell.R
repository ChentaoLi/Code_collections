# 打分
#########################
# AddModuleScore
# immune seruatobj
DefaultAssay(immune) <- "RNA"
cd_features <- list(c(
  'TNF',
  'CCL2',
  'CCL3',
  'CCL4',
  'CXCL10',
  'S100A8',
  'CXCL1'
))

Inscore <- AddModuleScore(immune,
                          features = cd_features,
                          ctrl = 100,
                          name = "CD_Features")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[27] <- 'Inflammatory_Score' 
# 其实构建了一个Inscore的seurat对象，对其进行可视化。
VlnPlot(Inscore,features = 'Inflammatory_Score', 
        pt.size = 0, adjust = 2,group.by = "orig.ident")
# UMAP
library(ggplot2)
mydata<- FetchData(Inscore,vars = c("UMAP_1","UMAP_2","Inflammatory_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Inflammatory_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
# 箱线图展示不同组织之间的评分差异
data<- FetchData(Inscore,vars = c("group","Inflammatory_Score"))
ggplot(data, aes(x=group,y=`Inflammatory_Score`)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=10,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL,title = "Regulation of necroptotic process")+ geom_jitter(col="#00000033", pch=19,cex=2.5, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(color = factor(group)))+
  NoLegend()+theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行,标题居中

#########################
# AUCell
# BiocManager::install("AUCell",force = TRUE)
library(AUCell)
library(clusterProfiler)
# 对细胞表达矩阵排列，下载GSEA基因集文件，网址：http://www.gsea-msigdb.org/gsea/downloads.jsp，选择自己需要关注的板块。进行评分计算。
cells_rankings <- AUCell_buildRankings(immune@assays$RNA@data) 
Hallmarker <- read.gmt("h.all.v7.5.1.symbols.gmt") 
geneSets <- lapply(unique(Hallmarker$term), function(x){print(x);Hallmarker$gene[Hallmarker$term == x]})
names(geneSets) <- unique(Hallmarker$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

# 选定某一个需要关注的通路，进行可视化。
##set gene set of interest here for plotting
geneSet <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
immune$AUC  <- aucs

library(ggraph)
ggplot(data.frame(immune@meta.data, immune@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "TNFA_SIGNALING_VIA_NFKB")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))


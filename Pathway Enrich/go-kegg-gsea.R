options(clusterProfiler.download.method = "wget")

library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
#install.packages('R.utils')
#require(devtools)
#install_version("vctrs", version = "0.4.1",
#                repos = "http://cran.us.r-project.org")

#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
# info <- read.xlsx( "gene.xlsx", rowNames = F,colNames = T)
CPXV008 <- read.xlsx("D:/OneDrive - International Campus, Zhejiang University/桌面/显著性差异--红色-1.xlsx", sheet = 1)
CPXV011 <- read.xlsx("D:/OneDrive - International Campus, Zhejiang University/桌面/显著性差异--红色-1.xlsx", sheet = 2)
CPXV198 <- read.xlsx("D:/OneDrive - International Campus, Zhejiang University/桌面/显著性差异--红色-1.xlsx", sheet = 3)

CPXV008$Log2FoldChange <- log2(CPXV008$`LZ022/LZ032`) # 记得更管
CPXV008$gene_symbol <- CPXV008$Gene.Name
info <- data.frame(gene_symbol = CPXV008$gene_symbol, Log2FoldChange = CPXV008$Log2FoldChange)

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

# go
GO<-enrichGO(gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 1,
              qvalueCutoff = 1)

R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)

names(info) <- c('SYMBOL','Log2FoldChange')#,'pvalue','padj'
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)#GSEA富集分析

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 20,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)

name_GO <- paste0("1.GO.png")
ggsave(filename = name_GO, barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free"), dpi = 300)
name_KEGG <- paste0("2.KEGG.png")
ggsave(filename = name_KEGG, barplot(KEGG,showCategory = 20,title = 'KEGG Pathway'), dpi = 300)

enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE

enrichplot::heatplot(GO,showCategory = 20)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 20)

GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")

write.table(KEGG$ID, file = "F:\\GO_KEGG_GSEA分析\\kegg_txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
browseKEGG(KEGG,"hsa05166")#选择其中的hsa05166通路进行展示



install
GO_CC<-enrichGO( gene$ENTREZID,#GO富集分析CC模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_CC)#GO-CC功能网络图
GO_MF<-enrichGO( gene$ENTREZID,#GO富集分析MF模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_MF)#GO-MF功能网络图

genedata<-data.frame(ID=info$gene_symbol,logFC=info$log2FoldChange)
write.table(GO$ONTOLOGY, file = "F:\\GO_KEGG_GSEA分析\\go.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[2103:2112,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[2410:2419,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)


GOCircle(circ_BP) #弦表图
GOCircle(circ_CC) 
GOCircle(circ_MF) 


chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 


ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,1:30)#30是根据ridgeplot中有30个富集通路得到的
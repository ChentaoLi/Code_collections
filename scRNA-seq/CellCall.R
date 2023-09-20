# library(devtools)
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("ShellyCoder/cellcall")
library(cellcall)

# 将GM和BM分开做互作，可以看看不同状态下细胞互作之间的区别
GM_immune <- subset(immune, group=="GM")
test <- CreateObject_fromSeurat(Seurat.object= GM_immune, #seurat对象
                                slot="counts", 
                                cell_type="celltype", #细胞类型
                                data_source="UMI",
                                scale.factor = 10^6, 
                                Org = "Homo sapiens") #物种信息
mt <- TransCommuProfile(object = test,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="mean",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')


#有多少细胞类型就设置多少个颜色
cell_color <- data.frame(color=c("#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500"), stringsAsFactors = FALSE)
rownames(cell_color) <- c("Macrophage","T cell","mDC","Neutrophil","Mast")
#绘制互作图
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.5, #细胞类型多的话设置小点，不然图太大画不出来
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)


#可视化互作受配体关系
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "45",  
             main="score")

# 转录因子分析
# 比如我想关注pDC细胞的转录因子:
pDC.tf <- names(mt@data$gsea.list$pDC@geneSets)
pDC.tf
# Draw the TF enrichment plot:
getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID = c("Rb1","Runx2","Tcf7l2"), 
            myCelltype="pDC", fc.list=mt@data$fc.list,
            selectedGeneID = mt@data$gsea.list$pDC@geneSets$Rb1[1:10],
            mycol = NULL)

# References: https://github.com/ShellyCoder/cellcall
# https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab638/6332819
# https://mp.weixin.qq.com/s?__biz=Mzg5OTYzMzY5Ng==&mid=2247484336&idx=1&sn=a63656793e9e744be901ff44e41e871a&chksm=c05104fff7268de9487ec4fd0b522666ee56bafb4804b30ccee7c134afcb80dbdc6a704d7ac9&scene=21#wechat_redirect
# https://zhuanlan.zhihu.com/p/445131235?utm_id=0
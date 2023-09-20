# GSVA Pathway Enrichment
library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(msigdbr) # msigdbr_species() #可以看到，这个包涵盖了20个物种

# 如果需要提取亚群，整体进行GSVA很费时间
immuneT <- subset(scRNA, celltype=="T cells")#提取我们需要分析的细胞类型
immuneT <- as.matrix(immuneT@assays$RNA@counts)#提取count矩阵
meta <- immuneT@meta.data[,c("orig.ident", "sex", "age", "stim", "samples")]#分组信息，为了后续作图

mouse <- msigdbr(species = "Mus musculus")
table(mouse$gs_cat) #查看目录，与MSigDB一样，包含9个数据集

mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", #GO在C5
                      subcategory = "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol，也可以选择ID，根据自己数据需求来，主要为了方便
mouse_GO_bp_Set = mouse_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)#后续gsva要求是list，所以将他转化为list

gsva.result <- gsva(expr = immuneT, 
               gset.idx.list = mouse_GO_bp_Set,
               kcdf="Poisson", #查看帮助函数选择合适的kcdf方法 
               parallel.sz = 5)
write.table(gsva.result, 'gsva.result.xls', row.names=T, col.names=NA, sep="\t")

es <- data.frame(t(gsva.result),stringsAsFactors=F)
scRNA <- AddMetaData(immuneT, es)
saveRDS(scRNA, file = "scRNA.rds")

# 接着差异分析可以用limma包，类似于转录组芯片数据分析流程。
group <- c(rep("control", 50), rep("test", 71)) %>% as.factor()#设置分组，对照在前
design <- model.matrix(~ 0 + group) #构建比较矩阵
colnames(design) <- levels(group)
fit = lmFit(es, desigN)
fit2 <- eBayes(fit)
diff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
write.csv(diff, file = "Diff.csv")

# 最后对差异的感兴趣的通路进行可视化：
up <- c("GOBP_EGG_ACTIVATION",
        "GOBP_TENDON_DEVELOPMENT",
        "GOBP_SOMITE_SPECIFICATION",
        "GOBP_THREONINE_CATABOLIC_PROCESS",
        "GOBP_REGULATION_OF_GLUTAMATE_RECEPTOR_CLUSTERING",
        "GOBP_NEGATIVE_CHEMOTAXIS",
        "GOBP_NEGATIVE_REGULATION_OF_FAT_CELL_PROLIFERATION",
        "GOBP_REGULATION_OF_T_HELPER_17_CELL_LINEAGE_COMMITMENT",
        "GOBP_REGULATION_OF_ANTIMICROBIAL_HUMORAL_RESPONSE")
down <- c("GOBP_DETERMINATION_OF_PANCREATIC_LEFT_RIGHT_ASYMMETRY",
          "GOBP_MITOTIC_DNA_REPLICATION",
          "GOBP_EOSINOPHIL_CHEMOTAXIS",
          "GOBP_NEUTROPHIL_MEDIATED_CYTOTOXICITY",
          "GOBP_POTASSIUM_ION_EXPORT_ACROSS_PLASMA_MEMBRANE",
          "GOBP_REGULATION_OF_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
          "GOBP_REGULATION_OF_SEQUESTERING_OF_ZINC_ION",
          "GOBP_ENDOTHELIN_RECEPTOR_SIGNALING_PATHWAY",
          "GOBP_PRE_REPLICATIVE_COMPLEX_ASSEMBLY_INVOLVED_IN_CELL_CYCLE_DNA_REPLICATION",
          "GOBP_ESTABLISHMENT_OF_PLANAR_POLARITY_OF_EMBRYONIC_EPITHELIUM")
TEST <- c(up,down)
diff$ID <- rownames(diff) 
Q <- diff[TEST,]
group1 <- c(rep("treat", 9), rep("control", 10)) 
df <- data.frame(ID = Q$ID, score = Q$t,group=group1 )
# 按照t score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列
ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="t value of GSVA score")+
  scale_fill_manual(values = c("#008020","#08519C"))#设置颜色

# References: https://www.bioconductor.org/packages/release/bioc/html/GSVA.html


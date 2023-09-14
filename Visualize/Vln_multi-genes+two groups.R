HCD4T <- "CD4T_human_ref_v1.rds"
HCD8T <- "CD8T_human_ref_v1.rds"
MAPC <- "APC_atlas_v1_SPICA.rds"
HDC <- "DC_human_ref_v1.rds"
MTIL <- "ref_TILAtlas_mouse_v1.rds"
MCD8T_LCMV <- "ref_CD8_LCMV_mouse_v2.rds"
MCD4T_LCMV <- "ref_LCMV_CD4_mouse_v1.rds"
MCD8T_LCMV_v1 <- "ref_CD8_LCMV_mouse_v1.rds"

Vln <- function(data, gene_list, group, split){
  # source("../Code_collections/Visualize/singlecell_gene_test.R")
  source("D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/Visualize/singlecell_gene_test.R")
  A <- singlecell_gene_test(data, genes.use = gene_list, 
                            group.by = group, comp = split)
  anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
  anno_sig <- A$sig
  plots_violins <- VlnPlot(data, 
                           cols = c("limegreen", "navy"),
                           pt.size = 0,
                           group.by = "orig.ident",
                           features = gene_list, 
                           ncol = 3, 
                           log = FALSE,
                           combine = FALSE)+
    geom_boxplot(width=.2, col="black",fill="white")
  for(i in 1:length(plots_violins)) {
    data <- plots_violins[[i]]$data
    colnames(data)[1] <- 'gene'
    plots_violins[[i]] <- plots_violins[[i]] + 
      theme_classic() + 
      theme(axis.text.x = element_text(size = 10,color="black"),
            axis.text.y = element_text(size = 10,color="black"),
            axis.title.y= element_text(size=12,color="black"),
            axis.title.x = element_blank(),
            legend.position='none')+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
      scale_x_discrete(labels = c("Female","Male"))+
      geom_signif(annotations = anno_sig[i],
                  y_position = max(data$gene)+0.5,
                  xmin = 1,
                  xmax = 2,
                  tip_length = 0)
  }
  CombinePlots(plots_violins)
}

Vln(SerautObj, c("Ripk3", "Gzmb", "Tox"), "functional.cluster", "GeneName_label")
SerautObj <- PROJECTIL
genes.use <- c("Ripk3", "Gzmb", "Tox")
group.by<- "functional.cluster"
comp <- "GeneName_label"
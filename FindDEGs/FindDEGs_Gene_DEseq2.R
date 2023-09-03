# ProjecTIL
HCD4T <- "CD4T_human_ref_v1.rds"
HCD8T <- "CD8T_human_ref_v1.rds"
MAPC <- "APC_atlas_v1_SPICA.rds"
HDC <- "DC_human_ref_v1.rds"
MTIL <- "ref_TILAtlas_mouse_v1.rds"
MCD8T_LCMV <- "ref_CD8_LCMV_mouse_v2.rds"
MCD4T_LCMV <- "ref_LCMV_CD4_mouse_v1.rds"
MCD8T_LCMV_v1 <- "ref_CD8_LCMV_mouse_v1.rds"

FindDEGs <- function(GeneName, ProjecTIL_path){
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  PROJECTIL <- readRDS(ProjecTIL_path)
  PROJECTIL_expression <- PROJECTIL@assays$RNA[GeneName, ]
  PROJECTIL_GeneName_expression <- as.vector(PROJECTIL@assays$RNA[GeneName, ])
  PROJECTIL@meta.data$GeneName_expression <- PROJECTIL_GeneName_expression
  label_high <- paste0(GeneName, " High")
  label_low <- paste0(GeneName, " Low")
  PROJECTIL@meta.data$GeneName_label <- ifelse(PROJECTIL@meta.data$GeneName_expression > 0, label_high, label_low)
  PROJECTIL_subsets <- SplitObject(PROJECTIL, split.by = "functional.cluster")
  PROJECTIL_subsets_with_GeneName <- lapply(PROJECTIL_subsets, function(subset) {
  SplitObject(subset, split.by = "GeneName_label")
  })
  library("DESeq2")
  results_list <- list()
  for (subset_name in names(PROJECTIL_subsets_with_GeneName)) {
    subset_data <- PROJECTIL_subsets_with_GeneName[[subset_name]]
    if (is.null(subset_data[[label_high]])) {
      print(paste("Subset", subset_name, "is NULL"))
    } else {
      subset_expression <- cbind(as.matrix(subset_data[[label_low]]@assays[["RNA"]]@data),as.matrix(subset_data[[label_high]]@assays[["RNA"]]@data))
      subset_meta <- rbind(subset_data[[label_low]]@meta.data,subset_data[[label_high]]@meta.data)
      dds <- DESeqDataSetFromMatrix(countData = round(2^subset_expression),
                                    colData = subset_meta,
                                   design = ~ GeneName_label)
      dds <- DESeq(dds)
      results <- results(dds)
      results_list[[subset_name]] <- results
   }}
  library(openxlsx)
  wb <- createWorkbook()
  df_names <- names(results_list)
  for (i in seq_along(results_list)) {
   df <- results_list[i][[1]]
   df$Gene <- rownames(df)
   df <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue) & abs(df$log2FoldChange) > 0.5 & df$pvalue < 0.05, ]
   sheet_name <- df_names[i]
   addWorksheet(wb, sheetName = sheet_name)
   df_new = as.data.frame(df)
   writeDataTable(wb, sheet = i, x = df_new, colNames = TRUE, rowNames = TRUE)
  }
  a <- gsub(pattern = "D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/ProjectTIL/", replacement = "", x = ProjecTIL_path)
  b <- gsub(pattern = ".rds", replacement = "", x = a)
  filename <- paste0(b, "_DEseq2_", GeneName, ".xlsx")
  saveWorkbook(wb, filename)

  library(readxl)
  library(UpSetR)
  data_frames <- list()
  gene_names_list <- list()
  sheet_names <- excel_sheets(filename)
  for (sheet_name in sheet_names) {
    sheet_data <- read_excel(filename, sheet = sheet_name)
    data_frames[[sheet_name]] <- sheet_data
    gene_names_list[[sheet_name]] <- sheet_data$Gene
  }
  
  max_length <- max(sapply(gene_names_list, length))
  filled_gene_names_list <- lapply(gene_names_list, function(x) {
    if (length(x) < max_length) {
     c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })
  gene_names_df <- as.data.frame(filled_gene_names_list)
  colnames(gene_names_df) <- sheet_names
  unique_genes <- unique(unlist(gene_names_df))
  gene_matrix <- matrix(0, nrow = length(unique_genes), ncol = ncol(gene_names_df))
  rownames(gene_matrix) <- unique_genes
  colnames(gene_matrix) <- colnames(gene_names_df)
  for (i in 1:nrow(gene_matrix)) {
   for (j in 1:ncol(gene_matrix)) {
      if (rownames(gene_matrix)[i] %in% gene_names_df[, j]) {
        gene_matrix[i, j] <- 1
     }
    }
  }
  gene_matrix <- gene_matrix[rowSums(gene_matrix) > 2, ]
  a <- gsub(pattern = ".xlsx", replacement = "", x = filename)
  filename2 <- paste0(a, "_UpsetR", ".csv")
  write.csv(gene_matrix, filename2)
  gene_matrix <- as.data.frame(gene_matrix)
  fig_name <- paste0(a, "_UpsetR", ".png")
  upset_data <- upset(gene_matrix, nset = 20, nintersects = 100, 
                      mainbar.y.label = "Intersection Size", sets.x.label = "DEG counts", 
                      mb.ratio = c(0.4, 0.6),
                      set_size.show = TRUE)
  
  upset_data
  png(file = fig_name, res = 300, width = 1200, height = 960)
  print(upset_data)
  dev.off()
  gene_list <- colnames(gene_matrix)

  library(scCustomize) 
  library(dittoSeq)
  library(ggplot2)
  fig_name <- paste0(a, "_vln", ".png")
  png(file = fig_name, res = 300, width = 1200, height = 960)
  VlnPlot(PROJECTIL, features = gene_list,
              stack=T,pt.size=0,
              flip = T,
              add.noise = T,
              split.by = 'GeneName_label',
              split.plot = T)+#横纵轴不标记任何东西
  theme(axis.text.y = element_blank(), #不显示坐标刻度
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 45),
        legend.position = 'top')
  dev.off()
}

# FindDEGs("Ripk3", MCD8T_LCMV)


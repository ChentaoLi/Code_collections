library(openxlsx)
scRNA <- readRDS("data.rds")
scRNA$CellType <- Idents(scRNA)
# 定义CellType和method的组合 (1st cluster name and 2nd cluster name)
cell_type_values <- unique(scRNA@meta.data$CellType)

# 在每个亚群里面比较不同组
for (cell_type in cell_type_values) {
  cell_selection <- subset(scRNA, cells = colnames(scRNA)[scRNA$CellType == cell_type])
  Idents(cell_selection) <- cell_selection$method
  DGE_cell_selection <- FindAllMarkers(cell_selection, logfc.threshold = 0.25)
  
  # save DEG
  result_excel <- createWorkbook()
  clusters <- unique(DGE_cell_selection$cluster)
  for (cluster_val in clusters) {
    subset_data <- DGE_cell_selection[DGE_cell_selection$cluster == cluster_val, ]
    addWorksheet(result_excel, sheetName = as.character(cluster_val))
    writeData(result_excel, sheet = as.character(cluster_val), x = subset_data)
  }
  result_dir <- "DEG_results"
  result_filename <- paste0(result_dir, "/", cell_type, ".xlsx")
  saveWorkbook(result_excel, result_filename, overwrite = TRUE)
}

# 在每个亚群里面比较指定两个组
result_excel <- createWorkbook()
for (cell_type in cell_type_values) {
  cell_selection <- subset(scRNA, cells = colnames(scRNA)[scRNA$CellType == cell_type])
  Idents(cell_selection) <- cell_selection$method
  DGE_cell_selection <- FindMarkers(cell_selection, ident.1 = "PD1-19bbz", ident.2 = "LV-19bbz_PD1-KO" , logfc.threshold = 0.25)
  addWorksheet(result_excel, sheetName = cell_type)
  writeData(result_excel, sheet = cell_type, x = DGE_cell_selection)
}
result_filename <- paste0("DEG_results/", "DEGs", "_LV-ET_DGE_results.xlsx")
saveWorkbook(result_excel, result_filename, overwrite = TRUE)
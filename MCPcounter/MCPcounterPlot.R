MCPcounterPlot <- function(genename){
  load("D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/MCPcounter/MCPcounter_exp.rdata")
  load("D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/MCPcounter/MCPcounter.rdata")
  a <- exp[genename, ]
  results <- as.data.frame(results)
  MCPcounter_matrix <- as.matrix(rbind(results, a))
  sorted_matrix <- MCPcounter_matrix[, order(MCPcounter_matrix[genename, ], decreasing = TRUE)]
  # results <- t(scale(t(sorted_matrix[-nrow(sorted_matrix),])))
  results <- sorted_matrix[-nrow(sorted_matrix),]
  Gene_exp <- as.data.frame(sorted_matrix[genename,])
  colnames(Gene_exp) <- genename

  library(pheatmap)
  color_breaks <- c(seq(-3, 3, length.out = 101))
  ann_colors <- list(genename = c("red","white",'blue'))
  filename <- paste0(genename, "_MCPcounterPlot.png")
  # results <- t(scale(t(results)))
  pheatmap(results,
           breaks = color_breaks,
           scale = "row", 
           cluster_col = FALSE,
           cutree_rows=4,
           annotation_col = Gene_exp,
           annotation_colors = ann_colors,
           show_colnames = FALSE, 
           border_color = "black",
           filename = filename,
           cellwidth = 0.01, 
           cellheight = 40, 
           )

}

# MCPcounterPlot("GSDMD")
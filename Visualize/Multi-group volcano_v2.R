readAllSheets <- function(xlsx_path) {
  library(openxlsx)
  library(readxl)
  all_sheets <- excel_sheets(xlsx_path)
  data_frames <- list()
  sheet_names <- list()
  
  for (sheet_name in all_sheets) {
    df <- read.xlsx(xlsx_path, sheet = sheet_name)
    data_frames[[sheet_name]] <- df
    sheet_names[[sheet_name]] <- sheet_name
  }
  
  return(list(data_frames = data_frames, sheet_names = sheet_names))
}
createMultiVolcanoPlot <- function(result, output_pdf_path) {
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(openxlsx)
  library(tidyr)
  library(dplyr)
  cluster_names <- result$sheet_names
  data_frames <- result$data_frames
  # Get the maximum and minimum values from each data frame in data_frames
  max_values <- sapply(data_frames, function(df) max(df$log2FoldChange, na.rm = TRUE))
  min_values <- sapply(data_frames, function(df) min(df$log2FoldChange, na.rm = TRUE))
  
  data_frames <- Map(function(df, cluster) {
    df %>%
      mutate(cluster = cluster)
  }, data_frames, cluster = cluster_names)
  data_frames <- Map(function(df, group) {
    df %>%
      mutate(group = group)
  }, data_frames, group = 1:length(data_frames))
  
  get_top_genes <- function(df, criterion = "log2FoldChange", n = 10) {
    return(df[order(-abs(df[[criterion]]))[1:n], ])
  }
  
  all_top_genes <- do.call(rbind, lapply(data_frames, get_top_genes))
  
  for (i in seq_along(data_frames)) {
    data_frames[[i]]$size <- ifelse(data_frames[[i]]$Gene %in% all_top_genes$Gene, 2, 1)
  }
  
  all_data <- do.call(rbind, data_frames)
  dt <- filter(all_data, size == 1)
  
  p <- ggplot() +
    geom_jitter(data = dt,
                aes(x = group, y = log2FoldChange, color = change),
                size = 0.4,
                width = 0.4) +
    geom_jitter(data = all_top_genes,
                aes(x = group, y = log2FoldChange, color = change),
                size = 1,
                width = 0.4)
  
  # Set the values for dfbar and dfbar1 based on max_values and min_values
  dfbar <- data.frame(x = 1:length(max_values), y = max_values)
  dfbar1 <- data.frame(x = 1:length(min_values), y = min_values)
  
  p1 <- ggplot() +
    geom_col(data = dfbar,
             mapping = aes(x = x, y = y),
             fill = "#dcdcdc", alpha = 0.6) +
    geom_col(data = dfbar1,
             mapping = aes(x = x, y = y),
             fill = "#dcdcdc", alpha = 0.6)
  
  p2 <- p1 +
    geom_jitter(data = dt,
                aes(x = group, y = log2FoldChange, color = change),
                size = 0.4,
                width = 0.4) +
    geom_jitter(data = all_top_genes,
                aes(x = group, y = log2FoldChange, color = change),
                size = 1,
                width = 0.4)
  num_sheets <- length(data_frames)
  dfcol <- data.frame(x = 1:(num_sheets), y = 0, label = 1:(num_sheets))
  library(viridis)
  mycol <- viridis_pal(option = "D")(nrow(dfcol))
  average_log2FoldChange <- mean(dt$log2FoldChange) + 1
  p3 <- p2 + geom_tile(data = dfcol,
                       aes(x = x, y = y, fill = mycol),
                       height = average_log2FoldChange, # 调整中间色块的大小
                       color = "black",
                       alpha = 0.9,
                       show.legend = FALSE)
  p4 <- p3 +
    geom_text_repel(
      data = all_top_genes,
      aes(x = group, y = log2FoldChange, label = Gene),
      force = 0.8,
      min.segment.length = 0.5,
      segment.color = "black",
      arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last"),
      max.overlaps = 200
    )
  
  p5 <- p4 +
    scale_color_manual(name = NULL,
                       values = c("black", "#91D1C27F", "red"))
  
  p6 <- p5 +
    labs(x = "Volcano Plots of DEGs in different clusters", y = "log2FoldChange") +
    geom_text(data = dfcol,
              aes(x = x, y = y, label = cluster_names),
              size = 4,
              color = "black")
  
  p7 <- p6 +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 13, color = "black", face = "bold"),
      axis.line.y = element_line(color = "black", size = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1, 0),
      legend.text = element_text(size = 15)
    )
  
  # Save the plot as a PDF
  pdf(output_pdf_path, width = 14, height = 10)
  print(p7)
  dev.off()
}

# result <- readAllSheets("DEseq2-results.xlsx")
# data_frames <- result$data_frames
# sheet_names <- result$sheet_names
# createMultiVolcanoPlot(result, "multi_volcano.pdf")

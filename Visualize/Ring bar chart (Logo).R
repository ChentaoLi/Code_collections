library(tidyverse)
library(countrycode)
library(ggflags)
library(ggtext)

# Data
wins <- read_csv("data.csv")
# Data washing
wins_by_cat <- wins %>% 
  group_by(Nationality, Category) %>%  # 按国籍和类别分组
  summarise(cat_total = n()) %>%  # 计算每组的总数
  group_by(Category) %>%  # 再按类别分组
  mutate(rank = rank(-cat_total, ties.method = "min")) %>%  # 对每组进行排名
  arrange(Category, desc(cat_total)) %>%  # 排序
  filter(rank <= 10) %>%  # 只保留前10名
  ungroup() %>%  # 取消分组
  mutate(Nationality = case_when(  # 修改国籍名称
    Nationality == "United Kingdom" ~ "UK",
    Nationality == "United States" ~ "USA",
    Nationality == "Soviet Union" ~ "Russia",
    TRUE ~ Nationality
  )) %>% 
  mutate(code = countrycode(Nationality, "country.name", "iso2c") %>% tolower())  # 转换为国家代码

# Add emepty rows for splicing the figure
empty_nrow = 3
to_add <- data.frame(matrix(NA, nrow = empty_nrow * n_distinct(wins_by_cat$Category), ncol = ncol(wins_by_cat)))
colnames(to_add) <- colnames(wins_by_cat)
to_add$Category <- rep(unique(wins_by_cat$Category), each = empty_nrow)

# 添加空行到数据集
wins_by_cat_space <- wins_by_cat %>% 
  rbind(to_add) %>% 
  arrange(Category, desc(cat_total)) %>% 
  mutate(seq_id = row_number())

# 计算坐标标签
label_data <- wins_by_cat_space
nBar <- nrow(label_data)
angle <- 90 - 360 * (label_data$seq_id - 0.5) / nBar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

# Set colors
bg_color <- "#13293d"
line_color <- "#e9f1f2"
font_color <- "white"
bar_pal <- c("Men" = "#679289", "Wheelchair Men" = "#c9cba3",
             "Women" = "#ee2e31", "Wheelchair Women" = "#f4c095")

# 构建极坐标坐标图
polar_barplot <- ggplot(wins_by_cat_space) +
  geom_bar(aes(x = seq_id, y = cat_total, fill = Category),
           stat = "identity", position = "dodge",
           show.legend = FALSE, alpha = 0.9, color = line_color, linewidth = 0.5) +
  coord_polar() +
  scale_y_continuous(limits = c(-8, 20)) +
  scale_fill_manual(values = bar_pal, guide = "none") +
  scale_color_manual(values = bar_pal, guide = "none") +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(panel.background = element_rect(fill = bg_color, colour = bg_color),
        plot.background = element_rect(fill = bg_color, colour = bg_color),
        text = element_text(color = line_color),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

# 构建分组条带和标签数据
base_data <- wins_by_cat_space %>% 
  group_by(Category) %>% 
  summarize(start = min(seq_id), end = max(seq_id) - empty_nrow) %>% 
  rowwise() %>%
  mutate(title = mean(c(start, end)),
         printname = case_when(
           Category == "Men" ~ "MEN",
           Category == "Women" ~ "WOMEN",
           Category == "Wheelchair Men" ~ "Wheelchair\nMEN",
           Category == "Wheelchair Women" ~ "Wheelchair\nWOMEN"
         ),
         angle = 360 * (title) / nrow(wins_by_cat_space),
         angle = ifelse(angle < 180, angle - 90, angle + 90))

# 添加分组标签和条带
polar_barplot +
  # 添加数值标签
  geom_text(data = label_data, 
            aes(x = seq_id, y = cat_total - 1.5, label = ifelse(rank <= 3, cat_total, ""), hjust = hjust), 
            color = font_color, fontface = "bold",
            alpha = 1, size = 3, angle = label_data$angle, inherit.aes = FALSE) +
  # 添加国旗
  ggflags::geom_flag(data = label_data, 
                     aes(x = seq_id, y = cat_total + 1.25, country = code), size = 3.5) +
  # 添加国家名称
  geom_text(data = label_data, 
            aes(x = seq_id, y = cat_total + 2.5, label = Nationality, hjust = hjust), 
            color = font_color, fontface = "bold",
            alpha = 0.7, size = 3, angle = label_data$angle, inherit.aes = FALSE) +
  # 添加分组条带
  geom_segment(data = base_data, 
               aes(x = start + 1, y = max(wins_by_cat_space$cat_total, na.rm = TRUE) + 0.5, 
                   xend = end, yend = max(wins_by_cat_space$cat_total, na.rm = TRUE) + 0.5,
                   color = Category), 
               alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
  # 添加分组名称
  geom_text(data = base_data, 
            aes(x = title + 0.5, y = max(wins_by_cat_space$cat_total, na.rm = TRUE) - 1.5, 
                label = printname, angle = angle, colour = Category), 
            hjust = c(1, 0.5, 0.5, 0.5), 
            alpha = 1, size = 3.25, fontface = "bold", inherit.aes = FALSE)

# References: https://mp.weixin.qq.com/s/SNfYZvLeUyLL06XnaSq5tA
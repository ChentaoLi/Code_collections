table(scRNA$orig.ident)#查看各组细胞数
prop.table(table(Idents(scRNA)))
table(Idents(scRNA), scRNA$orig.ident)#各组不同细胞群细胞数
# 计算各组样本不同细胞群比例
Cellratio = prop.table(table(Idents(scRNA), scRNA$TIME), margin = 2)#TIME should be changed!!!
Cellratio = as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = 'Sample', y = 'Ratio') +
  scale_fill_manual(values = allcolour) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
# References: https://blog.csdn.net/qq_42090739/article/details/124246219
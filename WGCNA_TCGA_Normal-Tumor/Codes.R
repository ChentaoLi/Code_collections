# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)
library(data.table)
a1=fread( 'TCGA-CHOL.htseq_counts.tsv.gz' , 
          data.table = F)
dim(a1)
a1[1:4,1:4]
a1[(nrow(a1)-5):nrow(a1),1:4]
dim(a1)
# all data is then log2(x+1) transformed.

#length(unique(a1$AccID))
#length(unique(a1$GeneName))

mat= a1[,2:ncol(a1)] 
mat[1:4,1:4] 
mat=mat[1:(nrow(a1)-4),]
mat=ceiling(2^(mat)-1) #log2(x+1) transformed.
mat[1:4,1:4] 

rownames(mat) = gsub('[.][0-9]+','',a1$Ensembl_ID[1:(nrow(a1)-4)])
keep_feature <- rowSums (mat > 1) > 1
colSums(mat)/1000000
table(keep_feature)
mat <- mat[keep_feature, ]
mat[1:4,1:4] 
mat=mat[,  colSums(mat)/1000000 >10]
dim(mat) 
colnames(mat) 
ensembl_matrix=mat
colnames(ensembl_matrix) 
ensembl_matrix[1:4,1:4]

library(AnnoProbe)
head(rownames(ensembl_matrix))
ids=annoGene(rownames(ensembl_matrix),'ENSEMBL','human')
head(ids)
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
symbol_matrix= ensembl_matrix[match(ids$ENSEMBL,
                                    rownames(ensembl_matrix)),]

rownames(symbol_matrix) = ids$SYMBOL
#symbol_matrix = ensembl_matrix
symbol_matrix[1:4,1:4]

# write.table( colnames(symbol_matrix),file = 'group.txt',
#              quote = F,
#              row.names = F,col.names = F)
library(stringr)
symbol_matrix[1:4,1:4]
colnames(symbol_matrix)
# group_list=ifelse( grepl('PLVX',colnames(symbol_matrix)),'control','case' )
group_list=ifelse(substring(colnames(symbol_matrix),14,15)=='11',
                  'control','case' )
table(group_list)
group_list = factor(group_list,levels = c('control','case' ))
group_list
# save(symbol_matrix, group_list,
#      file='symbol_matrix.Rdata')
# load(  file='symbol_matrix.Rdata')
colnames(symbol_matrix)

save(symbol_matrix,group_list,file = 'symbol_matrix.Rdata')

source('scripts/step2-qc-counts.R')
source('scripts/step3-deg-deseq2.R')
source('scripts/step3-deg-edgeR.R')
source('scripts/step3-deg-limma-voom.R') 
source('scripts/step4-qc-for-deg.R') 
source('scripts/step5-anno-by-GSEA.R') 
source('scripts/step5-anno-by-ORA.R')  









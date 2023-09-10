# install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
# library(devtools)
# install_github("ebecht/MCPcounter",ref="master", subdir="Source")
# library(easyTCGA)
# tcga_expr_file <- "D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/download_from_xena/tcga_RSEM_gene_tpm.gz"
# tcga_clin_file <- "D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/download_from_xena/Survival_SupplementalTable_S1_20171025_xena_sp.gz"
# getpancancer_xena(tcga_expr_file = tcga_expr_file, tcga_clin_file = tcga_clin_file, type = "tcga")

library(curl)
library(MCPcounter)
library(org.Hs.eg.db)
library(stringi)
library(clusterProfiler)
library(tinyarray)
load("D:/OneDrive - International Campus, Zhejiang University/Dry_Lab/Code_collections/TCGA/TCGA_pancancer_expr.rdata")
Ensembl_ID <- tcga_expr$gene_id
Ensembl_ID <- gsub("\\..*", "",  Ensembl_ID)
tcga_expr$gene_id <- Ensembl_ID
rownames(tcga_expr) <- tcga_expr$gene_id
tcga_expr$gene_id <- NULL
exp = trans_exp(tcga_expr)

probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
featuresType <- c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[2]
results <- MCPcounter.estimate(exp,featuresType=featuresType,
                              probesets=probesets,
                              genes=genes)
write.csv(results,'MCPcounter.csv')
save(results,file = 'MCPcounter.rdata')

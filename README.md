# A shared code collection by Chentao Li

> **Note**
> Hello and welcome to this website. It contains some codes I wrote or adapted from others for transcriptome/proteome data analysis. 
> Please feel free to contact me if you have any suggestions or need help!

## Content
### Data format conversion
The code in here is provided for converting different single cell data formats now.

### Data_TCGA-Pancancer
The code here is provided for organizing TCGA data for easy access and analysis. Some of the publicly available analysis methods and websites are collected in "BioInfor_Tools_Datasets_Collections.xlsx", please feel free to choose the ones you need.

### DEGs
The code in here is provided for finging DEGs based on DEseq2. edgeR and limma.

### FindDEGs-scRNA
This has code written in it based on ProjecTILs atlas. This code has been built in as a function. Classify each subpopulation in the atlas into high and low expression subpopulations according to the expression of a gene and perform differential gene analysis (DEseq2 and limma). It also outputs gene sets and automatically finds intersections and plots Upset plots.

### Visualize
The code in here is provided for visualizing.

#### Multi-group volcano.R and V2
The code can be sourced to plot DEGs in different groups. The V1 code is only for plotting 8 up-regulated genes, and the V2 is for 10 most significant changed ones. 



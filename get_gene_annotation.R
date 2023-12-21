if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

install.packages('clusterProfiler')
library(clusterProfiler)

# function to get the gene symbol, entrezid and gene name
# input: data frame
# output: data frame with gene annotation information
library(AnnotationDbi)
install.packages('org.Hs.eg.db')
library(org.Hs.eg.db)

get_gene_annotation <- function(df){
  
  ids <- bitr(df, fromType = "ENSEMBL",
              toType = c("SYMBOL", "ENTREZID", "GENENAME"),
              OrgDb = 'org.Hs.eg.db')
}
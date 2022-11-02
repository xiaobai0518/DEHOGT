
####load the dependent packages in R
library(MASS)
library(DESeq2)
library(DESeqAnalysis)
library(ROSE)
library(edgeR)
library(multcomp)
library(MLmetrics)
library(VennDiagram)
library(pROC)
library(progress)
library(dplyr)
library(tidyverse)
library(ggplot2)



# Input: Read count data, treatment 
# Output: identified DE genes
# Implement DESeq method
deseq_func <- function(data, treatment, padj = TRUE, pval_thre=0.05, l2fc=FALSE, l2fc_thre=NULL){
  
  ### input:
  ### data: N (total number of genes) by S (total number of subjects) matrix
  ### treatment: S(total number of subjects)-length integers indicating treatment received
  ### padj: TRUE or FALSE, adjust p values by BH procedure
  ### pval_thre: float, threshold of p values for identification of DE genes
  ### l2f: TRUE or FALSE, identify DE genes with/without log2foldchange values
  ### l2f_val: float, threshold of DE genes identification for log2foldchange
  
  dds = DESeqDataSetFromMatrix(data, DataFrame(treatment), design=~treatment)
  
  dds_results = DESeq(dds)
  dds_results = results(dds_results)
  pvals = dds_results$pvalue
  if (padj){
    pvals = p.adjust(pvals, method="BH")
  }
  log2fold = dds_results$log2FoldChange
  idx = (pvals < pval_thre) * 1
  
  if (l2fc){
    
    idx = idx * (log2fold > l2fc_thre)
    
  }
  
  norm_factors = estimateSizeFactors(dds)$sizeFactor
  
  output = list(DE_idx = idx, pvals = pvals, log2fc = log2fold, norm_factors = norm_factors)
  
  return(output)
}


# Implement edgeR method
edgeR_func <- function(data, treatment, padj = TRUE, pval_thre=0.05, l2fc=FALSE, l2fc_thre=NULL){
  
  ### input:
  ### data: N (total number of genes) by S (total number of subjects) matrix
  ### treatment: S(total number of subjects)-length integers indicating treatment received
  ### padj: TRUE or FALSE, adjust p values by BH procedure
  ### pval_thre: float, threshold of p values for identification of DE genes
  ### l2f: TRUE or FALSE, identify DE genes with/without log2foldchange values
  ### l2f_val: float, threshold of DE genes identification for log2foldchange
  
  edgedata = DGEList(counts = data, group = treatment)
  edgedata = calcNormFactors(edgedata)
  edgedata = estimateCommonDisp(edgedata, verbose=T)
  edgedata = estimateTagwiseDisp(edgedata)
  edge_results = exactTest(edgedata, pair=c(1, 2))
  pvals = edge_results$table$PValue
  if (padj){
    pvals = p.adjust(pvals, method="BH")
  }
  log2fold = edge_results$table$logFC
  idx = (pvals < pval_thre) * 1
  
  if (l2fc){
    
    idx = idx * (log2fold > l2fc_thre)
  } 
  
  norm_factors = edgedata@.Data[[2]]$norm.factors
  
  output = list(DE_idx = idx, pvals = pvals, log2fc = log2fold, norm_factors = norm_factors)
  
  return(output)
}

# Implement proposed DEHOGT method

glm_func <- function(data, treatment, norm_factors=NULL, dist = "qpois", padj = TRUE, pval_thre=0.05, l2fc=FALSE, l2fc_thre=NULL){
  
  ### input:
  ### data: N (total number of genes) by S (total number of subjects) matrix
  ### treatment: S(total number of subjects)-length integers indicating treatment received
  ### padj: TRUE or FALSE, adjust p values by BH procedure
  ### pval_thre: float, threshold of p values for identification of DE genes
  ### l2f: TRUE or FALSE, identify DE genes with/without log2foldchange values
  ### l2f_val: float, threshold of DE genes identification for log2foldchange
  ### dist: working distribution for read count data: quasi-Poisson (qpois) and negative Binominal (negbin) 
  dds = DESeqDataSetFromMatrix(data, DataFrame(treatment), ~treatment)
  
  if (is.null(norm_factors)){
    norm_factors = rep(1, ncol(data))
  }
  
  
  result_GLM = matrix(0, nrow(data), 4)
  
  for(i in 1: nrow(data)){
    
    genewise_data = data[i, ]
    genewise_dataframe = data.frame(read = genewise_data, treatment=treatment)
    
    if (dist == "qpois"){
      model = glm(read ~ treatment + offset(log(norm_factors)), family = quasipoisson(), data = genewise_dataframe)
    } else if (dist == "negbin"){
      model = glm.nb(read ~ treatment + offset(log(norm_factors)), data = genewise_dataframe)
    }
    # result = summary(glht(model, linfct=matrix(c(-1, 1), 1)))
    result_GLM[i, 1: 2] = c(summary(model)$coefficients[2, 1], summary(model)$coefficients[2, 4])
    # result_GLM[i, 1: 2] = c(result$test$coefficients, result$test$pvalues[[1]])
    
  }
  
  result_GLM[, 3] = p.adjust(result_GLM[, 2], method="BH")
  result_GLM[, 4] = log(exp(abs(result_GLM[, 1])), base=2)
  if (padj){
    pvals = result_GLM[, 3]
  } else{
    pvals =result_GLM[, 2]
  }
  
  log2fold = result_GLM[, 4]
  
  idx = 1 * (pvals < pval_thre)
  
  if (l2fc){
    idx = idx * (log2fold > l2fc_thre)
  } 
  
  output = list(DE_idx = idx, pvals = pvals, log2fc = log2fold)
  
  return(output)
}
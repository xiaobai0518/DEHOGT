library(MASS)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library("edgeR")
#install.packages("MASS")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages(
  pkgs = "DESeqAnalysis",
  repos = c(
    "r.acidgenomics.com",
    BiocManager::repositories()
  )
)

#install.packages("matrixStats")
#install.packages("tictoc")
library(DESeq2)
library(DESeqAnalysis)
library(openxlsx)
library(tibble)
library(data.table)
library(dplyr)
#install.packages('Rmpfr')
library(Rmpfr)
library(multcomp)
#install.packages("pROC")
library(pROC)
library(edgeR)

rqpois <- function(n, mu, theta){
  return(rnbinom(n = n, mu = mu, size = mu/(theta - 1)))
}


#Read count generating under simulation setting 2: high discrepancy of expression levels

sim_data <- function(N, N_DE, S_1, S_2, dispersion=c(0.5, 1), dist="negbin", seed=NULL){
  
  # N: number of genes
  # N_DE: number of differential expressed genes between two subsets of samples
  # S_1: number of the first subset of samples
  # S_2: number of the second subset of samples
  
  # potential issues:
  # The variance of negative binomial is too large
  
  if (! is.null(seed)){
    set.seed(seed)
  }
  
  # Some constants 
  S = S_1 + S_2 # total number of samples
  
  phi_g = runif(N, min=dispersion[1], max=dispersion[2]) # dispersion parameter
  
  DE_gene_idx = sample(1: N, N_DE) # index of differential expressed upregulated genes
  rest_gene_idx = setdiff(1: N, DE_gene_idx)
  
  M = runif(S, 600, 800)
  mean = matrix(0, N, S)
  
  for (j in 1: S_1){
    mean[, j] = runif(N, min=1, max=M[j])
  }
  
  for (j in (S_1 + 1): S){
    for (i in 1: N){
      if (i %in% DE_gene_idx){
        mean[i, j] = rexp(1, rate=1/100) + M[j]
      } else{
        mean[i, j] = runif(1, min=1, max=M[j])
      }
    }
  }
  
  data = matrix(0, N, S)
  
  for (i in 1: N){
    if (dist == "negbin"){
      data[i, 1: S_1] = rnegbin(S_1, mu = mean[i, 1: S_1], theta = phi_g[i])
      data[i, (S_1 + 1): S] = rnegbin(S_2, mu = mean[i, (S_1 + 1): S], theta = phi_g[i])
    } else if (dist == "qpois"){
      data[i, 1: S_1] = rqpois(S_1, mu = mean[i, 1: S_1], theta = phi_g[i])
      data[i, (S_1 + 1): S] = rqpois(S_2, mu = mean[i, (S_1 + 1): S], theta = phi_g[i])
    }
    
  }
  data[data == 0] = 1
  data = ceiling(data)
  
  output = list(data, mean, DE_gene_idx, phi_g)
  names(output) = list("data", "mean",  "DE_idx", "dispersion")
  return(output)
}


# Implement DESeq2 method
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
  #dds = DESeqDataSetFromMatrix(data, DataFrame(treatment), ~treatment)
  
  if (is.null(norm_factors)){
    norm_factors = rep(1, ncol(data))
  }
  
  
  result_GLM = matrix(0, nrow(data), 4)
  
  for(i in 1: nrow(data)){
    
    genewise_data = data[i, ]
    genewise_dataframe = data.frame(read = genewise_data, treatment=treatment)
    
    if (dist == "qpois"){
      try({model = glm(read ~ treatment + offset(log(norm_factors)), family = quasipoisson(), data = genewise_dataframe)},silent = TRUE)
    } else if (dist == "negbin"){
      try({model = glm.nb(read ~ treatment + offset(log(norm_factors)), data = genewise_dataframe)},silent = TRUE)
      #model = glm.nb(read ~ treatment + offset(log(norm_factors)), data = genewise_dataframe)
    }
    #result = summary(glht(model, linfct=matrix(c(-1, 1), 1)))
    #result = summary(glht(model, mcp(treatment="Tukey")))
    #result_GLM[i, 1: 2] = c(coef(result), result$test$pvalues[1])
    
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


####perform comparison on simulation data 


Dispersion = matrix(c(0.1,0.2,0.05,0.1,0.02,0.05),nrow=3,byrow = TRUE)

FNR<-array(0, dim = c(3,4,10))
AUC<-array(0, dim = c(3,4,10))


for(i in 1:3)
  
{
  
  for (iter in 1: 10){
    
    
    set.seed(iter)
    
    
    ## generate data and true DE genes
    output = sim_data(12500, 2500, 6, 6, dispersion = Dispersion[i,])
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))
    
    
    # DESeq2 method
    
    output1 = deseq_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    #result1 = idx_criterion(idx_true, output1$DE_idx)
    FNR[i,1,iter] = length(intersect(which(idx_true==1),union(which(output1$DE_idx==0),which(is.na(output1$DE_idx)))))/sum(idx_true)
    AUC[i,1,iter] = auc(idx_true, output1$pvals)
    
    #auc(idx_true[which(!is.na(output1$pvals))], output1$pvals[which(!is.na(output1$pvals))])
    
    
    # edgeR method
    
    output2 = edgeR_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    #result2 = idx_criterion(idx_true, output2$DE_idx)
    #FN_2 = length(intersect(which(idx_true==1),which(output2$DE_idx==0)))/sum(idx_true)
    FNR[i,2,iter] = length(intersect(which(idx_true==1),union(which(output2$DE_idx==0),which(is.na(output2$DE_idx)))))/sum(idx_true)
    AUC[i,2,iter] = auc(idx_true, output2$pvals)
    
    # Proposed GLM method
    
    
    
    output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5, norm_factors=output2$norm_factors,dist = "negbin")
    #output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5,dist = "negbin")
    #result3 = idx_criterion(idx_true, output3$DE_idx)
    #FN_3 = length(intersect(which(idx_true==1),which(output3$DE_idx==0)))/sum(idx_true)
    FNR[i,3,iter] = length(intersect(which(idx_true==1),union(which(output3$DE_idx==0),which(is.na(output3$DE_idx)))))/sum(idx_true)
    
    #length(intersect(which(idx_true==0),which(output3$DE_idx==1)))/10000
    
    AUC[i,3,iter] = auc(idx_true, output3$pvals)
    
    # Limma
    
    design <- model.matrix(~ 0 + factor(c(rep("Control", 6), rep("Treatment", 6))))
    colnames(design) <- c("Control", "Treatment")
    #dgeListObj <- DGEList(data)
    #dgeListObj <- calcNormFactors(dgeListObj)
    #voomObj <- voom(dgeListObj, design)
    voomObj <- voom(data, design)
    fit <- lmFit(voomObj, design)
    fit <- eBayes(fit,proportion = 0.1)
    limma_res <- topTable(fit, coef=2, number=Inf,sort.by="none",p.value=0.05,lfc = 1.5)
    # To get the adjusted p-values and log fold change
    #limma_adj_p_values <- limma_res$adj.P.Val
    #limma_logFC <- limma_res$logFC
    #intersect(which(limma_adj_p_values<0.05),which(abs(limma_logFC)>1.5))
    FNR[i,4,iter] = length(intersect(which(idx_true==1),setdiff(seq(1,12500,1),as.numeric(row.names(limma_res)))))/sum(idx_true)
    
    #length(intersect(which(idx_true==0),as.numeric(row.names(limma_res))))/10000
    
    limma_res_2 <- topTable(fit, coef=2, number=Inf,sort.by="none")
    limma_adj_p_values <- limma_res_2$adj.P.Val
    AUC[i,4,iter] = auc(idx_true, limma_adj_p_values)
    
    
    print(iter)
    
  }
  
}


FNR_plot_data_NB <- data.frame(method = c(rep(c("DESeq2","EdgeR","GLM (NB)", "Limma"),3)),
                               setting = c(rep(1,4),rep(2,4),rep(3,4)), 
                               FNR = c(0.8922, 0.7304, 0.4476, 0.7839, 0.8601, 0.8563, 0.5966, 0.9501, 0.8955, 0.9624, 0.7537, 0.8857),
                               sd = c(0.0101, 0.010, 0.010, 0.014, 0.0070, 0.011, 0.0136, 0.0060, 0.0052, 0.0060, 0.0086, 0.0147))
rank_plot_data$num_ID = rep(-seq(-60,-1,2),3)

# 
ggplot(FNR_plot_data_NB, aes(x = setting, y = FNR)) +
  geom_errorbar(
    aes(ymin = FNR-sd, ymax = FNR+sd,color = method),
    position = position_dodge(.2), width = .7
  ) + scale_x_continuous(breaks=seq(1,3,1),labels=c(TeX(r'($0.02 < \theta_{g}^{NB} < 0.05$)'), TeX(r'($0.05 < \theta_{g}^{NB} < 0.1$)'), TeX(r'($0.1 < \theta_{g}^{NB} < 0.2$)')),limits = c(0.5, 3.5))  +
  theme_bw() + labs(title = "False negative rate comparison under negative binomial count", y = "False negative rate", x = "") + 
  geom_point(aes(color = method), position = position_dodge(.2), size = 2.5) +  
  scale_color_manual(values=c("DESeq2"="purple", "EdgeR"="blue", "GLM (NB)"="red","Limma"="orange")) + coord_cartesian(ylim=c(0,1)) +
  theme(axis.text.x= element_text(size=12.5),axis.text.y= element_text(size=13),axis.title = element_text(size=15),plot.title = element_text(size=13),
        legend.title = element_text(size=13), legend.text = element_text(size=12))





AUC_plot_data_NB <- data.frame(method = c(rep(c("DESeq2","EdgeR","GLM (NB)", "Limma"),3)),
                               setting = c(rep(1,4),rep(3,4),rep(5,4)), AUC = c(0.5026, 0.4873, 0.4878, 0.5169, 0.4907, 0.5046, 0.4928, 0.4482, 0.4826, 0.4911, 0.4670, 0.3463))
rank_plot_data$num_ID = rep(-seq(-60,-1,2),3)
ggplot(AUC_plot_data_NB, aes(fill=method, y=AUC, setting)) + 
  geom_bar(position="dodge", stat="identity",width=1) + 
  scale_x_continuous(breaks=seq(1,6,2),labels=c(TeX(r'($0.02 < \theta_{g}^{NB} < 0.05$)'), TeX(r'($0.05 < \theta_{g}^{NB} < 0.1$)'), TeX(r'($0.1 < \theta_{g}^{NB} < 0.2$)')),limits = c(0, 6))  +
  coord_cartesian(ylim=c(0,0.75)) + 
  labs(title = "AUC comparison under negative binomial count", y = "AUC", x = "") +  theme_bw() +
  theme(axis.text.x= element_text(size=14),axis.text.y= element_text(size=13),axis.title = element_text(size=15),plot.title = element_text(size=13),
        legend.title = element_text(size=13), legend.text = element_text(size=12)) + 
  scale_fill_manual(values=c("purple", "blue", "red","orange"))





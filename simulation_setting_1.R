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


## functions for performance measurement 
## criterion for choosing DE genes 

idx_criterion <- function(idx_true, idx_pred, FPR=T, FNR=T, PP=T){
  
  output = list()
  tab <- ConfusionMatrix(idx_pred, idx_true)
  if (ncol(tab) == 2){
    TP = tab[2, 2]
    FP = tab[1, 2]
    TN = tab[1, 1]
    FN = tab[2, 1]
  } else if ((ncol(tab) == 1) && (colnames(tab) == "0")){
    TP = 0
    FP = 0
    TN = tab[1, 1]
    FN = tab[2, 1]
  } else if ((ncol(tab) == 1) && (colnames(tab) == "1")){
    TP = tab[2, 1]
    FP = tab[1, 1]
    TN = 0
    FN = 0
  }
  if (FPR){
    fpr = FP / (FP + TN)
    output = append(output, fpr)
  }
  if (FNR){
    fnr = FN / (TP + FN)
    output = append(output, fnr)
  }
  if (PP) {
    pp = TP + FP
    output = append(output, pp)
  }
  
  return(output)
}

## Venn Diagram

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}




#Simulation setting 1 for negative Binominal:
#In the following, we set the working distribution for proposed DEHOGT 
#method as quasi-Poisson. To obtain the results from working distribution being 
#negative binomial, please replace output3 as follows: 
#output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5,dist = "negbin")


#dispersion: 0.1 - 0.2, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)

pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(0.1, 0.2))
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx)
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")

#Simulation setting 1 for negative Binominal:
#dispersion: 0.2 - 0.5, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)


pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(0.2, 0.5))
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx)
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")

#Simulation setting 1 for negative Binominal: 
#dispersion: 0.5 - 1, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)


pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(0.5, 1))
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx )
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")



#Simulation setting 1 for quasi-Poisson: 
#dispersion: 1 - 5, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)

pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(1, 5), dist = "qpois")
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx)
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")



#Simulation setting 1 for quasi-Poisson: 
#dispersion: 5 - 10, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)

pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(5, 10), dist = "qpois")
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx)
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")


#Simulation setting 1 for quasi-Poisson: 
#dispersion: 10 - 20, repeat 50 times, log2foldchange adjusted at 1.5, pvalue cutoff 0.05
fpr_mat = matrix(rep(0, 50 * 3), nrow=50)
fnr_mat = matrix(rep(0, 50 * 3), nrow=50)
pp_mat = matrix(rep(0, 50 * 3), nrow=50)
auc_mat = matrix(rep(0, 50 * 3), nrow=50)

pb <- progress_bar$new(total = 50)
pb$tick(0)
for (iter in 1: 50){

    pb$tick()
    set.seed(iter)

    ## generate data and true DE genes
    output = sim_data(12500, 2500, 0, 6, 6, dispersion = c(10, 20), dist = "qpois")
    data = output$data
    idx_true = rep(0, 12500)
    idx_true[output$DE_up_idx] = 1
    treatment = factor(c(rep(1, 6), rep(2, 6)))

    # DESeq method

    output1 = deseq_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result1 = idx_criterion(idx_true, output1$DE_idx)
    fpr_mat[iter, 1] = result1[[1]]
    fnr_mat[iter, 1] = result1[[2]]
    pp_mat[iter, 1] = result1[[3]]
    auc_mat[iter, 1] = auc(idx_true, output1$pvals)
    # edgeR method

    output2 = edgeR_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result2 = idx_criterion(idx_true, output2$DE_idx)
    fpr_mat[iter, 2] = result2[[1]]
    fnr_mat[iter, 2] = result2[[2]]
    pp_mat[iter, 2] = result2[[3]]
    auc_mat[iter, 2] = auc(idx_true, output2$pvals)
    # Proposed GLM method

    output3 = glm_func(data, treatment, pval_thre = 0.05,l2fc = T, l2fc_thre = 1.5)
    result3 = idx_criterion(idx_true, output3$DE_idx)
    fpr_mat[iter, 3] = result3[[1]]
    fnr_mat[iter, 3] = result3[[2]]
    pp_mat[iter, 3] = result3[[3]]
    auc_mat[iter, 3] = auc(idx_true, output3$pvals)
}

write.matrix(fpr_mat, "setting1_fpr.txt")
write.matrix(fnr_mat, "setting1_fnr.txt")
write.matrix(pp_mat, "setting1_pp.txt")
write.matrix(auc_mat, "setting1_auc.txt")


# visualization
# output4 denotes results from DEHOGT with negative binomial as working distribution

par(mfrow=c(1, 2))

plot(roc(idx_true, 1 - output1$pvals), col="purple", main="NB with overdispersion")
plot(roc(idx_true, 1 - output2$pvals), col="blue", add=TRUE)
plot(roc(idx_true, 1 - output3$pvals), col="orange", add=TRUE)
plot(roc(idx_true, 1 - output4$pvals), col="red", add=TRUE)
legend(0.4, 0.4, legend=c("DESeq", "edgeR", "DEHOGT-QP", "DEHOGT-NB"), col=c("purple", "blue", "orange", "red"), lty=1, bg="white")


plot(roc(idx_true, 1 - output1$pvals), col="purple", main="QP with overdispersion")
plot(roc(idx_true, 1 - output2$pvals), col="blue", add=TRUE)
plot(roc(idx_true, 1 - output3$pvals), col="orange", add=TRUE)
plot(roc(idx_true, 1 - output4$pvals), col="red", add=TRUE)
legend(0.4, 0.4, legend=c("DESeq", "edgeR", "DEHOGT-QP", "DEHOGT-NB"), col=c("purple", "blue", "orange", "red"), lty=1, bg="white")





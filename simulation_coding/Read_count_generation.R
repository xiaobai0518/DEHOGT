
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

#Generate quasi-Possion distribution random number

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
    mean[, j] = runif(N, min=0, max=M[j])
  }
  
  for (j in (S_1 + 1): S){
    for (i in 1: N){
      if (i %in% DE_gene_idx){
        mean[i, j] = rexp(1, rate=1/100) + M[j]
      } else{
        mean[i, j] = runif(1, min=0, max=M[j])
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




#Read count generating under simulation setting 2: low discrepancy of expression levels

sim_data2 <- function(N, N_DE, S_1, S_2, dispersion = c(0.5, 1), dist="negbin", seed=NULL){
  
  if (! is.null(seed)){
    set.seed(seed)
  }
  
  DE_gene_idx = sample(1: N, N_DE)
  sample_data = matrix(0, N, S_1 + S_2)
  phi = runif(N, min = dispersion[1], max = dispersion[2])
  
  for (i in 1: N){
    
    count_mean = ceiling(runif(S_1 + S_2, min = 2, max = 500))
    
    if (i %in% DE_gene_idx){
      mu = mean(count_mean[1: S_1]) + 1
      count_mean[(S_1 + 1): (S_1 + S_2)] = count_mean[(S_1 + 1): (S_1 + S_2)] + ceiling(mu * (exp(1.5) - 1))
    }
    
    if (dist == "negbin"){
      sample_data[i, ] = rnegbin(S_1 + S_2, mu = count_mean, phi[i])
    } else if (dist == "qpois"){
      sample_data[i, ] = rqpois(S_1 + S_2, mu = count_mean, phi[i])
    }
    
  }
  
  #sample_data = sample_data[rowSums(sample_data) > 0, ] + 1
  sample_data = ceiling(sample_data)
  
  output = list(sample_data, DE_gene_idx)
  names(output) = list("data", "DE_idx")
  
  return(output)
}





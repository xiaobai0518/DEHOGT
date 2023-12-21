
library(edgeR)
library(limma)
library(DESeq2)
library(MASS) # For glm.nb function

data <- read.delim("GSE116583_transplant.am.htseq.all.rpkm.txt", header=TRUE)

# Setting up count data
countData <- data[, 2:ncol(data)]
countData <- round(countData) 
samples <- colnames(countData)


control_group = ifelse(grepl("Naive", samples), 1, 0)
treatment_2h = ifelse(grepl("02H", samples), 1, 0)
treatment_24h = ifelse(grepl("24H", samples), 1, 0)
X = cbind(control_group,treatment_2h,treatment_24h)

X_factor = factor(Matrix::rowSums(X%*%diag(c(1:3))))

library(dplyr)
library(multcomp)


#edgeR

library(edgeR)
# Create a DGEList object
treatment = as.factor(c(1,1,1,1,2,2,2,2,3,3,3,3)) 
dge_edgeR <- DGEList(counts = countData,group = treatment)
colnames(X) <- c("Control", "Treatment_2h", "Treatment_24h")
design <- model.matrix(~0+X)
colnames(design) <- c("Control", "Treatment_2h", "Treatment_24h")

#dge_edgeR <- estimateDisp(dge_edgeR, design)
dge_edgeR <- estimateDisp(dge_edgeR)
#dge_edgeR = calcNormFactors(dge_edgeR)
dge_edgeR = estimateCommonDisp(dge_edgeR, verbose=T)
dge_edgeR = estimateTagwiseDisp(dge_edgeR)

edge_results = exactTest(dge_edgeR, pair=c(1,2))
p_values_2h_vs_Control = edge_results$table$PValue
p_values_2h_vs_Control = p.adjust(p_values_2h_vs_Control, method="BH")
logfold_2h_vs_Control = edge_results$table$logFC

edgeR2hgeneselected<- data$Symbol[intersect(which(p_values_2h_vs_Control<0.05),which(abs(logfold_2h_vs_Control)>1.5))]

edge_results = exactTest(dge_edgeR, pair=c(1,3))
p_values_24h_vs_Control = edge_results$table$PValue
p_values_24h_vs_Control = p.adjust(p_values_24h_vs_Control, method="BH")
logfold_24h_vs_Control = edge_results$table$logFC
edgeR24hgeneselected<- data$Symbol[intersect(which(p_values_24h_vs_Control<0.05),which(abs(logfold_24h_vs_Control)>1.5))]


#####DEHOGT


result_GLM_control_2h = matrix(0,nrow(data),2)
result_GLM_control_24h = matrix(0,nrow(data),2)
#K <- matrix(c(-1, 1, 0),nrow=1)

for (i in 1:nrow(countData)) {
  genewise_data <- countData[i,]
  genewise_dataframe <- data.frame(read = as.numeric(genewise_data), treatment = X_factor)
  b = genewise_dataframe[,1]
  if(sd(b[1:4]) + sd(b[5:8]) + sd(b[9:12]) ==0)
  {result_GLM_control_2h[i,] <- c(0, 1)
  next}
  else
    #{a = glm.nb(read ~ treatment - 1, data = genewise_dataframe)
    #a1 <- summary(glht(a, linfct = K)) 
    #result_GLM_control_2h[i,] <- c(coef(a1), a1$test$pvalues[1])}
  {
    model = glm.nb(read ~ treatment, data = genewise_dataframe)
    result_GLM_control_2h[i, 1: 2] = c(summary(model)$coefficients[2, 1], summary(model)$coefficients[2, 4])
    result_GLM_control_24h[i, 1: 2] = c(summary(model)$coefficients[3, 1], summary(model)$coefficients[3, 4])
  }
  
}

adjusted_pvalues <- p.adjust(result_GLM_control_2h[, 2], method = "BH")
#log2fold chang and adjusted p-values
result_GLM_control_2h_final <- data.frame(
  log2fold =  log(exp(abs(result_GLM_control_2h[, 1])), base=2),       
  pvalue = adjusted_pvalues
)
selected_genes_glm_2h <- data$Symbol[intersect(which(result_GLM_control_2h_final$pvalue < 0.05),which(abs(result_GLM_control_2h_final$log2fold)>1.5))]


adjusted_pvalues <- p.adjust(result_GLM_control_24h[, 2], method = "BH")
#log2fold chang and adjusted p-values
result_GLM_control_24h_final <- data.frame(
  log2fold =  log(exp(abs(result_GLM_control_24h[, 1])), base=2),       
  pvalue = adjusted_pvalues
)
selected_genes_glm_24h <- data$Symbol[intersect(which(result_GLM_control_24h_final$pvalue < 0.05),which(abs(result_GLM_control_24h_final$log2fold)>1.5))]


#DESeq2 

library(DESeq2)

colData <- DataFrame(
  condition = X_factor,
  groups = colnames(countData)
)

design <- model.matrix(~0+X_factor)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)
sizeFactors(dds) <- 1
dds <- DESeq(dds)

# Define contrasts
contrast_2h_vs_Control <- c(1, -1, 0)
contrast_24h_vs_Control <- c(1, 0, -1)

results_2h_vs_Control <- results(dds, contrast = contrast_2h_vs_Control)
results_24h_vs_Control <- results(dds, contrast = contrast_24h_vs_Control)

DESEq2_topGenes_2h_vs_Control <- data$Symbol[intersect(intersect(which(!is.na(results_2h_vs_Control$padj)),which(results_2h_vs_Control$padj < 0.05)), which(abs(results_2h_vs_Control$log2FoldChange)>1.5))] 
DESEq2_topGenes_24h_vs_Control <- data$Symbol[intersect(intersect(which(!is.na(results_24h_vs_Control$padj)),which(results_24h_vs_Control$padj < 0.05)), which(abs(results_24h_vs_Control$log2FoldChange)>1.5))] 


# Count the number of genes selected by each method
method_counts1 <- data.frame(
  Method = c("edgeR", "GLM", "DESeq2"),
  Genes_2h_vs_Control = c(length(edgeR2hgeneselected),length(selected_genes_glm_2h), 
                          length(DESEq2_topGenes_2h_vs_Control)
  ), 
  Genes_24h_vs_Control = c(length(edgeR24hgeneselected), length(selected_genes_glm_24h),
                           length(DESEq2_topGenes_24h_vs_Control)
  )
)
# Print the summary table
print(method_counts1)

setdiff(selected_genes_glm_2h,union(edgeR2hgeneselected,DESEq2_topGenes_2h_vs_Control))
GLM_select = setdiff(selected_genes_glm_24h,union(edgeR24hgeneselected,DESEq2_topGenes_24h_vs_Control))
write.csv(x=GLM_select, file="unique_DE_gene_GLM.csv")


#Venn Diagram 
#install.packages("VennDiagram")
update.packages("VennDiagram")
library(VennDiagram)
gene_sets <- list(
  edgeR = edgeR24hgeneselected,
  GLM = selected_genes_glm_24h,
  Deseq2 = DESEq2_topGenes_24h_vs_Control
)

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("EdgeR", "GLM", "Deseq2"),
  filename = NULL,  # To display the plot in R console
  #output = TRUE, 
  col = c("red", "green","orange"),
  #category.cex = 1,
  cex = 2,
  #fill = c("red", "green", "blue"),
  cat.cex=0.8,
  #cat.col = c(1,2,3),
  #cat.col = c("red","green","orange"),
  width.prop = 0.5
)
grid.newpage()
grid.draw(venn.plot)







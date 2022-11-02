

#####install and load dependent packages
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
BiocManager::install("edgeR")

library(DESeq2)
library(DESeqAnalysis)
library(openxlsx)
library(tibble)
library(data.table)
library(dplyr)
library(edgeR)
library(MASS)
library(Rmpfr)
library(VennDiagram)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(multcomp)

######

##load read count data and sample information data
Non_normalized_Read_Counts_combined_replicates_ <- read.csv("RNA_Read_Counts.csv")
Raw_count_data = Non_normalized_Read_Counts_combined_replicates_

Sample_Information <- read.csv("Sample Information.csv")
Sample_Information_coldata = Sample_Information[,3:5]
rownames(Sample_Information_coldata) = Sample_Information[,2]
Sample_Information_coldata$TreatTime <- factor(paste0(Sample_Information_coldata$Treatment,Sample_Information_coldata$Timepoint))

count = Raw_count_data[,3:45]
rownames(count) = Raw_count_data$GeneID
dds <- DESeqDataSetFromMatrix(count, colData = Sample_Information_coldata, design = ~ Treatment)
dds <- DESeq(dds)


####extract sample index for each treatment
corthigh_group = which(colnames(Raw_count_data)%in%c('corth2_S20_L001','corth3_S12_L001','corth4_S19_L001',
                                                        'corthwo2_S11_L001','corthwo3_S13_L001','corthwo4_S4_L001'))
cortlow_group = which(colnames(Raw_count_data)%in%c('cortl2_S14_L001','cortl3_S19_L001','cortl4_S13_L001','cortlwo2_S18_L001','cortlwo3_S6_L001','cortlwo4_S2_L001'))
cortvechhigh_group = which(colnames(Raw_count_data)%in%c('cvh1_S21_L001','cvh2_S16_L001','cvh3_S14_L001','cvhwo1_S1_L001'))
cortvechlow_group = which(colnames(Raw_count_data)%in%c('cvl1_S18_L001','cvl2_S11_L001','cvl4_S7_L001','cvlwo2_S1_L001'))
dexhigh_group = which(colnames(Raw_count_data)%in%c('dexh2_S15_L001','dexh3_S20_L001','dexh4_S7_L001','dexhwo2_S12_L001','dexhwo3_S17_L001','dexhwo4_S5_L001'))
dexlow_group = which(colnames(Raw_count_data)%in%c('dexl1_S16_L001','dexl2_S15_L001','dexl3_S9_L001','dexlwo1_S22_L001','dexlwo2_S3_L001','dexlwo3_S21_L001','dexl1_S16_L001'))
dexvechhigh_group = which(colnames(Raw_count_data)%in%c('dvh1_S4_L001','dvh2_S8_L001','dvh4_S10_L001','dvhwo4_S9_L001'))
dexvechlow_group = which(colnames(Raw_count_data)%in%c('dvl1_S3_L001','dvl2_S6_L001','dvl3_S2_L001','dvlwo3_S17_L001'))

####extract sample index for time point 3
Time3_group = which(colnames(Raw_count_data)%in%c("corth2_S20_L001","corth3_S12_L001", "corth4_S19_L001", "cortl2_S14_L001", "cortl3_S19_L001",
                                                     "cortl4_S13_L001","cvh1_S21_L001",
                                                     "cvh2_S16_L001","cvh3_S14_L001","cvl1_S18_L001","cvl2_S11_L001","cvl4_S7_L001",   
                                                     "dexh2_S15_L001","dexh3_S20_L001","dexh4_S7_L001","dexl1_S16_L001","dexl2_S15_L001", 
                                                     "dexl3_S9_L001","dvh1_S4_L001","dvh2_S8_L001","dvh4_S10_L001","dvl1_S3_L001",   
                                                     "dvl2_S6_L001","dvl3_S2_L001"))
####extract sample index for time point 6
Time6_group = which(colnames(Raw_count_data)%in%c('corthwo2_S11_L001','corthwo3_S13_L001','corthwo4_S4_L001','cortlwo2_S18_L001','cortlwo3_S6_L001','cortlwo4_S2_L001','cvhwo1_S1_L001',
                                                     'cvlwo2_S1_L001','dexhwo2_S12_L001','dexhwo3_S17_L001','dexhwo4_S5_L001','dexlwo1_S22_L001','dexlwo2_S3_L001','dexlwo3_S21_L001','dvhwo4_S9_L001','dvlwo3_S17_L001'))

####extract sample index for control group
control_group = which(colnames(Raw_count_data)%in%c('Ctl1_S8_L001','ctl3_S10_L001','ctl4_S5_L001'))



###generate treatment index dummy matrix
corthigh_group_3 = rep(0,43)
corthigh_group_3[intersect(corthigh_group,Time3_group)-1] = 1
corthigh_group_6 = rep(0,43)
corthigh_group_6[intersect(corthigh_group,Time6_group)-1] = 1

cortlow_group_3 = rep(0,43)
cortlow_group_3[intersect(cortlow_group,Time3_group)-1] = 1
cortlow_group_6 = rep(0,43)
cortlow_group_6[intersect(cortlow_group,Time6_group)-1] = 1

cortvechhigh_group_3 = rep(0,43)
cortvechhigh_group_3[intersect(cortvechhigh_group,Time3_group)-1] = 1
cortvechhigh_group_6 = rep(0,43)
cortvechhigh_group_6[intersect(cortvechhigh_group,Time6_group)-1] = 1

cortvechlow_group_3 = rep(0,43)
cortvechlow_group_3[intersect(cortvechlow_group,Time3_group)-1] = 1
cortvechlow_group_6 = rep(0,43)
cortvechlow_group_6[intersect(cortvechlow_group,Time6_group)-1] = 1

dexhigh_group_3 = rep(0,43)
dexhigh_group_3[intersect(dexhigh_group,Time3_group)-1] = 1
dexhigh_group_6 = rep(0,43)
dexhigh_group_6[intersect(dexhigh_group,Time6_group)-1] = 1

dexlow_group_3 = rep(0,43)
dexlow_group_3[intersect(dexlow_group,Time3_group)-1] = 1
dexlow_group_6 = rep(0,43)
dexlow_group_6[intersect(dexlow_group,Time6_group)-1] = 1

dexvechhigh_group_3 = rep(0,43)
dexvechhigh_group_3[intersect(dexvechhigh_group,Time3_group)-1] = 1
dexvechhigh_group_6 = rep(0,43)
dexvechhigh_group_6[intersect(dexvechhigh_group,Time6_group)-1] = 1

dexvechlow_group_3 = rep(0,43)
dexvechlow_group_3[intersect(dexvechlow_group,Time3_group)-1] = 1
dexvechlow_group_6 = rep(0,43)
dexvechlow_group_6[intersect(dexvechlow_group,Time6_group)-1] = 1

control_group_3 = rep(0,43)
control_group_3[control_group-1] = 1

X = cbind(corthigh_group_3,corthigh_group_6,cortlow_group_3,cortlow_group_6,cortvechhigh_group_3,cortvechhigh_group_6,
          cortvechlow_group_3,cortvechlow_group_6,dexhigh_group_3,dexhigh_group_6,dexlow_group_3,dexlow_group_6,
          dexvechhigh_group_3,dexvechhigh_group_6,dexvechlow_group_3,dexvechlow_group_6,control_group_3)

X[,which(colnames(X) == 'cortvechhigh_group_3')]
X[,which(colnames(X) == 'cortvechhigh_group_6')]
X[,which(colnames(X) == 'cortvechlow_group_3')]
X[,which(colnames(X) == 'cortvechlow_group_6')]
X[,which(colnames(X) == 'dexvechhigh_group_3')]
X[,which(colnames(X) == 'dexvechhigh_group_6')]
X[,which(colnames(X) == 'dexvechlow_group_6')]

########obtain samplewise normalization constants via TMM method
X_factor = factor(Matrix::rowSums(X%*%diag(c(1:17))))
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
treatmentwise_factor_2 = edgedata@.Data[[2]]$norm.factors




#####treatment comparison: dex high 3 vs control 3
result_dexh3_con = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(rep(0,8),1,rep(0,7),-1),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read ~ treatment -1 + offset(log(treatmentwise_factor_2)), data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_dexh3_con[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_dexh3_con[,2], method = "BH")

result_dexh3_con_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_dexh3_con[,1], pvalue = p)

result_dexh3_con_2 = dplyr::filter(result_dexh3_con_1,pvalue<0.05)

result_dexh3_con_3 = dplyr::filter(result_dexh3_con_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method
c_1 = result_dexh3_con_3[order(result_dexh3_con_3[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"dexhigh_vs_control_time3.csv", row.names = FALSE)


#load the results of DE genes from DESeq2 method
TreatTime_dex.high3_Vs_control <- read.csv("Previous_output/DESeq2_2Strand_Filtered50/dex.high3_Vs_control3.csv")

#edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_dexh3_con = exactTest(edgedata, pair=c(9, 17))
log2fold2 = edge_dexh3_con$table$logFC

edge_dexh3_con$table$PValue = p.adjust(edge_dexh3_con$table$PValue, method="BH")

edgeR_select_dexh3_con = dplyr::filter(edge_dexh3_con$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 12
venn.diagram(
  x = list(result_dexh3_con_3$GeneID,TreatTime_dex.high3_Vs_control$X,rownames(edgeR_select_dexh3_con)),
  #list(GLM = result_dexh3_con_3$GeneID, DESeq = TreatTime_dex.high3_Vs_control$X, EdgeR = rownames(edgeR_select_dexh3_con)),
  category.names = c('GLM (NB)',"DESeq","EdgeR"),
  filename = 'dexh3_control_3.png',
  output = TRUE,
  imagetype="png" ,
  height = 380 , 
  width = 380 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.2,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Dex high 3 vs Control",
  main.cex = 0.5,
  sub.cex = 0.5 
)

#draw the p-value rank plot in figure 19
overlap_gene_2 = intersect(intersect(result_dexh3_con_3$GeneID,TreatTime_dex.high3_Vs_control$X),rownames(edgeR_select_dexh3_con))

a_1 = get_gene_annotation(overlap_gene_2)
overlap_gene_2 = a_1$ENSEMBL
overlap_gene_reorder_2 = overlap_gene_2[order(result_dexh3_con_3[which(result_dexh3_con_3$GeneID%in%overlap_gene_2),3])]
overlap_gene_reorder_2 = overlap_gene_reorder_2[1:30]


glm_select_2 = result_dexh3_con_3[which(result_dexh3_con_3$GeneID%in%overlap_gene_reorder_2),c(1,3)]
overlap_gene_reorder_2 = glm_select_2[order(glm_select_2[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank_2 = rank(glm_select_2[order(glm_select_2[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name_2 = get_gene_annotation(overlap_gene_reorder_2)

deq_select_2 = TreatTime_dex.high3_Vs_control[which(TreatTime_dex.high3_Vs_control$X%in%overlap_gene_reorder_2),c(1,6)]
deq_order_2 = match(overlap_gene_reorder_2,deq_select_2[,1])
deq_rank_2 = rank(deq_select_2[deq_order_2,2],ties.method = 'first')

edgeR_select_dexh3_con$names = rownames(edgeR_select_dexh3_con)
edge_select_dexh3_con = edgeR_select_dexh3_con[which(edgeR_select_dexh3_con$names%in%overlap_gene_reorder_2),c(3,4)]
edge_order_2 = match(overlap_gene_reorder_2,edge_select_dexh3_con[,2])
edge_rank_2 = rank(edge_select_dexh3_con[edge_order_2,1],ties.method = 'first')

rank_plot_data_2 <- data.frame(ID = rep(overlap_gene_reorder_2,3),rank = c(glm_rank_2,deq_rank_2,edge_rank_2),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data_2$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data_2, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name_2$SYMBOL) + coord_flip() +
  labs(title = "Dex high 3 vs control", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))

#####check the genes not in the overlap 
a = result_dexh3_con_3$GeneID
b = TreatTime_dex.high3_Vs_control$X
c = rownames(edgeR_select_dexh3_con)

GLM_uni = setdiff(a, union(b,c))
Des_uni = setdiff(b, union(a,c))
Edge_uni = setdiff(c, union(a,b))

GLM_uni_p = rep(0,length(GLM_uni))
for(i in 1:length(GLM_uni))
{
  GLM_uni_p[i] = result_dexh3_con_3[which(result_dexh3_con_3$GeneID== GLM_uni[i]),3]
}
a_1 = get_gene_annotation(GLM_uni[order(rank(GLM_uni_p,ties.method = 'first'))])

Des_uni_p = rep(0,length(Des_uni))
for(i in 1:length(Des_uni))
{
  Des_uni_p[i] = TreatTime_dex.high3_Vs_control[which(TreatTime_dex.high3_Vs_control$X== Des_uni[i]),6]
}
b_1 = get_gene_annotation(Des_uni[order(rank(Des_uni_p,ties.method = 'first'))])


Edge_uni_p = rep(0,length(Edge_uni))
for(i in 1:length(Edge_uni))
{
  Edge_uni_p[i] = edgeR_select_dexh3_con[which(rownames(edgeR_select_dexh3_con)== Edge_uni[i]),3]
}
c_1 = get_gene_annotation(Edge_uni[order(rank(Edge_uni_p,ties.method = 'first'))])





#####treatment comparison: dex high 6 vs control 3
result_dexhigh6_control = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(rep(0,9), 1, 0, 0, 0, 0, 0,0,-1),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read/treatmentwise_factor_2 ~ treatment - 1, data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_dexhigh6_control[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_dexhigh6_control[,2], method = "BH")

result_dexhigh6_control_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_dexhigh6_control[,1], pvalue = p)

result_dexhigh6_control_2 = filter(result_dexhigh6_control_1,pvalue<0.05)

result_dexhigh6_control_3 = filter(result_dexhigh6_control_2,abs(log2fold)>1.5)


GLM_edge_R_NB_result = result_dexhigh6_control_3

#obtain and store the names for the top 30 DE genes selected by proposed method
c_1 = GLM_edge_R_NB_result[order(GLM_edge_R_NB_result[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(GLM_edge_R_NB_result[order(GLM_edge_R_NB_result[,3])[1:30],1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"dexhigh_vs_control_time6.csv", row.names = FALSE)

#load the results of DE genes from DESeq2 method
TreatTime_dex.high6_vs_control <- read.csv("dex.high6_Vs_control3.csv")
#edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_results = exactTest(edgedata, pair=c(10, 17))
log2fold2 = edge_results$table$logFC

edge_results$table$PValue = p.adjust(edge_results$table$PValue, method="BH")
edgeR_select = filter(edge_results$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 13
library(VennDiagram)
library(RColorBrewer)
venn.diagram(
  x = list(result_dexhigh6_control_3$GeneID, TreatTime_dex.high6_vs_control$X, rownames(edgeR_select)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'dexhigh6_control.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #col=c("blue", 'red', '#00CC00'),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Dex High 6 vs Control"
)

#draw the p-value rank plot in figure 20

overlap_gene = intersect(intersect(result_dexhigh6_control_3$GeneID,TreatTime_dex.high6_vs_control$X),rownames(edgeR_select))
a= get_gene_annotation(overlap_gene)
overlap_gene = a$ENSEMBL
overlap_gene_reorder = overlap_gene[order(result_dexhigh6_control_3[which(result_dexhigh6_control_3$GeneID%in%overlap_gene),3])]
overlap_gene_reorder = overlap_gene_reorder[1:30]

glm_select = result_dexhigh6_control_3[which(result_dexhigh6_control_3$GeneID%in%overlap_gene_reorder),c(1,3)]
overlap_gene_reorder = glm_select[order(glm_select[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank = rank(glm_select[order(glm_select[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name = get_gene_annotation(overlap_gene_reorder)


deq_select = TreatTime_dex.high6_vs_control[which(TreatTime_dex.high6_vs_control$X%in%overlap_gene_reorder),c(1,6)]
deq_order = match(overlap_gene_reorder,deq_select[,1])
deq_rank = rank(deq_select[deq_order,2],ties.method = 'first')

edgeR_select$names = rownames(edgeR_select)
edge_select = edgeR_select[which(edgeR_select$names%in%overlap_gene_reorder),c(3,4)]
edge_order = match(overlap_gene_reorder,edge_select[,2])
edge_rank = rank(edge_select[edge_order,1],ties.method = 'first')

rank_plot_data <- data.frame(ID = rep(overlap_gene_reorder,3),rank = c(glm_rank,deq_rank,edge_rank),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name$SYMBOL) + coord_flip() +
  labs(title = "Dex high 6 vs control ", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))

#####check the genes not in the overlap 
nonoverlap_GLM_1 = GLM_edge_R_NB_result[which(GLM_edge_R_NB_result$GeneID%in%setdiff(GLM_edge_R_NB_result$GeneID,rownames(edgeR_select))),]
nonoverlap_GLM_1[order(rank(nonoverlap_GLM_1$pvalue,ties.method = 'first')),]
##top 20 genes
nonoverlap_GLM_1_name = get_gene_annotation(nonoverlap_GLM_1[order(rank(nonoverlap_GLM_1$pvalue,ties.method = 'first'))[1:20],1])
nonoverlap_edge_1 = edgeR_select[which(rownames(edgeR_select)%in%setdiff(rownames(edgeR_select),GLM_edge_R_NB_result$GeneID)),]
nonoverlap_edge_1[order(rank(nonoverlap_edge_1$PValue,ties.method = 'first')),]
nonoverlap_edge_1_name = get_gene_annotation(rownames(nonoverlap_edge_1[order(rank(nonoverlap_edge_1$PValue,ties.method = 'first')),])[1:20])


#####treatment comparison: cort high 3 vs control 3
result_corthigh3_control = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(1,rep(0,15),-1),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read/treatmentwise_factor_2 ~ treatment - 1, data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_corthigh3_control[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_corthigh3_control[,2], method = "BH")

result_corthigh3_control_1 = data.frame(GeneID = rownames(counts(dds)), log2fold =  log(exp(1),base = 2)*result_corthigh3_control[,1], pvalue = p)

result_corthigh3_control_2 = filter(result_corthigh3_control_1,pvalue<0.05)

result_corthigh3_control_3 = filter(result_corthigh3_control_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method
GLM_edge_R_NB_result_2 = result_corthigh3_control_3

c_1 = GLM_edge_R_NB_result_2[order(GLM_edge_R_NB_result_2[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"corthigh_vs_control_time3.csv", row.names = FALSE)


#load the results of DE genes from DESeq2 method
TreatTime_cort.high3_vs_control <- read.csv("Previous_output/DESeq2_2Strand_Filtered50/cort.high3_Vs_control3.csv")

# edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_results = exactTest(edgedata, pair=c(1, 17))
log2fold2 = edge_results$table$logFC

edge_results$table$PValue = p.adjust(edge_results$table$PValue, method="BH")
edgeR_select_2 = filter(edge_results$table,abs(logFC)>1.5&PValue<0.05)
intersect(result_corthigh3_control_3$GeneID,rownames(edgeR_select_2))


#draw the VennDiagram plot for figure 14
venn.diagram(
  x = list(result_corthigh3_control_3$GeneID, TreatTime_cort.high3_vs_control$X, rownames(edgeR_select_2)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'corthigh3_control.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Cort high 3 vs Control "
)


#draw the p-value rank plot in figure 21
overlap_gene_2 = intersect(intersect(result_corthigh3_control_3$GeneID,TreatTime_cort.high3_vs_control$X),rownames(edgeR_select_2))
a_1 = get_gene_annotation(overlap_gene_2)
overlap_gene_2 = a_1$ENSEMBL

overlap_gene_reorder_2 = overlap_gene_2[order(result_corthigh3_control_3[which(result_corthigh3_control_3$GeneID%in%overlap_gene_2),3])]
overlap_gene_reorder_2 = overlap_gene_reorder_2[1:30]


glm_select_2 = result_corthigh3_control_3[which(result_corthigh3_control_3$GeneID%in%overlap_gene_reorder_2),c(1,3)]
overlap_gene_reorder_2 = glm_select_2[order(glm_select_2[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank_2 = rank(glm_select_2[order(glm_select_2[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name_2 = get_gene_annotation(overlap_gene_reorder_2)

deq_select_2 = TreatTime_cort.high3_vs_control[which(TreatTime_cort.high3_vs_control$X%in%overlap_gene_reorder_2),c(1,6)]
deq_order_2 = match(overlap_gene_reorder_2,deq_select_2[,1])
deq_rank_2 = rank(deq_select_2[deq_order_2,2],ties.method = 'first')

edgeR_select_2$names = rownames(edgeR_select_2)
edge_select_2 = edgeR_select_2[which(edgeR_select_2$names%in%overlap_gene_reorder_2),c(3,4)]
edge_order_2 = match(overlap_gene_reorder_2,edge_select_2[,2])
edge_rank_2 = rank(edge_select_2[edge_order_2,1],ties.method = 'first')

rank_plot_data_2 <- data.frame(ID = rep(overlap_gene_reorder_2,3),rank = c(glm_rank_2,deq_rank_2,edge_rank_2),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data_2$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data_2, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name_2$SYMBOL) + coord_flip() +
  labs(title = "Cort high 3 vs control ", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))



#####check the genes not in the overlap 
nonoverlap_GLM_corth3_con = GLM_edge_R_NB_result_2[which(GLM_edge_R_NB_result_2$GeneID%in%setdiff(GLM_edge_R_NB_result_2$GeneID,rownames(edgeR_select_2))),]
nonoverlap_GLM_corth3_con[order(rank(nonoverlap_GLM_corth3_con$pvalue,ties.method = 'first')),]
##top 20 genes
nonoverlap_GLM_corth3_con_name = get_gene_annotation(nonoverlap_GLM_corth3_con[order(rank(nonoverlap_GLM_corth3_con$pvalue,ties.method = 'first')),1][1:20])

nonoverlap_edge_corth3_con = edgeR_select_2[which(rownames(edgeR_select_2)%in%setdiff(rownames(edgeR_select_2),GLM_edge_R_NB_result_2$GeneID)),]
nonoverlap_edge_corth3_con[order(rank(nonoverlap_edge_corth3_con$PValue,ties.method = 'first')),]
nonoverlap_edge_corth3_con_name = get_gene_annotation(rownames(nonoverlap_edge_corth3_con[order(rank(nonoverlap_edge_corth3_con$PValue,ties.method = 'first')),])[1:20])




#####treatment comparison: cort high 6 vs control
result_corth6_con = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(0,1,rep(0,14),-1),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read ~ treatment -1 + offset(log(treatmentwise_factor_2)), data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_corth6_con[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_corth6_con[,2], method = "BH")

result_corth6_con_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_corth6_con[,1], pvalue = p)

result_corth6_con_2 = dplyr::filter(result_corth6_con_1,pvalue<0.05)

result_corth6_con_3 = dplyr::filter(result_corth6_con_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method
c_1 = result_corth6_con_3[order(result_corth6_con_3[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"corthigh_vs_control_time6.csv", row.names = FALSE)

#load the results of DE genes from DESeq2 method
TreatTime_cort.high6_Vs_control <- read.csv("cort.high6_Vs_control3.csv")

#edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_corth6_con = exactTest(edgedata, pair=c(2, 17))
log2fold2 = edge_corth6_con$table$logFC

edge_corth6_con$table$PValue = p.adjust(edge_corth6_con$table$PValue, method="BH")

edgeR_select_corth6_con = dplyr::filter(edge_corth6_con$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 15
venn.diagram(
  x = list(result_corth6_con_3$GeneID,TreatTime_cort.high6_Vs_control$X,rownames(edgeR_select_corth6_con)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'corth6_control_3.png',
  output = TRUE ,
  imagetype="png" ,
  height = 380 , 
  width = 380 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Cort high 6 vs Control",
  main.cex = 0.5,
  sub.cex = 0.5 
)


#####check the genes not in the overlap 
a = result_corth6_con_3$GeneID
b = TreatTime_cort.high6_Vs_control$X
c = rownames(edgeR_select_corth6_con)

GLM_uni = setdiff(a, union(b,c))
Des_uni = setdiff(b, union(a,c))
Edge_uni = setdiff(c, union(a,b))

GLM_uni_p = rep(0,length(GLM_uni))
for(i in 1:length(GLM_uni))
{
  GLM_uni_p[i] = result_corth6_con_3[which(result_corth6_con_3$GeneID== GLM_uni[i]),3]
}
a_1 = get_gene_annotation(GLM_uni[order(rank(GLM_uni_p,ties.method = 'first'))])

Des_uni_p = rep(0,length(Des_uni))
for(i in 1:length(Des_uni))
{
  Des_uni_p[i] = TreatTime_cort.high6_Vs_control[which(TreatTime_cort.high6_Vs_control$X== Des_uni[i]),6]
}
b_1 = get_gene_annotation(Des_uni[order(rank(Des_uni_p,ties.method = 'first'))])


Edge_uni_p = rep(0,length(Edge_uni))
for(i in 1:length(Edge_uni))
{
  Edge_uni_p[i] = edgeR_select_corth6_con[which(rownames(edgeR_select_corth6_con)== Edge_uni[i]),3]
}
c_1 = get_gene_annotation(Edge_uni[order(rank(Edge_uni_p,ties.method = 'first'))])



#####treatment comparison: dex vehicle high 3 vs dex high 3
result_dexvh_dexh = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(rep(0,8),1,rep(0,3),-1,rep(0,4)),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read/treatmentwise_factor_2 ~ treatment - 1, data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_dexvh_dexh[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_dexvh_dexh[,2], method = "BH")

result_dexvh_dexh_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_dexvh_dexh[,1], pvalue = p)

result_dexvh_dexh_2 = filter(result_dexvh_dexh_1,pvalue<0.05)

result_dexvh_dexh_3 = filter(result_dexvh_dexh_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method
c_1 = result_dexvh_dexh_3[order(result_dexvh_dexh_3[,3])[1:383],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"dexvechhigh_vs_dexhigh_time3.csv", row.names = FALSE)

#load the results of DE genes from DESeq2 method
TreatTime_dexv.high3_vs_dex.high3 <- read.csv("dex.high3_Vs_dex.vehicle.high3.csv")

# edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_dexvh_dexh = exactTest(edgedata, pair=c(9, 13))
log2fold2 = edge_dexvh_dexh$table$logFC

edge_dexvh_dexh$table$PValue = p.adjust(edge_dexvh_dexh$table$PValue, method="BH")

edgeR_select_dexvh_dexh = filter(edge_dexvh_dexh$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 16
library(scales)
venn.diagram(
  x = list(result_dexvh_dexh_3$GeneID, TreatTime_dexv.high3_vs_dex.high3$X, rownames(edgeR_select_dexvh_dexh)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'dexvh_dexh_3.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Dex vehicle high 3 vs Dex high 3 ",
  main.cex = 0.5,sub.cex = 1.5 
)


#draw the p-value rank plot in figure 22
overlap_gene_2 = intersect(intersect(result_dexvh_dexh_3$GeneID,TreatTime_dexv.high3_vs_dex.high3$X),rownames(edgeR_select_dexvh_dexh))
a_1 = get_gene_annotation(overlap_gene_2)
overlap_gene_2 = a_1$ENSEMBL
overlap_gene_reorder_2 = overlap_gene_2[order(result_dexvh_dexh_3[which(result_dexvh_dexh_3$GeneID%in%overlap_gene_2),3])]
overlap_gene_reorder_2 = overlap_gene_reorder_2[1:30]


glm_select_2 = result_dexvh_dexh_3[which(result_dexvh_dexh_3$GeneID%in%overlap_gene_reorder_2),c(1,3)]
overlap_gene_reorder_2 = glm_select_2[order(glm_select_2[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank_2 = rank(glm_select_2[order(glm_select_2[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name_2 = get_gene_annotation(overlap_gene_reorder_2)

deq_select_2 = TreatTime_dexv.high3_vs_dex.high3[which(TreatTime_dexv.high3_vs_dex.high3$X%in%overlap_gene_reorder_2),c(1,6)]
deq_order_2 = match(overlap_gene_reorder_2,deq_select_2[,1])
deq_rank_2 = rank(deq_select_2[deq_order_2,2],ties.method = 'first')

edgeR_select_dexvh_dexh$names = rownames(edgeR_select_dexvh_dexh)
edge_select_dexvh_dexh = edgeR_select_dexvh_dexh[which(edgeR_select_dexvh_dexh$names%in%overlap_gene_reorder_2),c(3,4)]
edge_order_2 = match(overlap_gene_reorder_2,edge_select_dexvh_dexh[,2])
edge_rank_2 = rank(edge_select_dexvh_dexh[edge_order_2,1],ties.method = 'first')

rank_plot_data_2 <- data.frame(ID = rep(overlap_gene_reorder_2,3),rank = c(glm_rank_2,deq_rank_2,edge_rank_2),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data_2$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data_2, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name_2$SYMBOL) + coord_flip() +
  labs(title = "Dex vehicle high 3 vs Dex high 3", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))

#####check the genes not in the overlap 

nonoverlap_GLM_dexvdex_con = result_dexvh_dexh_3[which(result_dexvh_dexh_3$GeneID%in%setdiff(result_dexvh_dexh_3$GeneID,rownames(edgeR_select_dexvh_dexh))),]
nonoverlap_GLM_dexvdex_con[order(rank(nonoverlap_GLM_dexvdex_con$pvalue,ties.method = 'first')),]
##top 20 genes
nonoverlap_GLM_dexvdex_con_name = get_gene_annotation(nonoverlap_GLM_dexvdex_con[order(rank(nonoverlap_GLM_dexvdex_con$pvalue,ties.method = 'first')),1][1:20])

nonoverlap_edge_dexvdex_con = edgeR_select_dexvh_dexh[which(rownames(edgeR_select_dexvh_dexh)%in%setdiff(rownames(edgeR_select_dexvh_dexh),result_dexvh_dexh_3$GeneID)),]
nonoverlap_edge_dexvdex_con[order(rank(nonoverlap_edge_dexvdex_con$PValue,ties.method = 'first')),]
nonoverlap_edge_dexvdex_con_name = get_gene_annotation(rownames(nonoverlap_edge_dexvdex_con[order(rank(nonoverlap_edge_dexvdex_con$PValue,ties.method = 'first')),])[1:20])



#####treatment comparison: dex vehcile low 3 vs dex low 3
result_dexvl_dexl = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(rep(0,10),1,rep(0,4),-1,rep(0,1)),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read/treatmentwise_factor_2 ~ treatment - 1, data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_dexvl_dexl[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_dexvl_dexl[,2], method = "BH")

result_dexvl_dexl_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_dexvl_dexl[,1], pvalue = p)

result_dexvl_dexl_2 = filter(result_dexvl_dexl_1,pvalue<0.05)

result_dexvl_dexl_3 = filter(result_dexvl_dexl_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method

c_1 = result_dexvl_dexl_3[order(result_dexvl_dexl_3[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"dexvechlow_vs_dexlow_time3.csv", row.names = FALSE)

#load the results of DE genes from DESeq2 method
TreatTime_dexv.low3_vs_dex.low3 <- read.csv("dex.low3_Vs_dex.vehicle.low3.csv")

# edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_dexvl_dexl = exactTest(edgedata, pair=c(11, 16))
log2fold2 = edge_dexvl_dexl$table$logFC

edge_dexvl_dexl$table$PValue = p.adjust(edge_dexvl_dexl$table$PValue, method="BH")

edgeR_select_dexvl_dexl = filter(edge_dexvl_dexl$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 17
library(scales)
venn.diagram(
  x = list(result_dexvl_dexl_3$GeneID, TreatTime_dexv.low3_vs_dex.low3$X, rownames(edgeR_select_dexvl_dexl)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'dexvl_dexl_3.png',
  output = TRUE ,
  imagetype="png" ,
  height = 380 , 
  width = 380 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Dex vehicle low 3 vs Dex low 3 ",
  main.cex = .5,
  sub.cex = 0.5 
)

#draw the p-value rank plot in figure 23
overlap_gene_2 = intersect(intersect(result_dexvl_dexl_3$GeneID,TreatTime_dexv.low3_vs_dex.low3$X),rownames(edgeR_select_dexvl_dexl))
a_1 = get_gene_annotation(overlap_gene_2)
overlap_gene_2 = a_1$ENSEMBL
overlap_gene_reorder_2 = overlap_gene_2[order(result_dexvl_dexl_3[which(result_dexvl_dexl_3$GeneID%in%overlap_gene_2),3])]
overlap_gene_reorder_2 = overlap_gene_reorder_2[1:30]


glm_select_2 = result_dexvl_dexl_3[which(result_dexvl_dexl_3$GeneID%in%overlap_gene_reorder_2),c(1,3)]
overlap_gene_reorder_2 = glm_select_2[order(glm_select_2[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank_2 = rank(glm_select_2[order(glm_select_2[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name_2 = get_gene_annotation(overlap_gene_reorder_2)

deq_select_2 = TreatTime_dexv.low3_vs_dex.low3[which(TreatTime_dexv.low3_vs_dex.low3$X%in%overlap_gene_reorder_2),c(1,6)]
deq_order_2 = match(overlap_gene_reorder_2,deq_select_2[,1])
deq_rank_2 = rank(deq_select_2[deq_order_2,2],ties.method = 'first')

edgeR_select_dexvl_dexl$names = rownames(edgeR_select_dexvl_dexl)
edge_select_dexvl_dexl = edgeR_select_dexvl_dexl[which(edgeR_select_dexvl_dexl$names%in%overlap_gene_reorder_2),c(3,4)]
edge_order_2 = match(overlap_gene_reorder_2,edge_select_dexvl_dexl[,2])
edge_rank_2 = rank(edge_select_dexvl_dexl[edge_order_2,1],ties.method = 'first')

rank_plot_data_2 <- data.frame(ID = rep(overlap_gene_reorder_2,3),rank = c(glm_rank_2,deq_rank_2,edge_rank_2),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data_2$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data_2, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name_2$SYMBOL) + coord_flip() +
  labs(title = "Dex vehicle low 3 vs Dex low 3", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))

#####check the genes not in the overlap 

nonoverlap_GLM_dexvdex_con = result_dexvl_dexl_3[which(result_dexvl_dexl_3$GeneID%in%setdiff(result_dexvl_dexl_3$GeneID,rownames(edgeR_select_dexvl_dexl))),]
nonoverlap_GLM_dexvdex_con[order(rank(nonoverlap_GLM_dexvdex_con$pvalue,ties.method = 'first')),]
##top 20 genes
nonoverlap_GLM_dexvdex_con_name = get_gene_annotation(nonoverlap_GLM_dexvdex_con[order(rank(nonoverlap_GLM_dexvdex_con$pvalue,ties.method = 'first')),1][1:20])

nonoverlap_edge_dexvdex_con = edgeR_select_dexvh_dexh[which(rownames(edgeR_select_dexvh_dexh)%in%setdiff(rownames(edgeR_select_dexvh_dexh),result_dexvh_dexh_3$GeneID)),]
nonoverlap_edge_dexvdex_con[order(rank(nonoverlap_edge_dexvdex_con$PValue,ties.method = 'first')),]
nonoverlap_edge_dexvdex_con_name = get_gene_annotation(rownames(nonoverlap_edge_dexvdex_con[order(rank(nonoverlap_edge_dexvdex_con$PValue,ties.method = 'first')),])[1:20])

a = result_dexvl_dexl_3$GeneID
b = TreatTime_dexv.low3_vs_dex.low3$X
c = rownames(edgeR_select_dexvl_dexl)

GLM_uni = setdiff(a, union(b,c))
Des_uni = setdiff(b, union(a,c))
Edge_uni = setdiff(c, union(a,b))
GLM_uni_p = rep(0,length(GLM_uni))
for(i in 1:length(GLM_uni))
{
  GLM_uni_p[i] = result_dexvl_dexl_3[which(result_dexvl_dexl_3$GeneID== GLM_uni[i]),3]
}
a_1 = get_gene_annotation(GLM_uni[order(rank(GLM_uni_p,ties.method = 'first'))])

Des_uni_p = rep(0,length(Des_uni))
for(i in 1:length(Des_uni))
{
  Des_uni_p[i] = TreatTime_dexv.low3_vs_dex.low3[which(TreatTime_dexv.low3_vs_dex.low3$X== Des_uni[i]),6]
}
b_1 = get_gene_annotation(Des_uni[order(rank(Des_uni_p,ties.method = 'first'))])


Edge_uni_p = rep(0,length(Edge_uni))
for(i in 1:length(Edge_uni))
{
  Edge_uni_p[i] = edgeR_select_dexvl_dexl[which(rownames(edgeR_select_dexvl_dexl)== Edge_uni[i]),3]
}
c_1 = get_gene_annotation(Edge_uni[order(rank(Edge_uni_p,ties.method = 'first'))])




#####treatment comparison: cort vehcile high 3 vs cort high 3
result_cortvh_corth = matrix(0,20052,2)
#DEHOGT
K <- matrix(c(1,rep(0,3),-1,rep(0,12)),1)
for(i in 1:20052)
{
  genewise_data = counts(dds)[i,]
  
  genewise_dataframe = data.frame(read = genewise_data, treatment = X_factor)
  a = glm.nb(read/treatmentwise_factor_2 ~ treatment - 1, data = genewise_dataframe)
  a1 =  summary(glht(a, linfct = K)) 
  result_cortvh_corth[i,] = c(coef(a1), a1$test$pvalues[1])
}

p = p.adjust(result_cortvh_corth[,2], method = "BH")

result_cortvh_corth_1 = data.frame(GeneID = rownames(counts(dds)), log2fold = log(exp(1),base = 2)*result_cortvh_corth[,1], pvalue = p)

result_cortvh_corth_2 = filter(result_cortvh_corth_1,pvalue<0.05)

result_cortvh_corth_3 = filter(result_cortvh_corth_2,abs(log2fold)>1.5)

#obtain and store the names for the top 30 DE genes selected by proposed method
c_1 = result_cortvh_corth_3[order(result_cortvh_corth_3[,3])[1:30],]
a = get_gene_annotation(c_1[,1])
b = pmatch(c_1[,1],a$ENSEMBL)
d = c_1[b[!is.na(b)],]
d$SYMBOL = a$SYMBOL
write.csv(d,"cortvechigh_vs_corthigh_time3.csv", row.names = FALSE)

#load the results of DE genes from DESeq2 method
TreatTime_cort.high3_Vs_cort.vehicle.high3 <- read.csv("cort.high3_Vs_cort.vehicle.high3.csv")


# edgeR method
edgedata = DGEList(counts = counts(dds), group = X_factor)
edgedata = calcNormFactors(edgedata)
edgedata = estimateCommonDisp(edgedata, verbose=T)
edgedata = estimateTagwiseDisp(edgedata)
edge_cortvh_corth = exactTest(edgedata, pair=c(1, 5))
log2fold2 = edge_cortvh_corth$table$logFC

edge_cortvh_corth$table$PValue = p.adjust(edge_cortvh_corth$table$PValue, method="BH")

edgeR_select_cortvh_corth = filter(edge_cortvh_corth$table,abs(logFC)>1.5&PValue<0.05)

#draw the VennDiagram plot for figure 18
venn.diagram(
  x = list(result_cortvh_corth_3$GeneID, TreatTime_cort.high3_Vs_cort.vehicle.high3$X,rownames(edgeR_select_cortvh_corth)),
  category.names = c('GLM (NB)' ,"DESeq" , "EdgeR"),
  filename = 'cortvh_corth_3.png',
  output = TRUE ,
  imagetype="png" ,
  height = 380 , 
  width = 380 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  col=c("red", 'blue', '#00CC00'),
  fill = c(alpha("red",0.7), alpha('blue',0.7), alpha('#00CC00',0.7)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', '#00CC00'),
  rotation = 1,
  main = "Cort vehicle high 3 vs Cort high 3 ",
  main.cex = 0.5,
  sub.cex = 0.5 
)

#draw the p-value rank plot in figure 24
overlap_gene_2 = intersect(intersect(result_cortvh_corth_3$GeneID,TreatTime_cort.high3_Vs_cort.vehicle.high3$X),rownames(edgeR_select_cortvh_corth))

a_1 = get_gene_annotation(overlap_gene_2)
overlap_gene_2 = a_1$ENSEMBL
overlap_gene_reorder_2 = overlap_gene_2[order(result_cortvh_corth_3[which(result_cortvh_corth_3$GeneID%in%overlap_gene_2),3])]
overlap_gene_reorder_2 = overlap_gene_reorder_2[1:30]

glm_select_2 = result_cortvh_corth_3[which(result_cortvh_corth_3$GeneID%in%overlap_gene_reorder_2),c(1,3)]
overlap_gene_reorder_2 = glm_select_2[order(glm_select_2[,2]),1]
#glm_order_2 = match(overlap_gene_reorder_2,glm_select_2[,1])
glm_rank_2 = rank(glm_select_2[order(glm_select_2[,2]),2],ties.method = 'first')
#gene_name_2 = a_1$SYMBOL[glm_order_2]
gene_name_2 = get_gene_annotation(overlap_gene_reorder_2)

deq_select_2 = TreatTime_cort.high3_Vs_cort.vehicle.high3[which(TreatTime_cort.high3_Vs_cort.vehicle.high3$X%in%overlap_gene_reorder_2),c(1,6)]
deq_order_2 = match(overlap_gene_reorder_2,deq_select_2[,1])
deq_rank_2 = rank(deq_select_2[deq_order_2,2],ties.method = 'first')

edgeR_select_cortvh_corth$names = rownames(edgeR_select_cortvh_corth)
edge_select_cortvh_corth = edgeR_select_cortvh_corth[which(edgeR_select_cortvh_corth$names%in%overlap_gene_reorder_2),c(3,4)]
edge_order_2 = match(overlap_gene_reorder_2,edge_select_cortvh_corth[,2])
edge_rank_2 = rank(edge_select_cortvh_corth[edge_order_2,1],ties.method = 'first')

rank_plot_data_2 <- data.frame(ID = rep(overlap_gene_reorder_2,3),rank = c(glm_rank_2,deq_rank_2,edge_rank_2),method = c(rep("GLM (NB)",30),rep("DESeq",30),rep("EdgeR",30)))
rank_plot_data_2$num_ID = rep(-seq(-60,-1,2),3)
ggplot(rank_plot_data_2, aes(fill=method, y=rank, num_ID)) + 
  geom_bar(position="dodge", stat="identity",width=1.4) + 
  scale_x_continuous(breaks=-seq(-60,-1,2),labels=gene_name_2$SYMBOL) + coord_flip() +
  labs(title = "Cort vehcile high 3 vs Cort high 3", y = "Rank of adjusted pvalue", x = "") + 
  theme(axis.text.x = element_text(size = 0.1)) + theme_bw() +
  scale_fill_manual(values=c("DESeq"="blue", "EdgeR"="green", "GLM (NB)"="red"))


#####check the genes not in the overlap 
a = result_cortvh_corth_3$GeneID
b = TreatTime_cort.high3_Vs_cort.vehicle.high3$X
c = rownames(edgeR_select_cortvh_corth)

GLM_uni = setdiff(a, union(b,c))
Des_uni = setdiff(b, union(a,c))
Edge_uni = setdiff(c, union(a,b))
GLM_uni_p = rep(0,length(GLM_uni))
for(i in 1:length(GLM_uni))
{
  GLM_uni_p[i] = result_cortvh_corth_3[which(result_cortvh_corth_3$GeneID== GLM_uni[i]),3]
}
a_1 = get_gene_annotation(GLM_uni[order(rank(GLM_uni_p,ties.method = 'first'))])

Des_uni_p = rep(0,length(Des_uni))
for(i in 1:length(Des_uni))
{
  Des_uni_p[i] = TreatTime_cort.high3_Vs_cort.vehicle.high3[which(TreatTime_cort.high3_Vs_cort.vehicle.high3$X== Des_uni[i]),6]
}
b_1 = get_gene_annotation(Des_uni[order(rank(Des_uni_p,ties.method = 'first'))])


Edge_uni_p = rep(0,length(Edge_uni))
for(i in 1:length(Edge_uni))
{
  Edge_uni_p[i] = edgeR_select_cortvh_corth[which(rownames(edgeR_select_cortvh_corth)== Edge_uni[i]),3]
}
c_1 = get_gene_annotation(Edge_uni[order(rank(Edge_uni_p,ties.method = 'first'))])





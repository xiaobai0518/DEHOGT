# Differentially Expressed Heterogeneous Overdispersion Genes Testing for Count Data
## Authors: 

Yubai Yuan#1, Qi Xu#3, Wani Agaz#2, Jan Dahrendorff#2, Chengqi Wang#2, Janelle Donglasan#2, Sarah Burgan#2, Zachary Graham#2, Monica Uddin#2,
Derek E. Wildman#2, Annie Qu#3 
            
1: Department of Statistics, The Pennsylvania State University   
2: Genomics Program, College of Public Health, University of South Florida \
3: Department of Statistics, University of California, Irvine

## Description of Repository

This github repository includes the necessary dataset and codes for generate the simualation results and real application results in paper
**Differentially Expressed Heterogeneous Overdispersion Genes Testing for Count Data**. 


## How to run the codes

Software: R version 4.1.2 \
Please install the following R packages:

1. DESeq2
2. DESeqAnalysis
3. openxlsx
4. tibble
5. data.table
6. dplyr
7. edgeR
8. MASS
9. Rmpfr
10. VennDiagram
11. RColorBrewer
12. scales
13. ggplot2
14. multcomp

To run the simulation codes, first load all functions in **Read_count_generation.R and DE_gene_detection_methods.R**.

### Simulation 1

Run the codes in file **simulation_setting_1.R** 

### Simulation 2

Run the codes in file **simulation_setting_2.R** 

### Experiments on Microglia RNA-seq read count data

Load function in file **get_gene_annotation.R**.\
Load data under **Microglia RNA-seq read count/DE_genes_from_DESeq2**.\
Load **RNA_Read_Counts.csv** under **Microglia RNA-seq read count**.\
Load **Sample Information.csv** under **Microglia RNA-seq read count**.

Run the code in file **Microglia RNA experiment code/microglia_cell_rna_experiment.R**.\
The outputs includes the seclected genes from DESeq2[1], edgeR[2], and DEHOGT, and Figures 12 to 24. \
The experiments compare DE genes from DESeq2, edgeR, and DEHOGT under 7 different treatment pairs:

1. dex high 3 vs control
2. dex high 6 vs control
3. cort high 3 vs control 
4. cort high 6 vs control
5. dex vehicle high 3 vs dex high 3
6. dex vehcile low 3 vs dex low 3
7. cort vehcile high 3 vs cort high 3

Reference: \
[1]: Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. 
     Genome biology, 15(12), 1-21 \
[2]: Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of 
     digital gene expression data. bioinformatics, 26(1), 139-140.     
     
### Contact:
Dr. Yubai Yuan: yvy5509@psu.edu

















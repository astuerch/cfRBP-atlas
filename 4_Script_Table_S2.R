##### Table S2: determine dependencies of all miRNAs for miNRA-biogenesis RBPs
library(dplyr)
library(DESeq2)

#use only miRNA reads for this analysis
####assign mirna dependencies to RBPs

ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/exceRpt_miRNA_ReadCounts_allmiRNA_ordered_cleaned_batchadj.txt")
colData_batch <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/colData_batch_allmiRNA.txt")
colData_batch$genotype <- factor(colData_batch$genotype)
colData_batch$origin <- factor(colData_batch$origin)

#run deseq for WT vs KO
#select specific KO from data
colData_cell <- subset(colData_batch, genotype %in% c('WT', 'RBP_1') & origin == 'Cell')
colData_cm <- subset(colData_batch, genotype %in% c('WT', 'RBP_1') & origin == 'CM')

columns_to_keep_cell <- c("Gene",colData_cell$Sample)
columns_to_keep_cm <- c("Gene",colData_cm$Sample)

batch_cell <- dplyr::select(ex_mapped, all_of(columns_to_keep_cell))
batch_cm <-  dplyr::select(ex_mapped, all_of(columns_to_keep_cm))

#prepare data for deseq
data4pca_cell <- round(batch_cell[,-c(1)] )
data4pca_cm <- round(batch_cm[,-c(1)] )

rownames(data4pca_cell) <- batch_cell$Gene
rownames(data4pca_cm) <- batch_cm$Gene

#run deseq2
dds_cell <- DESeqDataSetFromMatrix(countData = data4pca_cell,
                                   colData = colData_cell,
                                   design = ~ genotype  )

dds_cm <- DESeqDataSetFromMatrix(countData = data4pca_cm,
                                 colData = colData_cm,
                                 design = ~ genotype  )


dds_cell <- estimateSizeFactors(dds_cell)
dds_cm <- estimateSizeFactors(dds_cm)


#normalize data
ddsr_cell <- DESeq(dds_cell)
ddsr_cm <- DESeq(dds_cm)

#results
res_cell <- results(ddsr_cell, contrast=c("genotype", "RBP_1", "WT"))
res_cm <- results(ddsr_cm, contrast=c("genotype", "RBP_1", "WT"))

table_cell <- as.data.frame(res_cell)
table_cell <- cbind(gene=rownames(table_cell), table_cell)

table_cm <- as.data.frame(res_cm)
table_cm <- cbind(gene=rownames(table_cm), table_cm)

#save the data for each RBP
RBP_1_cell <- table_cell[,c(1,3,6,7)]
colnames(RBP_1_cell)[2] <- "RBP_1_log2FoldChange"

RBP_1_cm <- table_cm[,c(1,3,6,7)]
colnames(RBP_1_cm)[2] <- "RBP_1_log2FoldChange"

#join all the statistical data to define dependencies
df_list_cell <- list(RBP_1_cell, DGCR8_cell, XPO5_cell, RBP_1_cell, TRBP_cell, PACT_cell, PACT.TRBP_cell, RBP_1_cell, AGO1_cell,AGO2_cell,AGO3_cell,AGO12_cell,AGO23_cell,AGO123_cell,AGO1234_cell)

# Perform full join on list of data frames
miRNA_cell <- df_list_cell[[1]]
for (i in 2:length(df_list_cell)) {
  miRNA_cell <- full_join(miRNA_cell, df_list_cell[[i]], by = "gene")
}

write.table(miRNA_cell, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/miRNA_RBP_dependency_cell.txt", quote = F, row.names = F, sep = "\t")



df_list_cm <- list(RBP_1_cm, DGCR8_cm, XPO5_cm, RBP_1_cm, TRBP_cm, PACT_cm, PACT.TRBP_cm, RBP_1_cm, AGO1_cm,AGO2_cm,AGO3_cm,AGO12_cm,AGO23_cm,AGO123_cm,AGO1234_cm)

# Perform full join on list of data frames
miRNA_cm <- df_list_cm[[1]]
for (i in 2:length(df_list_cm)) {
  miRNA_cm <- full_join(miRNA_cm, df_list_cm[[i]], by = "gene")
}

write.table(miRNA_cm, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/miRNA_RBP_dependency_cm.txt", quote = F, row.names = F, sep = "\t")


#assign dependencies with pvalue and log2FC (log2FC < -1, FDR < 0.05)
#### .------> miRNA_table_dependency_2fc_final.xlsx

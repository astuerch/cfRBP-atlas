#### FIGURE S3B: correlation matrix among RBP KOs
#instCM.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)
library(tibble)

#upload data
table1 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/all_readcounts_final_batch.txt")
coldata_batch1 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/colData_batch_all_batches.txt")

#select CELL or CM data
colData_CM <- subset(coldata_batch1, origin == 'Cell')
columns_to_keep_CM <- c("Gene",colData_CM$Sample)

batch_CM <- dplyr::select(table1, all_of(columns_to_keep_CM))

#prepare data for normalization
batch1 <- lapply(batch_CM[-1], function(x) round(as.numeric(x) + 0.1))
batch1 <- as.data.frame(batch1)
batch1[is.na(batch1)] <- 0

#filter for rows with minimum 5 reads
keep1 <- rowMeans(batch1) >= 50
batch2 <- batch1[keep1,]

#normalize data
dds_batch <- DESeqDataSetFromMatrix(countData = batch2,
                                    colData = colData_CM ,
                                    design = ~ genotype  )

#vst <- vst(dds_batch, blind=FALSE)
rlog <- rlog(dds_batch, blind=FALSE)

# Calculate correlation matrix for all replicates
correlation_matrix <- cor(assay(rlog))

# Create the heatmap
pheatmap(correlation_matrix,
         color = colorRampPalette(brewer.pal(n = 7, name =  "Reds"))(100),
         main = "Correlation Heatmap",
         fontsize = 6,  # Adjust the font size if needed
         cellwidth = 7, cellheight = 7,  # Adjust CM size if needed
         cluster_rows = T, cluster_cols = T)  # Disable row and column clustering

#13x13
#corr_map_all_DEGs


###### make analysis with average counts per genotype (combine replicates for each KO)

# Transpose the data so that samples become rows
transposed_data <- as.data.frame(t(assay(rlog))) # Exclude the first column (assumed to be identifiers)

# Add row names as a column (to be used for merging)
transposed_data$genotype <- colData_CM$genotype  #be sure that the order is the same

# Calculate the mean for each combination of Gene and KO
averaged_data <- transposed_data %>%
  group_by(genotype) %>%
  summarize(across(everything(), mean))

rownames(averaged_data) <- averaged_data$genotype

# Calculate correlation matrix
data_cor <-  t(averaged_data[,-1])

colnames(data_cor) <- averaged_data$genotype
correlation_matrix2 <- cor(data_cor)

# Create the heatmap
pheatmap(correlation_matrix2,
         main = "Correlation Heatmap - Cell samples",
         color = colorRampPalette(brewer.pal(n = 7, name =  "Reds"))(100),
         fontsize = 9,  # Adjust the font size if needed
         cellwidth = 10, cellheight = 10,  # Adjust CM size if needed
         cluster_rows = T, cluster_cols = T)  # Disable row and column clustering

#corr_map_all_all-DEGs_average
#10x10



### FIGURE S3E AND TABLE S1: DESEQ2 ANALYSIS of RBP KO CELL AND CM - TO PERFORM WITH RAW COUNTS AND NOT NORMALIZED DATA!!

# Prepare data before analysis – TO DO ONLY ONCE!!!
#1. remove double annotation miRNAs (remove gencode annotation and keep mirBase annotation)
#2. adjust raw counts among batches
library(sva)

ex_mapped <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/RBP_KO_batch/exceRpt_all_ReadCounts_allRBPs_ordered_cleaned.txt")
colData_batch <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/RBP_KO_batch/colData_batch_allRBPs.txt")
colData_batch$genotype <- factor(colData_batch$genotype)
colData_batch$origin <- factor(colData_batch$origin)

rownames(ex_mapped) <- ex_mapped$Gene
ex_mapped[is.na(ex_mapped)] <- 0

#batch adjusted
adjusted_counts <- ComBat_seq(as.matrix(ex_mapped[,-1]), batch=colData_batch$batch)
write.table(as.data.frame(adjusted_counts), "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/exceRpt_all_ReadCounts_allRBPs_ordered_cleaned_batchgenoriadj.txt", quote = F, row.names = T, sep = "\t")


###
library(dplyr)
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

#CHOOSE DATA AND METADATA FILES: 1) FOR FINAL SELECTED RBP KOs, 2) FOR INITIAL SCREENING RBP KOs THAT WERE NOT SELECTED FOR FOLLOWUP
#1)
ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/exceRpt_all_ReadCounts_allRBPs_ordered_cleaned_batchgenoriadj.txt") #final selected RBP KOs
colData_batch <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/colData_batch_allRBPs.txt")  #final selected RBP KOs
#2)
ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/exceRpt_all_ReadCounts_merge_paper_cleaned_batchadj.txt") #final selected RBP KOs
colData_batch <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/colData_batch_SCREEN_BATCH.txt")  #final selected RBP KOs


colData_batch$genotype <- factor(colData_batch$genotype)
colData_batch$origin <- factor(colData_batch$origin)

#select specific KO from data to perform deseq
colData_cell <- subset(colData_batch, genotype %in% c('WT', 'RBP_1') & origin == 'Cell')
colData_cm <- subset(colData_batch, genotype %in% c('WT', 'RBP_1') & origin == 'CM')

columns_to_keep_cell <- c("Gene",colData_cell$Sample)
columns_to_keep_cm <- c("Gene",colData_cm$Sample)

batch_cell <- dplyr::select(ex_mapped, all_of(columns_to_keep_cell))
batch_cm <-  dplyr::select(ex_mapped, all_of(columns_to_keep_cm))

#prepare data for deseq
data4pca_cell <- batch_cell[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table
data4pca_cm <- batch_cm[,-c(1)] + 0.1
data4pca_cell[is.na(data4pca_cell)] <- 0.1
data4pca_cm[is.na(data4pca_cm)] <- 0.1

data4pca_cell <- round(data4pca_cell)
data4pca_cm <- round(data4pca_cm)

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


#filter for rows with minimum 10 reads
keep1 <- rowMeans(counts(ddsr_cell, normalized = TRUE)) >= 10
keep2 <- rowMeans(counts(ddsr_cm, normalized = TRUE)) >= 10
ddsr_cell1 <- ddsr_cell[keep1,]
ddsr_cm1 <- ddsr_cm[keep2,]

#results
res_cell <- results(ddsr_cell1, contrast=c("genotype", "RBP_1", "WT"))
res_cm <- results(ddsr_cm1, contrast=c("genotype", "RBP_1", "WT"))

table_cell <- as.data.frame(res_cell)
table_cell <- cbind(gene=rownames(table_cell), table_cell)

table_cm <- as.data.frame(res_cm)
table_cm <- cbind(gene=rownames(table_cm), table_cm)


#how many adj. pvalues were less thant 0.1?
sum(res_cm$padj < 0.05, na.rm=TRUE) #893
sum(res_cell$padj < 0.05, na.rm=TRUE) #893


# Select the top differentially expressed genes
top_genes_cell <- subset(res_cell, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.8)
top_genes_cm <- subset(res_cm, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.8 )

# Get the count data for the top differentially expressed genes and the full dataset
counts_cell <- data.frame(cbind(gene = row.names(counts(ddsr_cell1)),counts(ddsr_cell1, normalized = TRUE)))
counts_cell1 <- as.data.frame(counts_cell[counts_cell$gene %in% rownames(top_genes_cell), ])
final_counts_cell <- inner_join(table_cell, counts_cell1, by="gene" )
final_tot_counts_cell <- inner_join(table_cell, counts_cell, by="gene" )


counts_cm <- data.frame(cbind(gene = row.names(counts(ddsr_cm1)),counts(ddsr_cm1, normalized = TRUE)))
counts_cm1 <- as.data.frame(counts_cm[counts_cm$gene %in% rownames(top_genes_cm), ])
final_counts_cm <- inner_join(table_cm, counts_cm1, by="gene" )
final_tot_counts_cm <- inner_join(table_cm, counts_cm, by="gene" )


write.table(final_counts_cell, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/deseq_KO_results/diff_Cell_RBP_1_norm_counts_filtered.txt", quote = F, row.names = F, sep = "\t")
write.table(final_counts_cm, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig3/deseq_KO_results/diff_CM_RBP_1_norm_counts_filtered.txt", quote = F, row.names = F, sep = "\t")


########## VENN DIAGRAMS

# combine deseq results from cell and cm
volcano_cell <- table_cell
volcano_cm <- table_cm

volcano_cell1 <- cbind(volcano_cell, Sample = "Cell" )
volcano_cm1 <- cbind(volcano_cm, Sample = "CM")
volcano <- rbind(volcano_cell1,volcano_cm1)

#select only DE genes (total genes, up and down regulated)
pass_thresholds_cell <- subset(volcano$gene, abs(volcano$log2FoldChange) > 0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "Cell")
pass_thresholds_cm <- subset(volcano$gene, abs(volcano$log2FoldChange) > 0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "CM")


venn.diagram(list(pass_thresholds_cell, pass_thresholds_cm),
             category.names = c("Cell" , "CM "),
             fill = c("sienna1", "darkturquoise"),
             filename = '~/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/RBP_KO_batch/venn_diagram_DEGs/venn_diagramm_RBP_1.png',
             main = "RBP_1 KO",
             main.cex = 5,
             main.fontface = "bold",
             output=TRUE,
             imagetype="png" ,
             height = 5000 ,
             width = 5000  ,
             #inverted = TRUE,
             #ext.text = TRUE,
             resolution = 300,
             compression = "lzw",
             lwd = 2,
             # Numbers
             cex = 4,
             fontface = "italic",
             cat.cex = 4,
             cat.fontface = "bold",
             #cat.default.pos = "outer",
             #cat.dist = c(-0.45, 0.45) ,
             cat.pos = c(1, 1) ,#
             fontfamily = "sans",
             alpha = c(0.7, 0.7),  scaled = TRUE, rotation.degree = 180)



#select only up or down DE genes
up_thresholds_cell <- subset(volcano$gene, volcano$log2FoldChange > 0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "Cell") #upregulated in RBP KO
up_thresholds_cm <- subset(volcano$gene, volcano$log2FoldChange > 0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "CM")

down_thresholds_cell <- subset(volcano$gene, volcano$log2FoldChange < -0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "Cell") #downregulated in RBP KO
down_thresholds_cm <- subset(volcano$gene, volcano$log2FoldChange < -0.8 & -log10(volcano$padj) > 1.3 & volcano$Sample == "CM")

venn.diagram(list(up_thresholds_cell, up_thresholds_cm),
             category.names = c("Cell" , "CM "),
             fill = c("sienna1", "turquoise1"),
             filename = '~/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/RBP_KO_batch/venn_diagram_DEGs/venn_diagramm_RBP_1_upregulated.png',
             main = "RBP_1 KO upregulated genes",
             main.cex = 5,
             main.fontface = "bold",
             output=TRUE,
             imagetype="png" ,
             height = 5000 ,
             width = 5000  ,
             #inverted = TRUE,
             #ext.text = TRUE,
             resolution = 300,
             compression = "lzw",
             lwd = 2,
             # Numbers
             cex = 4,
             fontface = "italic",
             cat.cex = 4,
             cat.fontface = "bold",
             #cat.default.pos = "outer",
             #cat.dist = c(-0.45, 0.45) ,
             cat.pos = c(1, 1) ,#
             fontfamily = "sans",
             alpha = c(0.5, 0.5),  scaled = TRUE, rotation.degree = 180)


venn.diagram(list(down_thresholds_cell, down_thresholds_cm),
             category.names = c("Cell" , "CM "),
             fill = c("sienna3", "turquoise3"),
             filename = '~/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/RBP_KO_batch/venn_diagram_DEGs/venn_diagramm_RBP_1_downregulated.png',
             main = "RBP_1 KO downregulated genes",
             main.cex = 5,
             main.fontface = "bold",
             output=TRUE,
             imagetype="png" ,
             height = 5000 ,
             width = 5000  ,
             #inverted = TRUE,
             #ext.text = TRUE,
             resolution = 300,
             compression = "lzw",
             lwd = 2,
             # Numbers
             cex = 4,
             fontface = "italic",
             cat.cex = 4,
             cat.fontface = "bold",
             #cat.default.pos = "outer",
             #cat.dist = c(-0.45, 0.45) ,
             cat.pos = c(1, 1) ,#
             fontfamily = "sans",
             alpha = c(0.7, 0.7),  scaled = TRUE, rotation.degree = 180)


intersect_up <- intersect(up_thresholds_cell, up_thresholds_cm)
intersect_down <- intersect(down_thresholds_cell, down_thresholds_cm)





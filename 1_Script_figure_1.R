# Prepare data before analysis – TO DO ONLY ONCE!!!
#1. remove double annotation miRNAs (remove gencode annotation and keep mirBase annotation)
#2. adjust raw counts among batches
library(sva)

ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned.txt")
colData_batch2 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/colData_WT.txt")
colData_batch2$genotype <- factor(colData_batch2$genotype)
colData_batch2$origin <- factor(colData_batch2$origin)

rownames(ex_mapped) <- ex_mapped$Gene

#batch adjusted
adjusted_counts <- ComBat_seq(as.matrix(ex_mapped[,-1]), batch=colData_batch2$batch)
write.table(as.data.frame(adjusted_counts), "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt", quote = F, row.names = T, sep = "\t")



######### FIGURE 1D: differential analysis between WT cell and CM for all biotypes
library(dplyr)
library(DESeq2)
library(ggrepel)
library(ggplot2)


ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt")
colData_batch2 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/colData_WT.txt")
colData_batch2$genotype <- factor(colData_batch2$genotype)
colData_batch2$origin <- factor(colData_batch2$origin)

rownames(ex_mapped) <- ex_mapped$Gene

#run deseq2
data4pca2 <- round(ex_mapped[,-1])
dds_2 <- DESeqDataSetFromMatrix(countData = data4pca2,
                              colData = colData_batch2,
                              design = ~ origin  )

dds_2 <- estimateSizeFactors(dds_2)

#normalize data
ddsr_2 <- DESeq(dds_2)

#filter for rows with minimum 10 reads
keep2 <- rowMeans(counts(ddsr_2, normalized = TRUE)) >= 10
ddsr_3 <- ddsr_2[keep2,]

#results
res_cell <- results(ddsr_3, contrast=c("origin", "Cell", "CM"))
res_cm <- results(ddsr_3, contrast=c("origin", "CM", "Cell"))

table_cell <- as.data.frame(res_cell)
table_cell <- cbind(gene=rownames(table_cell), table_cell)

table_cm <- as.data.frame(res_cm)
table_cm <- cbind(gene=rownames(table_cm), table_cm)


#how many adj. pvalues were less thant 0.5?
sum(res_cm$padj < 0.05, na.rm=TRUE) #893


# Select the top differentially expressed genes
top_genes_cell <- subset(res_cell, -log10(padj) > 1.3 & log2FoldChange > 0.8)
top_genes_cm <- subset(res_cm, -log10(padj) > 1.3 & log2FoldChange > 0.8 )

# Get the count data for the top differentially expressed genes
counts2 <- data.frame(cbind(gene = row.names(counts(ddsr_3)),counts(ddsr_3, normalized = TRUE)))

count_cell <- as.data.frame(counts2[counts2$gene %in% rownames(top_genes_cell), ])
final_counts_cell <- inner_join(table_cell, count_cell, by="gene" )

count_cm <- as.data.frame(counts2[counts2$gene %in% rownames(top_genes_cm), ])
final_counts_cm <- inner_join(table_cm, count_cm, by="gene" )

final_counts_cm2 <- inner_join(table_cm, counts2, by="gene" )


write.table(final_counts_cell, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/WT_batch/diff_cell_all_counts_filtered.txt", quote = F, row.names = F, sep = "\t")
write.table(final_counts_cm, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/WT_batch/diff_cm_all_counts_filtered.txt", quote = F, row.names = F, sep = "\t")


#plot volcano
ggplot(table_cm, aes(x = log2FoldChange, y = -log10(padj),  stroke = 1.2)) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 0.8 & -log10(padj) > 1.3, ifelse(log2FoldChange > 0.8, "#62c3ef", "#e6a3c3"), "grey")),
             size = 1.5, alpha = 0.9)  +
  scale_color_identity() +
  labs(x = "log2 Fold Change", y = "-log10 FDR") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-0.8,0.8), linetype="dashed") +
  #scale_y_continuous(limits=c(0, 110)) +
  theme_classic(base_size = 20, base_family = "Helvetica") +
  labs( title = "All biotypes",
        #subtitle = "protein_coding",
        caption = "log2FC cutoff = 0.8, FDR cutoff = 0.05")


#WT_all_biotypes_batchadj
#10x10



######### FIGURE 1E: extraction of DEGs between WT cell and CM for all biotypes
library(dplyr)

diff_cell <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/diff_cell_all_counts_filtered.txt")
diff_CM <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/diff_cm_all_counts_filtered.txt")

split_values_cell <- strsplit(diff_cell$gene, ":")
diff_cell$gene_subgroup <- sapply(split_values_cell, function(x) x[2])
diff_cell$gene_name <- sapply(split_values_cell, function(x) x[1])

split_values_CM <- strsplit(diff_CM$gene, ":")
diff_CM$gene_subgroup <- sapply(split_values_CM, function(x) x[2])
diff_CM$gene_name <- sapply(split_values_CM, function(x) x[1])

deg_cell <- diff_cell %>%
  group_by(gene_subgroup) %>%
  summarise(unique_genes = n_distinct(gene_name))

deg_CM <- diff_CM %>%
  group_by(gene_subgroup) %>%
  summarise(unique_genes = n_distinct(gene_name))

#transfer numbers to stats_WT.prism


######### FIGURE 1F: selective release of exRNA or passive? check correlation in-out per biotype (use RPM data full libraries)
library(dplyr)
library(DESeq2)
library(ggrepel)
library(ggplot2)

# For small RNA has been used smallRNA seq for both intra and extracellular fractions.
# For protein coding RNA has been used longRNA seq for intracellular fraction, and smallRNA seq for extracellular fraction

### FOR SMALL RNA – to do for each biotype
data <-  read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt")

data[is.na(data)] <- 0

split_values1 <- strsplit(data$Gene, ":")
data$gene_subgroup <- sapply(split_values1, function(x) x[2]) #could give errors if you want to rerun the batch df generation
data$gene_name <- sapply(split_values1, function(x) x[1])

#filter for biotype
data_filt <- data[data$gene_subgroup == "miRNA",] ## CHANGE BIOTYPE
rownames(data_filt) <-  data_filt[,23]

#calculate rpm
table1_filtered <- data_filt[!rowSums(data_filt[,-c(1,22,23)]) == 0, ] #remove protein_codings not expressed
total_reads <- colSums(table1_filtered[,-c(1,22,23)])
rpm <- sweep(table1_filtered[,-c(1,22,23)], 2, total_reads, "/") * 1e6

rpm$WT_Cell <- rowMeans(rpm[, 1:10])
rpm$WT_CM <- rowMeans(rpm[, 11:20])

new_data <- rpm[,c(21:22)]

new_data$log_WT_Cell <- log(new_data$WT_Cell + 1)
new_data$log_WT_CM <- log(new_data$WT_CM + 1)

#calculate correlation coefficients
cor_coef1 <- cor(new_data$log_WT_Cell, new_data$log_WT_CM) #pearson test
r_squared1 <- cor_coef1^2
max_value <- max(c(max(new_data$log_WT_Cell), max(new_data$log_WT_CM)))


# Plot the data 
ggplot(new_data, aes(x = log_WT_Cell, y = log_WT_CM)) +
  geom_point() +  # Add the scatter points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  geom_ribbon(aes(ymin = log(0.5 * WT_Cell), ymax = log(2 * WT_Cell)), fill = "blue", alpha = 0.2) +
  labs(x = "log2(RPM+1) \n [mean of intracellular fraction]", y = "log2(RPM+1) \n [mean of extracellular fraction]", title = "miRNA") +
  coord_fixed(ratio = 1) +
  xlim(0, max_value) +  # Set x-axis limits
  ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  annotate("text", x = 2, y = 10, label = paste0("n = ", round(nrow(new_data), 2)),col = "blue", size = 6) +
  annotate("text", x = 2, y = 9, label = paste0("R = ", round(cor_coef1, 2)),col = "blue", size = 6)

#WT_miRNA
# 5x5


### FOR PROTEIN-CODING RNA, where longRNA seq has been used for intracellular fraction
data <-  read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt")
df_longRNA <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/Cell_longRNA_batch1_umidep_counts2_reordered.txt")

data[is.na(data)] <- 0

split_values1 <- strsplit(data$Gene, ":")
data$gene_subgroup <- sapply(split_values1, function(x) x[2]) #could give errors if you want to rerun the batch df generation
data$gene_name <- sapply(split_values1, function(x) x[1])

#filter for biotype
data_filt <- data[data$gene_subgroup == "protein_coding",]
rownames(data_filt) <-  data_filt[,23]

#for protein coding
data_filt <- inner_join(data_filt, df_longRNA[,c(1:4)], by = c("gene_name" = "Gene"))
data_filt[is.na(data_filt)] <- 0

#calculate rpm
table1_filtered <- data_filt[!rowSums(data_filt[,-c(1:11,22,23)]) == 0, ] #remove protein_codings not expressed
total_reads <- colSums(table1_filtered[,-c(1:11,22,23)])
rpm <- sweep(table1_filtered[,-c(1:11,22,23)], 2, total_reads, "/") * 1e6

rpm$WT_Cell <- rowMeans(rpm[, 11:13])
rpm$WT_CM <- rowMeans(rpm[, 1:10])

new_data <- rpm[,c(14:15)]

new_data$log_WT_Cell <- log(new_data$WT_Cell + 1)
new_data$log_WT_CM <- log(new_data$WT_CM + 1)

#calculate correlation coefficients
cor_coef1 <- cor(new_data$log_WT_Cell, new_data$log_WT_CM) #pearson test
r_squared1 <- cor_coef1^2
max_value <- max(c(max(new_data$log_WT_Cell), max(new_data$log_WT_CM)))


# Plot the data
ggplot(new_data, aes(x = log_WT_Cell, y = log_WT_CM)) +
  geom_point() +  # Add the scatter points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  geom_ribbon(aes(ymin = log(0.5 * WT_Cell), ymax = log(2 * WT_Cell)), fill = "blue", alpha = 0.2) +
  labs(x = "log2(RPM+1) \n [mean of intracellular fraction]", y = "log2(RPM+1) \n [mean of extracellular fraction]", title = "protein_coding") +
  coord_fixed(ratio = 1) +
  xlim(0, max_value) +  # Set x-axis limits
  ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  annotate("text", x = 7.5, y = 2, label = paste0("n = ", round(nrow(new_data), 2)),col = "blue", size = 6) +
  annotate("text", x = 7.5, y = 1, label = paste0("R = ", round(cor_coef1, 2)),col = "blue", size = 6)

#WT_protein_coding
# 5x5




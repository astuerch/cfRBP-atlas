##### FIGURE S4A and B: upset plot for miRNA dependencies
#generating upset plots to visualise shared DEGs among RBPs and find possible common pathways for exRNA biogenesis and export

#install.packages("UpSetR")
library(UpSetR)
library(dplyr)

#prepare list of DEGs for each RBP (simple file with RBP and dependent miRNAs)
#to repeat the same for Cell and CM data
dep_mirna <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/RBP_cm_dependency.txt")

# Grouping genes by KO and aggregating them into lists
grouped <- dep_mirna %>%
  group_by(RBP) %>%
  summarize(miRNA = list(unique(miRNA)))

# Create a list of gene sets for each KO
gene_sets <- grouped$miRNA

# Create an empty list to store the genes for each KO
gene_list <- list()

# Iterate through each gene set and create a list of genes for each KO
for (i in 1:length(gene_sets)) {
  gene_list[[grouped$RBP[i]]] <- gene_sets[[i]]
}

# Generate UpSet plot
upset(fromList(gene_list), sets.bar.color = "#f8ec1f",
      nsets = 50, order.by = c("freq"), nintersects = 100, cutoff = 10, point.size = 3.5, line.size = 1.2,
      mainbar.y.label = "miRNA Intersections in CM", sets.x.label = "RBP-miRNA dependency",
      text.scale = c(2, 2, 1.5, 1, 2, 2))


#upset_mirna_dep_cell
#8x13



##### FIGURE S4C: miRNA enrichment in cell or cm fractions
library(dplyr)
library(DESeq2)

#use only miRNA reads for this analysis
ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/exceRpt_miRNA_ReadCounts_allmiRNA_ordered_cleaned_batchadj.txt")
colData_batch <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/colData_batch_allmiRNA.txt")
colData_batch$genotype <- factor(colData_batch$genotype)
colData_batch$origin <- factor(colData_batch$origin)

wt_coldata <- subset(colData_batch, genotype %in% c('WT'))
columns_to_keep_wt <- c("Gene",wt_coldata$Sample)
batch_wt <- dplyr::select(ex_mapped, all_of(columns_to_keep_wt))
data4pca_wt <- round(batch_wt[,-1])
rownames(data4pca_wt) <- batch_wt$Gene

#run deseq2
dds_2 <- DESeqDataSetFromMatrix(countData = data4pca_wt,
                                colData = wt_coldata,
                                design = ~ origin  )

dds_2 <- estimateSizeFactors(dds_2)

#normalize data
ddsr_2 <- DESeq(dds_2)

#results
#res_cell <- results(dds_2, contrast=c("origin", "Cell", "CM"))
res_wt <- results(ddsr_2, contrast=c("origin", "CM", "Cell"))

#table_cell <- as.data.frame(res_cell)
#table_cell <- cbind(gene=rownames(table_cell), table_cell)

table_wt <- as.data.frame(res_wt)
table_wt <- cbind(gene=rownames(table_wt), table_wt)

counts_wt <- data.frame(cbind(gene = row.names(counts(ddsr_2)),counts(ddsr_2, normalized = TRUE)))
final_counts_wt <- inner_join(table_wt, counts_wt, by="gene" )

write.table(table_wt, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/miRNA_enrichment_WT_CM.txt", quote = F, row.names = F, sep = "\t")
write.table(final_counts_wt, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/miRNA_enrichment_counts_WT_CM.txt", quote = F, row.names = F, sep = "\t")

#### collect numbers and plot in prism (miRNA_dep_stats)




##### FIGURE S4E: miRNA expression in WT CM
# see prism file "chimeric_results.prism"



##### FIGURE S4F: miRNA-target correlations for top20 selected miRNAs
results_df_cm <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/correlation-results_cm_WT.txt")

results_df_cm <- results_df_cm %>%
  left_join(mean_cor_values, by = "miRNA")

ggplot(data = results_df_cm, aes(x = cor_value_cm)) +
  geom_density(aes(y = ..density.., fill = "CM"), alpha = 0.7) +
  scale_fill_manual(values = c("CM" = "#62c3ef")) +
  ggtitle("Comparison of Correlation Coefficients") +
  xlab("Correlation Coefficient") +
  ylab("Density") +
  theme_classic() +
  facet_wrap(~ miRNA, nrow = 7, ncol = 3, scales = "free", drop = TRUE) +
  geom_vline(aes(xintercept = mean_cor_value), color = "dodgerblue", linetype = "dashed", size = 1)
#top10_miRNA_correlation_distr_WT_CM
#10x6


##### FIGURE S4G: miRNA-target expression comparison in cell and cm of WT and AGO KO
# see prism file "chimeric_results.prism"


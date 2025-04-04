##### FIGURE 4A: heatmap of miRNA dependencies on RBPs
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

write.table(miRNA_cell, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_dependency_table/miRNA_RBP_dependency_cell.txt", quote = F, row.names = F, sep = "\t")



df_list_cm <- list(RBP_1_cm, DGCR8_cm, XPO5_cm, RBP_1_cm, TRBP_cm, PACT_cm, PACT.TRBP_cm, RBP_1_cm, AGO1_cm,AGO2_cm,AGO3_cm,AGO12_cm,AGO23_cm,AGO123_cm,AGO1234_cm)

# Perform full join on list of data frames
miRNA_cm <- df_list_cm[[1]]
for (i in 2:length(df_list_cm)) {
  miRNA_cm <- full_join(miRNA_cm, df_list_cm[[i]], by = "gene")
}

write.table(miRNA_cm, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_dependency_table/miRNA_RBP_dependency_cm.txt", quote = F, row.names = F, sep = "\t")






###### FIGURE 4B: scatterplot cell cm fraction mirna with dot color showing dependencies to RBPs
library(ggplot2)
library(dplyr)
library(readxl)

data <-  read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/deseq2_analysis/final_batch/smallRNA/WT_batch/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt")
label <- read_excel("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_dependency_table/miRNA_table_dependency_label.xlsx")

#calculate rpm
data[is.na(data)] <- 0
table1_filtered <- data[!rowSums(data[,-1]) == 0, ] #remove miRNAs not expressed

total_reads <- colSums(table1_filtered[,-1])
rpm <- sweep(table1_filtered[,-1], 2, total_reads, "/") * 1e6

split_values1 <- strsplit(table1_filtered$Gene, ":")
rpm$gene_subgroup <- sapply(split_values1, function(x) x[2]) #could give errors if you want to rerun the batch df generation
rpm$gene <- sapply(split_values1, function(x) x[1]) #could give errors if you want to rerun the batch df generation

rpm$WT_Cell <- rowMeans(rpm[, 1:10])
rpm$WT_CM <- rowMeans(rpm[, 11:20])

#filter for biotype
data_filt <- rpm[rpm$gene_subgroup == "miRNA",]

new_data <- data_filt[,c(22:24)]
new_data$log_WT_Cell <- log(new_data$WT_Cell + 1)
new_data$log_WT_CM <- log(new_data$WT_CM + 1)

plot <- inner_join(new_data, label, by="gene")

cor_coef1 <- cor(plot$log_WT_Cell, plot$log_WT_CM) #pearson test
r_squared1 <- cor_coef1^2
max_value <- max(c(max(plot$log_WT_Cell), max(plot$log_WT_CM)))

plot$enrichment <- factor(plot$enrichment,
                          levels = c("RBP_dependent", "RBP_independent", "inconclusive", "not expressed"))


# Plot the data with conditional coloring
ggplot(data = plot)+
  geom_point(aes(x = log_WT_Cell, y = log_WT_CM, fill = enrichment), alpha = 0.8, pch=21, size=2) +
  #geom_point(data = RBP_ko, aes(x = log_AGO_Cell, y = log_AGO_CM), color = "red", alpha = 0.5) +
  labs(x = "log2(RPM+1) \n [mean of intracellular fraction]", y = "log2(RPM+1) \n [mean of extracellular fraction]", title = "miRNA") +
  scale_fill_manual(values = c("#f8ec1a", "#ee3224","#fffef1",  "grey45" ), name = "RBP dependency") +  # Apply the custom color palette
  coord_fixed(ratio = 1) +
  xlim(0, max_value) +  # Set x-axis limits
  ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  annotate("text", x = 2, y = 10, label = paste0("n = ", round(nrow(plot), 2)),col = "black", size = 6) +
  annotate("text", x = 2, y = 9, label = paste0("R = ", round(cor_coef1, 2)),col = "black", size = 6)

#WT_miRNA_dependency
# 7x5




###### FIGURE 4D: scatterplot WT vs AGO1234 KO CM fraction to show decrease of protein_coding transcripts in mutant
###### check if AGO mutant has less mRNA in CM
data <-  read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig4/exceRpt_all_ReadCounts_AGO_ordered_cleaned_batchadj.txt")

data[is.na(data)] <- 0

#separate name and biotype
split_values1 <- strsplit(data$Gene, ":")
data$gene_subgroup <- sapply(split_values1, function(x) x[2]) #could give errors if you want to rerun the batch df generation
data$gene_name <- sapply(split_values1, function(x) x[1])

#filter for protein-coding RNA
data_filt <- data[data$gene_subgroup == "protein_coding",]
rownames(data_filt) <-  data_filt[,1]

#get only wt and ago1234 ko cm
data_filt_ago <- data_filt[,c(26:28, 47:49)]

#calculate rpm
table1_filtered <- data_filt_ago[!rowSums(data_filt_ago) == 0, ] #remove protein_codings not expressed
total_reads <- colSums(table1_filtered)
rpm <- sweep(table1_filtered, 2, total_reads, "/") * 1e6

#average replicates
rpm$WT <- rowMeans(rpm[, 1:3])
rpm$AGO_mut <- rowMeans(rpm[, 4:6])

new_data <- rpm[,c(7:8)]

new_data$log_WT <- log(new_data$WT + 1)
new_data$log_AGO_mut <- log(new_data$AGO_mut + 1)

#calculate correlation
cor_coef1 <- cor(new_data$log_WT, new_data$log_AGO_mut) #pearson test
r_squared1 <- cor_coef1^2
max_value <- max(c(max(new_data$log_WT), max(new_data$log_AGO_mut)))

# Plot the data with conditional coloring
ggplot(data = new_data)+
  geom_point(aes(x = log_WT, y = log_AGO_mut, color = log_WT - log_AGO_mut > 1), alpha = 0.8) + #
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "grey27")) +
  #geom_point(data = RBP_ko, aes(x = log_AGO_Cell, y = log_AGO_CM), color = "red", alpha = 0.5) +
  labs(x = "log2(RPM+1) \n [mRNA in WT CM]", y = "log2(RPM+1) \n [mRNA in AGO KO CM]", title = "protein_coding") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_fixed(ratio = 1) +
  xlim(0, max_value) +  # Set x-axis limits
  ylim(0, max_value) +  # Set y-axis limits
  theme_classic(base_size = 15, base_family = "Helvetica") +
  theme(legend.position = "none") +
  annotate("text", x = 2, y = 9.5, label = paste0("n = ", round(nrow(new_data), 2)),col = "black", size = 6) +
  annotate("text", x = 2, y = 9, label = paste0("R = ", round(cor_coef1, 2)),col = "black", size = 6)
  
#AGO_protein_coding
# 5x5



###### FIGURE 4E: correlation density distribution between miRNA and targets in CM vs randomized couples
library(data.table)

input <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/input_permutations.bed", header=FALSE)

regions1 <- read.delim("~/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/hg38_mrna_exons_no_miRNA_no_snRNA.bed", header=FALSE)
regions <- regions1 %>% select(V1,V2,V3) %>% rename(V1="chr", V2="start", V3="end")

#function to shuffle coordinates for permutations
shuffle_coords <- function(bed, regions) {
  set.seed(NULL)
  shuffled <- bed
  for (i in seq_len(nrow(bed))) {
    valid_regions <- regions[regions$chr == bed$chr[i],]
    chosen_region <- valid_regions[sample(nrow(valid_regions), 1),]
    shift_max <- chosen_region$end - bed$end[i]
    shift_min <- chosen_region$start - bed$start[i]
    shift <- sample(shift_min:shift_max, 1)
    shuffled$start[i] <- bed$start[i] + shift
    shuffled$end[i] <- bed$end[i] + shift
  }
  shuffled
}

#set.seed(123)  # for reproducibility
permuted_correlations2 <- data.frame()

# Using replicate to perform the operations 1000 times
permuted_correlations2 <- replicate(10, {
  correlations <- list()
  for (i in 1:100) {
    shuffled_bed <- shuffle_coords(input, regions)
    gr1 <- makeGRangesFromDataFrame(shuffled_bed, keep.extra.columns = TRUE)
    gr2 <- makeGRangesFromDataFrame(wt_cm1[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
    overlaps <- findOverlaps(gr1, gr2, maxgap = 5, ignore.strand = TRUE)
    df_gr1 <- as.data.frame(gr1[queryHits(overlaps)], row.names = NULL)
    df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)], row.names = NULL)
    shuffle_cm <- cbind(df_gr1, df_gr2)[, c(8:10, 12, 6, 7, 13)]
    
    if (nrow(shuffle_cm) > 2) {
      variance_miRNA <- var(shuffle_cm$coverage_miRNA)
      if (!is.na(variance_miRNA) && variance_miRNA != 0) {
        correlations[[i]] <- cor(shuffle_cm$coverage_miRNA, shuffle_cm$coverage)
      }
    } else {
      correlations[[i]] <- NA
    }
  }
  # Convert list of correlations to a dataframe to be returned by each replicate
  as.data.frame(unlist(correlations))
}, simplify = FALSE) %>% bind_rows() %>% drop_na()

#plot figure 4E
ggplot() +
  geom_density(data = permuted_correlations2, aes(x = `unlist(correlations)`, y = ..density.., fill = "Permutations"), alpha = 0.7) +
  geom_density(data = results_df_cm, aes(x = cor_value_cm, y = ..density.., fill = "CM"), alpha = 0.7) +
  scale_fill_manual(values = c("Permutations" = "grey96", "CM" = "#62c3ef")) +
  #geom_vline(aes(xintercept = mean(results_df_cell$cor_value_cell)), color = "deeppink", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = mean(results_df_cm$cor_value_cm)), color = "#62c3ef", linetype = "dashed", size = 1) +
  annotate("text", x = -1, y = max(density(results_df_cm$cor_value_cm)$y)-0.5, label = paste("Chi-squared test: \n p-value < 2.2e-16 "), hjust = 0, color = "grey6", size = 5) +
  ggtitle("1000 permutations") +
  xlab("Correlation Coefficient") +
  ylab("Density") +
  theme_classic()


#perform chi-square stats by diving the distribution in bins of 0.2
breaks <- pretty(results_df_cm$cor_value_cm, n = 10)
observed_counts <- cut(results_df_cm$cor_value_cm, breaks=breaks, right=TRUE, include.lowest=TRUE)
table_cm <- table(observed_counts)

breaks <- pretty(permuted_correlations3$`unlist(correlations)`, n = 10)
observed_counts <- cut(permuted_correlations3$`unlist(correlations)`, breaks=breaks, right=TRUE, include.lowest=TRUE)
table_perm <- table(observed_counts)

combined_matrix <- rbind(table_cm, table_perm)
# Chi-square test
chisq.test(combined_matrix)




###### FIGURE 4F: number of target fragments for top20 miRNAs, with related dependency to AGO KO in CM

#number of targets per miRNA in filtered data from all replicates
result3 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/No_targets_top20_dependencyAGO.txt", header=T)

# Ordering miRNAs by total count
result3 <- result3 %>%
  arrange(desc(unique_targets)) %>%
  mutate(miRNA = factor(miRNA, levels = unique(miRNA))) # Ensure factor levels follow the order

ggplot(result3, aes(x = miRNA, y = unique_targets, fill = subcategory)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(y = "Unique Targets", fill = "Subcategory") +
  scale_fill_manual(values = c("WT" = "#62c3ef", "AGO mutant" = "grey"),
                    breaks = c("WT", "AGO mutant")) +
  ggtitle("Fragments for top10 miRNAs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust text angle for better visibility

#fragments_top10miRNA_counts_CM
#3x5




###### FIGURE 4G: expression of target framgents between WT and AGO KO in cells
# see prism file "chimeric_results.prism"

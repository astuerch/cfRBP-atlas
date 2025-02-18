
### FIGURE S1A: correlation among replicates
data <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/exceRpt_all_ReadsPerMillion_ordered_cleaned.txt", header=T)
data[is.na(data)] <- 0

# Calculate correlation matrix replicates WT Cell and CM
correlation_matrix <- cor(data[,-1])

# Visualize as a heatmap using pheatmap
library(pheatmap)
library(RColorBrewer)


# Create the heatmap
# Set up color palette
my_palette <- colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100)

# Create correlation plot with values
pheatmap(correlation_matrix,
         color = my_palette,
         main = "Correlation Heatmap",
         fontsize = 6,
         cellwidth = 15, cellheight = 15,
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE,  # Display correlation values
         fontsize_number = 5,  # Adjust font size for correlation values
         fontsize_row = 6, fontsize_col = 6  # Adjust font size for row and column labels
)


### FIGURE S1B: fragment length distribution of different biotypes betwen WT Cell and CM
####### FIRST PART DONE IN EULER CLUSTER WITH BASH SCRIPTING

data <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/WT_all_read_length_summary.txt", quote="\"", comment.char="") #to do for each biotype


#normalization for biotypes
data$V1[data$V3 == "Cell"] <- (data$V1[data$V3 == "Cell"]/sum(data$V1[data$V3 == "Cell"]))  *1000000
data$V1[data$V3 == "CM"] <- (data$V1[data$V3 == "CM"]/sum(data$V1[data$V3 == "CM"]))  *1000000

#plot
ggplot(data, aes(x = V2, y = V1 , fill = V3)) +
  geom_area(alpha = 0.5, position = 'identity')+
  scale_fill_manual(values = c("#e6a3c3", "#62c3ef")) +  # Assigning fill colors
  #scale_x_continuous(breaks = unique(as.numeric(factor(data$V2))), labels = levels(factor(data$V2))) +
  scale_x_continuous(n.breaks = 14) +
  labs(title = "All biotypes",
       x = "Length / nt",
       y = "RPM") +
  theme_classic(base_family = "Arial")+
  theme(axis.text = element_text(size = 12),  # Adjust font size for axis text
        axis.title = element_text(size = 15),  # Adjust font size for axis titles
        plot.title = element_text(size = 24, face = "bold"))





### FIGURE S1D: heatmap DEGs betwen WT Cell and CM
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyverse)

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

# Select the top differentially expressed genes
top_genes_cell <- subset(res_cell, -log10(padj) > 1.3 & log2FoldChange > 0.8)
top_genes_cm <- subset(res_cm, -log10(padj) > 1.3 & log2FoldChange > 0.8 )

# Get the count data for the top differentially expressed genes
counts2 <- data.frame(cbind(gene = row.names(counts(ddsr_3)),counts(ddsr_3, normalized = TRUE)))

count_cell <- as.data.frame(counts2[counts2$gene %in% rownames(top_genes_cell), ])
count_cm <- as.data.frame(counts2[counts2$gene %in% rownames(top_genes_cm), ])


#extract DE gene names for cell and cm
rows_to_keep <- unique(union(count_cell[,1], count_cm[,1]))

#collect count data for selected genes
norm_counts <- as.data.frame(counts(ddsr_2, normalized=TRUE))
norm_counts <- cbind(gene=rownames(norm_counts), norm_counts)

rownames(norm_counts) <- norm_counts$gene

batch <- norm_counts[norm_counts$gene %in% rows_to_keep, ]
batch1 <-round(batch[,-1]+0.1)

#normalize data
dds_batch <- DESeqDataSetFromMatrix(countData = batch1,
                                    colData = colData_batch2,
                                    design = ~ origin  )

rld <- rlog(dds_batch, blind=FALSE)


#heatmap DE genes only
install.packages("viridis")
library(viridis)

out <- pheatmap(assay(rld),
                show_colnames = T,
                show_rownames = F,
                fontsize_col = 8,
                fontsize_row = 5,
                cluster_cols=F, cluster_rows=T,
                #cutree_cols = 7,
                #cutree_rows = 3,
                gaps_col = 10 ,
                scale="row",
                color = colorRampPalette(c("#151515", "#4A4A4A","#919191", "#C4C4C3", "#EFEFF2", "white","#C0A9C4", "#815EA1", "#724E8B", "#592F73", "#3C2250"))(100),#colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),
                #color = colorRampPalette(brewer.pal(11, "RdBu"))(100),
                border_color = F, treeheight_col=30,
                cellwidth = 10, cellheight = 0.2,
                clustering_distance_rows = "euclidean", clustering_method = "complete")

#10x10
#heatmap_diff_all_WT


### FIGURE S1E: volcano plots to see single genes differentially secreted <-  normalized by biotype counts
library(dplyr)
library(DESeq2)
library(ggrepel)
library(ggplot2)

data <-  read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/cell_cm_all_ReadCounts_ordered_cleaned_batchadj.txt")
data[is.na(data)] <- 0

colData_batch2 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig1/colData_WT.txt")
colData_batch2$genotype <- factor(colData_batch2$genotype)
colData_batch2$origin <- factor(colData_batch2$origin)


split_values1 <- strsplit(data$Gene, ":")
data$gene_subgroup <- sapply(split_values1, function(x) x[2]) #could give errors if you want to rerun the batch df generation
data$gene_name <- sapply(split_values1, function(x) x[1])

#filter for biotype
data_filt <- data[data$gene_subgroup == "protein_coding",] #CHANGE FOR EACH BIOTYPE
rownames(data_filt) <-  data_filt[,23]

#run deseq
dds_2 <- DESeqDataSetFromMatrix(countData = round(data_filt[,-c(1,22,23)]),
                                colData = colData_batch2,
                                design = ~ origin  )


dds_2 <- estimateSizeFactors(dds_2)

#normalize data
ddsr_2 <- DESeq(dds_2)

#filter for rows with minimum 10 reads
keep2 <- rowMeans(counts(ddsr_2, normalized = TRUE)) >= 10
ddsr_3 <- ddsr_2[keep2,]

#results
res_cm <- results(ddsr_3, contrast=c("origin", "CM", "Cell"))
table_cm <- as.data.frame(res_cm)
table_cm <- cbind(gene=rownames(table_cm), table_cm)

counts2 <- data.frame(cbind(gene = row.names(counts(ddsr_3)),counts(ddsr_3, normalized = TRUE)))

final_table <- full_join(table_cm, counts2, by = "gene")

#no of degs
subset(results(ddsr_3, contrast=c("origin", "Cell", "CM")), -log10(padj) > 1.3 & log2FoldChange > 0.8 )
subset(res_cm, -log10(padj) > 1.3 & log2FoldChange > 0.8 )



#plot volcano with gene names
ggplot(table_cm, aes(x = log2FoldChange, y = -log10(padj), label = gene , stroke = 1.2)) +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 0.8 & -log10(padj) > 1.3, ifelse(log2FoldChange > 0.8, "#62c3ef", "#e6a3c3"), "grey")),
             size = 1.5, alpha = 0.9)  +
  scale_color_identity() +
  labs(x = "log2 Fold Change", y = "-log10 FDR") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-0.8,0.8), linetype="dashed") +
  #scale_y_continuous(limits=c(0, 110)) +
  geom_text_repel(size = 3, segment.linetype = 1, segment.size = 0.2, min.segment.length = 0.1,  box.padding = 0.5, max.overlaps = 15) +
  theme_classic(base_size = 20, base_family = "Helvetica") +
  labs( title = "protein_coding",
        #subtitle = "protein_coding",
        caption = "log2FC cutoff = 0.8, FDR cutoff = 0.05")


#WT_protein_coding_batchadj
#10x10




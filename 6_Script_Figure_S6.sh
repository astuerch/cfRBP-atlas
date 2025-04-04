#### Figure S6B-C: differences between noPNK and PNK treatment in cell and cm
#to run in R
library(dplyr)
library(DESeq2)
library(sva)

#perform deseq2 for different rna biotypes
ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/exceRpt_all_ReadCounts.txt")
colData_batch2 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/colData_batch_IP_1.txt")
colData_batch2$genotype <- factor(colData_batch2$genotype)
colData_batch2$origin <- factor(colData_batch2$origin)
colData_batch2$treatment <- factor(colData_batch2$treatment)

rownames(ex_mapped) <- ex_mapped$Gene

#select specific data
colData_cell <- subset(colData_batch2, origin == "Cell" & genotype == "WT")
colData_cm <- subset(colData_batch2, origin == "CM" & genotype == "WT")

columns_to_keep_cell <- c("Gene",colData_cell$Sample)
columns_to_keep_cm <- c("Gene",colData_cm$Sample)

batch_cell <- dplyr::select(ex_mapped, all_of(columns_to_keep_cell))
batch_cm <-  dplyr::select(ex_mapped, all_of(columns_to_keep_cm))

#prepare data for deseq
data4pca_cell <- batch_cell[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table
data4pca_cm <- batch_cm[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table
data4pca_cell[is.na(data4pca_cell)] <- 0.1
data4pca_cm[is.na(data4pca_cm)] <- 0.1

data4pca_cell <- round(data4pca_cell)
data4pca_cm <- round(data4pca_cm)

rownames(data4pca_cell) <- batch_cell$Gene
rownames(data4pca_cm) <- batch_cm$Gene

#run deseq2
dds_cell <- DESeqDataSetFromMatrix(countData = data4pca_cell,
                                   colData = colData_cell,
                                   design = ~ treatment  )

dds_cm <- DESeqDataSetFromMatrix(countData = data4pca_cm,
                                 colData = colData_cm,
                                 design = ~ treatment  )

ddsr_cell <- DESeq(dds_cell)
ddsr_cm <- DESeq(dds_cm)


#filter for rows with minimum 10 reads
keep1 <- rowMeans(counts(ddsr_cell, normalized = TRUE)) >= 10
keep2 <- rowMeans(counts(ddsr_cm, normalized = TRUE)) >= 10
ddsr_cell1 <- ddsr_cell[keep1,]
ddsr_cm1 <- ddsr_cm[keep2,]

#results
res_cell <- results(ddsr_cell1, contrast=c("treatment", "PNK", "noPNK")) #noPNK as control = WT
res_cm <- results(ddsr_cm1, contrast=c("treatment", "PNK", "noPNK"))

table_cell <- as.data.frame(res_cell)
table_cell <- cbind(Geneid=rownames(table_cell), table_cell)

table_cm <- as.data.frame(res_cm)
table_cm <- cbind(Geneid=rownames(table_cm), table_cm)

# Select the top differentially expressed genes
top_genes_cell <- subset(res_cell, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.8)
top_genes_cm <- subset(res_cm, -log10(padj) > 1.3 & abs(log2FoldChange) > 0.8 )

# Get the count data for the top differentially expressed genes and the full dataset
counts_cell <- data.frame(cbind(Geneid = row.names(counts(ddsr_cell1)),counts(ddsr_cell1, normalized = TRUE)))
counts_cell1 <- as.data.frame(counts_cell[counts_cell$Geneid %in% rownames(top_genes_cell), ])
final_counts_cell <- inner_join(table_cell, counts_cell1, by="Geneid" )
final_tot_counts_cell <- inner_join(table_cell, counts_cell, by="Geneid" )

counts_cm <- data.frame(cbind(Geneid = row.names(counts(ddsr_cm1)),counts(ddsr_cm1, normalized = TRUE)))
counts_cm1 <- as.data.frame(counts_cm[counts_cm$Geneid %in% rownames(top_genes_cm), ])
final_counts_cm <- inner_join(table_cm, counts_cm1, by="Geneid" )
final_tot_counts_cm <- inner_join(table_cm, counts_cm, by="Geneid" )


write.table(final_tot_counts_cell, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/PNK_treatment/deseq_analysis/deseq2_Cell_PNK_treatment_norm_counts.txt", quote = F, row.names = F, sep = "\t")

write.table(final_tot_counts_cm, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/PNK_treatment/deseq_analysis/deseq2_CM_PNK_treatment_norm_counts.txt", quote = F, row.names = F, sep = "\t")



#venn diagramm different genes
library(VennDiagram)
library(RColorBrewer)

##selected genes with average (3 reps) bigger or equal to 5 norm reads (taken from final_tot_counts_)
noPNK <- read.table("~/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/2-PNK_treatment/intersection_PNK_noPNK/CM_noPNK_genes.txt", quote="\"", comment.char="")
PNK <- read.table("~/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/2-PNK_treatment/intersection_PNK_noPNK/CM_PNK_genes.txt", quote="\"", comment.char="")


noPNK_list <- noPNK$V1
PNK_list <- PNK$V1


venn.diagram(list(noPNK_list, PNK_list),
             category.names = c("no_PNK" , "PNK"),
             fill = c("skyblue", "#26C3EB"),
             filename = '~/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/2-PNK_treatment/intersection_PNK_noPNK/venn_diagramm_cm.png',
             main = "CM fraction",
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


#look for biotypes in different genes
unique_noPNK <- setdiff(noPNK_list, PNK_list)
write.table(unique_noPNK, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/2-PNK_treatment/intersection_PNK_noPNK/CM_unique_noPNK.txt", quote = F, row.names = F, sep = "\t")

# Extract non-intersecting genes in PNK_list
unique_PNK <- setdiff(PNK_list, noPNK_list)
write.table(unique_PNK, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/2-PNK_treatment/intersection_PNK_noPNK/CM_unique_PNK.txt", quote = F, row.names = F, sep = "\t")

#use this data to plot lost and gained fragments in cell and cm after PNK (biotype total prism file)



#### Figure S6E-F: fragment size comparison between RBP CLIP and RBP IP from CM, with motif enrichment
# to run in R
data <- read.table("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/read_length_fragmentomics/read_length_RBP_1_b.txt", quote="\"", comment.char="")

data2 <- read.table("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/read_length_fragmentomics/read_length_RBP_1_eCLIP.txt", quote="\"", comment.char="")

data$sample <- "RBP_1_IP"
data2$sample <- "RBP_1_eCLIP"

final <- rbind(data, data2)

#normalization for biotypes
final$V1[final$sample == "RBP_1_IP"] <- (final$V1[final$sample == "RBP_1_IP"]/sum(final$V1[final$sample == "RBP_1_IP"]))  *1000000
final$V1[final$sample == "RBP_1_eCLIP"] <- (final$V1[final$sample == "RBP_1_eCLIP"]/sum(final$V1[final$sample == "RBP_1_eCLIP"]))  *1000000

ggplot(final, aes(x = V2, y = V1, fill = sample)) +
  geom_area(alpha = 0.5, position = 'identity')+
  scale_fill_manual(values = c( "#86c774", "#4682B4")) +  # Assigning fill colors
  #scale_x_continuous(breaks = unique(as.numeric(factor(data$V2))), labels = levels(factor(data$V2))) +
  scale_x_continuous(n.breaks = 14) +
  labs(title = "RNA fragment sizes",
       x = "Length / nt",
       y = "RPM") +
  theme_classic(base_family = "Arial")+
  xlim(10,100)+
  theme(axis.text = element_text(size = 12),  # Adjust font size for axis text
        axis.title = element_text(size = 15),  # Adjust font size for axis titles
        plot.title = element_text(size = 24, face = "bold"))

#RBP_2_IPvsECLIP_a



#divide the two type of fragments bound by RBPs (<20 and >20nt)
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP

# Filter reads longer than 20 nt
awk '($3 - $2) > 20' filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_above20nt.bed
awk '($3 - $2) <= 20' filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_below20nt.bed

awk '($3 - $2) > 20' filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_above20nt.bed
awk '($3 - $2) <= 20' filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_below20nt.bed

awk '($3 - $2) > 20' WT_CM_input_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/WT_CM_input_smRNA_merged_sorted_MASTER_a_above20nt.bed
awk '($3 - $2) <= 20' WT_CM_input_smRNA_merged_sorted_MASTER_a.bed > ./read_length_fragmentomics/WT_CM_input_smRNA_merged_sorted_MASTER_a_below20nt.bed


#motif enrichment
seq2profile.pl GAATGD 0 RBP_1 > RBP_1.motif
seq2profile.pl AMAHWCA 0 RBP_2 > RBP_2.motif
seq2profile.pl VMAHWCA 0 RBP_2 > RBP_2_.motif

#RNAcompete: 359_12507992, RNCMPT00033: VMAHWCA, RNCMPT00172: AMAHWCA (RBP_2)
#other: CAUU
#RNAcompete: RNCMPT00076: GAAUGD, (RBP_1)
cd ./read_length_fragmentomics/

findMotifsGenome.pl filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_above20nt.bed hg38 outputDir_RBP_2_above20nt/ \
  -size given -p 4 -rna  \
  -bg WT_CM_input_smRNA_merged_sorted_MASTER_a_above20nt.bed  #-find RBP_2.motif -keepOverlappingBg

findMotifsGenome.pl filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_above20nt.bed hg38 outputDir_RBP_1_above20nt/ \
  -size given -p 4 -rna  \
  -bg WT_CM_input_smRNA_merged_sorted_MASTER_a_above20nt.bed  #-find RBP_2.motif -keepOverlappingBg

findMotifsGenome.pl filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_below20nt.bed hg38 outputDir_RBP_2_below20nt/ \
  -size given -p 4 -rna  \
  -bg WT_CM_input_smRNA_merged_sorted_MASTER_a_below20nt.bed  #-find RBP_2.motif -keepOverlappingBg

findMotifsGenome.pl filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_below20nt.bed hg38 outputDir_RBP_1_below20nt/ \
  -size given -p 4 -rna  \
  -bg WT_CM_input_smRNA_merged_sorted_MASTER_a_below20nt.bed  #-find RBP_2.motif -keepOverlappingBg



######## Figure S6G: intersection fragments bound to both RBPs (venn diagram)
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/

bedtools intersect -v -s -a filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed -b filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed | wc -l #unique_to_A.bed 7318
bedtools intersect -v -s -a filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed -b filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed | wc -l #unique_to_B.bed 2169
bedtools intersect -u -s -a filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed -b filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_merged.bed | wc -l #intersection 1070



######## Figure S6H: GO analysis of bound fragments
##GO analysis in R
awk '{print $8}' /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a_merged_annotation.txt | uniq > RBP_2_IP_MASTER_genes.txt

awk '{print $8}' /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a_merged_annotation.txt | uniq > RBP_1_IP_MASTER_genes.txt

#r
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


RBP_2_IP_genes <- read.table("~/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/RBP_2_IP_MASTER_genes.txt", quote="\"", comment.char="")
RBP_1_IP_genes <- read.table("~/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/RBP_1_IP_MASTER_genes.txt", quote="\"", comment.char="")

# Example: convert gene symbols to ENTREZ IDs
gene_ids_RBP_2 <- bitr(RBP_2_IP_genes$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_ids_RBP_1 <- bitr(RBP_1_IP_genes$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Example data
protein1_results <- as.data.frame(enrichGO(gene = gene_ids_RBP_2$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             readable = TRUE))
protein2_results <- as.data.frame(enrichGO(gene = gene_ids_RBP_1$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             readable = TRUE))

# Combine data and add protein labels
protein1_results$Protein <- "RBP_2"
protein2_results$Protein <- "RBP_1"
combined_results <- rbind(protein1_results, protein2_results)

top_n_per_type <- combined_results %>%
  group_by(Protein) %>%
  top_n(-10, p.adjust) %>%
  ungroup()

ggplot(top_n_per_type, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Protein, size = Count)) +
  geom_point(shape = 21, stroke = 1) +  # shape 21 is a filled circle with border; stroke controls border thickness
  labs(x = "GO Term", y = "-log10FDR", title = "Comparison of GO Enrichments") +
  theme_minimal(base_size = 18, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("RBP_2" = "darkblue", "RBP_1" = "cyan4")) +
  scale_size(range = c(3, 10))+
  guides(size = guide_legend(title = "Count")) +
  coord_flip()  # Flip coordinates for horizontal layout

#GO_BP_RBP_IP



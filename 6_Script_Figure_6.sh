
####### Figure 6B: RNA biotypes plots PNK vs noPNK treatment in WT cell and CM
#total reads for biotype found in excerpt_biotypes_raw.xlsx, plots made with prism ("biotype total")


####### Figure 6C: IGV tracks for IP enrichment compared to input and IgG control
##### better IGV profiles with RPM
cd /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/IGV_plots
#xattr -d com.apple.quarantine /Users/alessandra/bedGraphToBigWig

total_reads=$(samtools view -c IgG_IP_CM_smRNA_merged_sorted.bam)
bedtools genomecov -ibam IgG_IP_CM_smRNA_merged_sorted.bam -bg -scale $(echo "1000000 / $total_reads" | bc -l) > IgG_IP_CM_smRNA_merged_sorted_RPM.bedgraph
/Users/alessandra/bedGraphToBigWig IgG_IP_CM_smRNA_merged_sorted_RPM.bedgraph hg38.chrom.sizes.txt IgG_IP_CM_smRNA_merged_sorted_RPM.bw

total_reads=$(samtools view -c RBP_1_IP_CM_smRNA_merged_sorted.bam)
bedtools genomecov -ibam RBP_1_IP_CM_smRNA_merged_sorted.bam -bg -scale $(echo "1000000 / $total_reads" | bc -l) > RBP_1_IP_CM_smRNA_merged_sorted_RPM.bedgraph
/Users/alessandra/bedGraphToBigWig RBP_1_IP_CM_smRNA_merged_sorted_RPM.bedgraph hg38.chrom.sizes.txt RBP_1_IP_CM_smRNA_merged_sorted_RPM.bw

total_reads=$(samtools view -c RBP_2_IP_CM_smRNA_merged_sorted.bam)
bedtools genomecov -ibam RBP_2_IP_CM_smRNA_merged_sorted.bam -bg -scale $(echo "1000000 / $total_reads" | bc -l) > RBP_2_IP_CM_smRNA_merged_sorted_RPM.bedgraph
/Users/alessandra/bedGraphToBigWig RBP_2_IP_CM_smRNA_merged_sorted_RPM.bedgraph hg38.chrom.sizes.txt RBP_2_IP_CM_smRNA_merged_sorted_RPM.bw

total_reads=$(samtools view -c WT_CM_input_smRNA_merged_sorted.bam)
bedtools genomecov -ibam WT_CM_input_smRNA_merged_sorted.bam -bg -scale $(echo "1000000 / $total_reads" | bc -l) > WT_CM_input_smRNA_merged_sorted_RPM.bedgraph
/Users/alessandra/bedGraphToBigWig WT_CM_input_smRNA_merged_sorted_RPM.bedgraph hg38.chrom.sizes.txt WT_CM_input_smRNA_merged_sorted_RPM.bw




####### Figure 6D: fragment size distribution comparison between input and RBP IP
#to perform also for clip data for sup figures
cd /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP

#read lengths
awk '{print $3 - $2}' filtered_RBP_1_IP_CM_smRNA_merged_sorted_MASTER_a.bed | sort | uniq -c > read_length_RBP_1_a.txt

awk '{print $3 - $2}' /Volumes/AS_2023/scripts/eCLIP/bed_files/merged_cells/TARDBP.bed | sort | uniq -c > read_length_RBP_1_eCLIP.txt

awk '{print $3 - $2}' filtered_RBP_2_IP_CM_smRNA_merged_sorted_MASTER_a.bed | sort | uniq -c > read_length_RBP_2_a.txt

awk '{print $3 - $2}' /Volumes/AS_2023/scripts/eCLIP/bed_files/merged_cells/RBP_2.bed | sort | uniq -c > read_length_RBP_2_eCLIP.txt


#to run in R
data <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/read_length_fragmentomics/read_length_RBP_1_a.txt", quote="\"", comment.char="")

data2 <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/read_length_fragmentomics/read_length_RBP_2_a.txt", quote="\"", comment.char="")

data3 <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/read_length_fragmentomics/read_length_input_a.txt", quote="\"", comment.char="")

data$sample <- "RBP_1_IP"
data2$sample <- "RBP_2_IP"
data3$sample <- "Input_CM"

final <- rbind(data,data2, data3)

#normalization for sample
final$V1[final$sample == "RBP_1_IP"] <- (final$V1[final$sample == "RBP_1_IP"]/sum(final$V1[final$sample == "RBP_1_IP"]))  *1000000
final$V1[final$sample == "Input_CM"] <- (final$V1[final$sample == "Input_CM"]/sum(final$V1[final$sample == "Input_CM"]))  *1000000
final$V1[final$sample == "RBP_2_IP"] <- (final$V1[final$sample == "RBP_2_IP"]/sum(final$V1[final$sample == "RBP_2_IP"]))  *1000000

desired_order <- c("Input_CM", "RBP_1_IP", "RBP_2_IP")

# Convert 'sample' to a factor with the specified order
final$sample <- factor(final$sample, levels = desired_order)

ggplot(final, aes(x = V2, y = V1, color = sample, fill=sample)) +
  geom_area(alpha = 0.1, position = 'identity') +       # Transparent filled area
  geom_line(size = 1.2) +                               # Colored lines
  scale_fill_manual(values = c("#26C3EB", "#008D8DF4", "#00008D")) + # Assigning fill colors
  scale_color_manual(values = c("#26C3EB","#008D8DF4", "#00008D")) +# Assigning line colors
  scale_x_continuous(n.breaks = 14) +
  labs(title = "RNA Fragment Sizes",
       x = "Length / nt",
       y = "Density") +
  theme_classic(base_family = "Arial") +
  xlim(10, 100) +
  theme(
    axis.text = element_text(size = 12),       # Axis text size
    axis.title = element_text(size = 15),      # Axis title size
    plot.title = element_text(size = 24, face = "bold") # Title style
  )

#RBP_2_IPvsInput_a





###### Figure 6E: fragments annotation pie charts
#annotate remaining reads
cd /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/

Rscript merge_fragments_footprint.R WT_CM_input_smRNA_merged_sorted_MASTER_a.bed WT_CM_input_smRNA_merged_sorted_MASTER_a_merged1.bed

sortBed -i WT_CM_input_smRNA_merged_sorted_MASTER_a_merged1.bed > WT_CM_input_smRNA_merged_sorted_MASTER_a_merged.bed
rm WT_CM_input_smRNA_merged_sorted_MASTER_a_merged1.bed
 

bedtools intersect -s -f 1 -a WT_CM_input_smRNA_merged_sorted_MASTER_a_merged.bed \
                   -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/gencode.v45.annotation_utr_tRNAs_miRNAs.gtf -wa -wb \
| awk -F'\t' 'BEGIN {OFS="\t"} {
    gene_type = "";
    gene_name = "";

    # Extract transcript type and gene region from the annotation columns
    split($NF, fields, ";");
    for (i = 1; i <= length(fields); i++) {
        if (fields[i] ~ /gene_type/) {
            gsub(/.*gene_type "/, "", fields[i]);
            gsub(/".*/, "", fields[i]);
            gene_type = fields[i];
        }
        if (fields[i] ~ /gene_name/) {
            gsub(/.*gene_name "/, "", fields[i]);
            gsub(/".*/, "", fields[i]);
            gene_name = fields[i];
        }
    }

    # Output the original line with transcript type and gene region
    print $1, $2, $3, $4, $5, $6, $7, gene_name, gene_type, $10 ;
}' > WT_CM_input_smRNA_merged_sorted_MASTER_a_merged_annotation.txt


Rscript footprint_annotation.R /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/WT_CM_input_smRNA_merged_sorted_MASTER_a_merged_annotation.txt /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/WT_CM_input_smRNA_merged_sorted_MASTER_a_merged_annotation.txt > /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/WT_CM_input_smRNA_merged_sorted_MASTER_a_merged_annotation.log 2>&1

#do biotype pie chart in prism





######## Figure 6F: deseq analysis between RBP IP and input
#performed in R
library(dplyr)
library(DESeq2)
library(sva)

#perform deseq2 for different rna biotypes
ex_mapped <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/exceRpt_all_ReadCounts.txt")
colData_batch2 <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/colData_batch_IP_1.txt")
colData_batch2$genotype <- factor(colData_batch2$genotype)
colData_batch2$origin <- factor(colData_batch2$origin)
colData_batch2$treatment <- factor(colData_batch2$treatment)

#ex_mapped <- ex_mapped[-c(40963:40966),]
rownames(ex_mapped) <- ex_mapped$Gene

#select specific data for analysis
colData_IGF <- subset(colData_batch2, genotype %in% c("WT", "RBP_2") & origin %in% c("CM", "IP") & treatment == "PNK")
colData_TDP <- subset(colData_batch2, genotype %in% c("WT", "RBP_1") & origin %in% c("CM", "IP") & treatment == "PNK")
colData_IgG <- subset(colData_batch2, genotype %in% c("WT", "IgG") & origin %in% c("CM", "IP") & treatment == "PNK")

columns_to_keep_IGF <- c("Gene",colData_IGF$Sample)
columns_to_keep_TDP <- c("Gene",colData_TDP$Sample)
columns_to_keep_igg <- c("Gene",colData_IgG$Sample)

batch_IGF <- dplyr::select(ex_mapped, all_of(columns_to_keep_IGF))
batch_TDP <-  dplyr::select(ex_mapped, all_of(columns_to_keep_TDP))
batch_igg<-  dplyr::select(ex_mapped, all_of(columns_to_keep_igg))

#prepare data for deseq
data4pca_IGF <- batch_IGF[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table
data4pca_TDP <- batch_TDP[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table
data4pca_igg <- batch_igg[,-c(1)] + 0.1  #add 0.1 to all values in order to remove 0s from the table

data4pca_IGF[is.na(data4pca_IGF)] <- 0.1
data4pca_TDP[is.na(data4pca_TDP)] <- 0.1
data4pca_igg[is.na(data4pca_igg)] <- 0.1

data4pca_IGF <- round(data4pca_IGF)
data4pca_TDP <- round(data4pca_TDP)
data4pca_igg <- round(data4pca_igg)

rownames(data4pca_IGF) <- batch_IGF$Gene
rownames(data4pca_TDP) <- batch_TDP$Gene
rownames(data4pca_igg) <- batch_igg$Gene


#run deseq2
dds_igf <- DESeqDataSetFromMatrix(countData = data4pca_IGF,
                                   colData = colData_IGF,
                                   design = ~ genotype  )

dds_tdp <- DESeqDataSetFromMatrix(countData = data4pca_TDP, #[,-6]
                                 colData = colData_TDP, #[-6,]
                                 design = ~ genotype  )

dds_igg <- DESeqDataSetFromMatrix(countData = data4pca_igg,
                                 colData = colData_IgG,
                                 design = ~ genotype  )

ddsr_igf <- DESeq(dds_igf)
ddsr_tdp <- DESeq(dds_tdp)
ddsr_igg <- DESeq(dds_igg)


#filter for rows with minimum 10 reads
keep1 <- rowMeans(counts(ddsr_igf, normalized = TRUE)) >= 10
keep2 <- rowMeans(counts(ddsr_tdp, normalized = TRUE)) >= 10
keep3 <- rowMeans(counts(ddsr_igg, normalized = TRUE)) >= 10

ddsr_igf1 <- ddsr_igf[keep1,]
ddsr_tdp1 <- ddsr_tdp[keep2,]
ddsr_igg1 <- ddsr_igg[keep3,]

#results
res_igf <- results(ddsr_igf1, contrast=c("genotype", "RBP_2", "WT"))
res_tdp <- results(ddsr_tdp1, contrast=c("genotype", "RBP_1", "WT"))
res_igg <- results(ddsr_igg1, contrast=c("genotype", "IgG", "WT"))

table_igf <- as.data.frame(res_igf)
table_igf <- cbind(Geneid=rownames(table_igf), table_igf)

table_tdp <- as.data.frame(res_tdp)
table_tdp <- cbind(Geneid=rownames(table_tdp), table_tdp)

table_igg <- as.data.frame(res_igg)
table_igg <- cbind(Geneid=rownames(table_igg), table_igg)


# Select the top differentially expressed genes
top_genes_igf <- subset(res_igf, -log10(padj) > 1.3 & abs(log2FoldChange) > 1)
top_genes_tdp <- subset(res_tdp, -log10(padj) > 1.3 & abs(log2FoldChange) > 1 )
top_genes_igg <- subset(res_igg, -log10(padj) > 1.3 & abs(log2FoldChange) > 1 )


# Get the count data for the top differentially expressed genes and the full dataset
counts_igf <- data.frame(cbind(Geneid = row.names(counts(ddsr_igf1)),counts(ddsr_igf1, normalized = TRUE)))
counts_igf1 <- as.data.frame(counts_igf[counts_igf$Geneid %in% rownames(top_genes_igf), ])
final_counts_igf <- inner_join(table_igf, counts_igf1, by="Geneid" )
final_tot_counts_igf <- inner_join(table_igf, counts_igf, by="Geneid" )


counts_tdp <- data.frame(cbind(Geneid = row.names(counts(ddsr_tdp1)),counts(ddsr_tdp1, normalized = TRUE)))
counts_tdp1 <- as.data.frame(counts_tdp[counts_tdp$Geneid %in% rownames(top_genes_tdp), ])
final_counts_tdp <- inner_join(table_tdp, counts_tdp1, by="Geneid" )
final_tot_counts_tdp <- inner_join(table_tdp, counts_tdp, by="Geneid" )

counts_igg <- data.frame(cbind(Geneid = row.names(counts(ddsr_igg1)),counts(ddsr_igg1, normalized = TRUE)))
counts_igg1 <- as.data.frame(counts_igg[counts_igg$Geneid %in% rownames(top_genes_igg), ])
final_counts_igg <- inner_join(table_igg, counts_igg1, by="Geneid" )
final_tot_counts_igg <- inner_join(table_igg, counts_igg, by="Geneid" )


write.table(final_tot_counts_igf, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/deseq2_RBP_2_IP_vsInput_norm_counts.txt", quote = F, row.names = F, sep = "\t")

write.table(final_tot_counts_tdp, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/deseq2_RBP_1_IP_vsInput_norm_counts.txt", quote = F, row.names = F, sep = "\t")

write.table(final_tot_counts_igg, "/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/deseq2_IgG_IP_vsInput_norm_counts.txt", quote = F, row.names = F, sep = "\t")



##### scatterplots all data
ex_mapped[is.na(ex_mapped)] <- 0
rownames(ex_mapped) <- ex_mapped$Geneid
table1_filtered <- ex_mapped[!rowSums(ex_mapped[,-1]) == 0, ] #remove rna not expressed
total_reads <- colSums(table1_filtered[,-1])
rpm <- sweep(table1_filtered[,-1], 2, total_reads, "/") * 1e6 #calculate RPM for comparison
rpm <- rpm + 1

#subset the data for each RBP
data_igf <- rpm %>%  mutate(across(starts_with("L0724"), as.numeric)) %>%
  mutate(input = log(rowMeans(select(., L0724A10, L0724A11, L0724A12), na.rm = TRUE)),
         IP = log(rowMeans(select(., L0724A18, L0724A19, L0724A20), na.rm = TRUE)), group = "RBP_2")

data_tdp <- rpm %>%  mutate(across(starts_with("L0724"), as.numeric))  %>%
  mutate(input = log(rowMeans(select(., L0724A10, L0724A11, L0724A12), na.rm = TRUE)),
         IP = log(rowMeans(select(., L0724A21, L0724A22, L0724A23), na.rm = TRUE)), group = "RBP_1")

scatter <- rbind(data_igf[,22:24], data_tdp[,22:24])

# Create scatter plot
library(ggExtra)
library(ggplot2)

#single plot
ggplot(scatter[scatter$group == "RBP_2",]) +
  geom_point(aes(x = input, y = IP, color = IP - input > 2)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_classic(base_size = 15, base_family = "Helvetica") +
  labs(x = "Input \n log2(RPM+1)", y = "RBP_2 IP \n log2(RPM+1)") +
  coord_fixed(ratio = 1, xlim = c(0,10), ylim = c(0,10)) +
  theme(aspect.ratio = 1, legend.position="none")

#scatterplot_RBP_2_IP

ggplot(scatter[scatter$group == "RBP_1",]) +
  geom_point(aes(x = input, y = IP, color = IP - input > 2)) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_classic(base_size = 15, base_family = "Helvetica") +
  labs(x = "Input \n log2(RPM+1)", y = "RBP_1 IP \n log2(RPM+1)") +
  coord_fixed(ratio = 1, xlim = c(0,10), ylim = c(0,10)) +
  theme(aspect.ratio = 1, legend.position="none")

#scatterplot_RBP_1_IP




######## Figure 6G: newly detected RNA after RBP IP (not found in input and IgG)
#to do in R
wt_cm <- read.delim("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/PNK_treatment/deseq_analysis/deseq2_cm_noPNKvsPNK_treatment_norm_counts.txt", header=T)

## remove genes in IP found also in IgG
RBP_2_IP_genes <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/RBP_2_IP_DEGs.txt", quote="\"", comment.char="")
RBP_1_IP_genes <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/RBP_1_IP_DEGs.txt", quote="\"", comment.char="")
IgG_IP_genes <- read.table("/Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig6/RBP_IP/deseq_analysis/IgG_IP_DEGs.txt", quote="\"", comment.char="")

#remove any DEG from IP that is also present in IgG
filtered_RBP_1 <- anti_join(RBP_1_IP_genes, IgG_IP_genes, by = "V1")
filtered_RBP_2 <- anti_join(RBP_2_IP_genes, IgG_IP_genes, by = "V1")

#format RNA gene names and biotypes
split_values <- strsplit(filtered_RBP_1$V1, ":")
filtered_RBP_1$biotype <- sapply(split_values, function(x) x[2])
filtered_RBP_1$gene <- sapply(split_values, function(x) x[1])

split_values <- strsplit(filtered_RBP_2$V1, ":")
filtered_RBP_2$biotype <- sapply(split_values, function(x) x[2])
filtered_RBP_2$gene <- sapply(split_values, function(x) x[1])

#newly discovered fragments after IP (not present in input WT CM)
newdiscovery_RBP_1 <- filtered_RBP_1[!filtered_RBP_1$V1 %in% wt_tdp$Geneid, , drop = FALSE] #733
newdiscovery_RBP_2 <-filtered_RBP_2[!filtered_RBP_2$V1 %in% wt_tdp$Geneid, , drop = FALSE] #852

#format RNA gene names and biotypes
split_values <- strsplit(newdiscovery_RBP_1$V1, ":")
newdiscovery_RBP_1$biotype <- sapply(split_values, function(x) x[2])
newdiscovery_RBP_1$gene <- sapply(split_values, function(x) x[1])

split_values <- strsplit(newdiscovery_RBP_2$V1, ":")
newdiscovery_RBP_2$biotype <- sapply(split_values, function(x) x[2])
newdiscovery_RBP_2$gene <- sapply(split_values, function(x) x[1])

#extract bioype info to transfer in prism for plot
biotype_counts1 <- newdiscovery_RBP_1 %>% group_by(biotype) %>% summarize(count = n())
biotype_counts2 <- newdiscovery_RBP_2 %>% group_by(biotype) %>% summarize(count = n())

#transfer results to prism file "biotype total"





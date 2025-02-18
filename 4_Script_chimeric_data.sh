# create a unique bed file for chimeric data
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/GSE198251_RAW/

#intersect chimeric data from 4 experiments with our WT CM to find overlaps
bedtools intersect -a WT_CM_smRNA_AGObatch_merged_filtered_reliable.bed -b ./GSE198251_RAW/GSM5942090_Expt1_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942089_Expt1_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942093_Expt2_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942094_Expt2_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942095_Expt3_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942096_Expt3_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed ./GSE198251_RAW/GSM5942096_Expt3_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942097_Expt4_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz -s > WT_CM_chimericCLIP_all.bed #114436 -f 0.5



#annotation genes, biotypes –> use an annotation file cleaned from mRNA that overlap with other biotypes, in order to keep only high confidence mRNA
bedtools intersect -s -a WT_CM_chimericCLIP_all.bed \
                   -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/hg38_mrna_exons_no_miRNA_no_snRNA.bed -wa -wb \
| awk -F'\t' 'BEGIN {OFS="\t"} {
    gene_type = "";
    gene_name = "";

    # Extract transcript type and gene region from the annotation columns
    split($NF, fields, ";");
    for (i = 1; i <= length(fields); i++) {
        if (match(fields[i], /gene_type "([^"]+)"/, type)) {
                gene_type = type[1];
        }
        if (match(fields[i], /gene_name "([^"]+)"/, region)) {
                gene_name = region[1];
        }
    }

    # Output the original line with transcript type and gene region
    print $1, $2, $3, $4, $5, $6, $7, gene_name, gene_type, $15 ;
}' > WT_CM_chimericCLIP_all_highConf_mRNA_annotation.txt

library(dplyr)

#import annotated file as df
df <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/WT_CM_chimericCLIP_all_highConf_mRNA_annotation.txt", header=FALSE)

#remove redundant entries
df$V8[grepl("tRNA", df$V9)] <- df$V9[grepl("tRNA", df$V9)] 
df$V9[grepl("tRNA", df$V9)] <- "tRNA"
df$V10[grepl("tRNA", df$V9)] <- "ncRNA"
df$V10[grepl("snRNA", df$V9)] <- "ncRNA"
df$V10[grepl("snoRNA", df$V9)] <- "ncRNA"
df$V10[grepl("miRNA", df$V9)] <- "ncRNA"
df$V10[grepl("scaRNA", df$V9)] <- "ncRNA"
df$V10[grepl("piRNA", df$V9)] <- "ncRNA"
df$V10[grepl("misc_RNA", df$V9)] <- "ncRNA"
#df$V10[grepl("lncRNA", df$V9)] <- "ncRNA"
df$V10[grepl("Mt_rRNA", df$V9)] <- "ncRNA"
df$V10[grepl("vault_RNA", df$V9)] <- "ncRNA"

#remove redundant entries
df <- mutate(df, V8 = gsub("hsa-mir", "hsa-miR", V8))

unique_df <- unique(df)
unique_df2 <- unique_df[!(unique_df$V9 == "protein_coding_CDS_not_defined" | unique_df$V9 == "nonsense_mediated_decay" | unique_df$V9 == "retained_intron" | unique_df$V9 == "non_stop_decay" | unique_df$V9 == "" | unique_df$V9 == "artifact"), , drop = FALSE]

#get the best annotation for each peak
#selected using the following priority: CDS exon > 3′ UTR > 5′ UTR > protein-coding gene intron > noncod- ing RNA exon > noncoding RNA intron > intergenic.
   result_df <- unique_df2 %>%
    arrange(V1, V2, V3) %>%
    group_by(V1, V2, V3, V4, V5, V6, V7) %>%
    mutate(
      rank_V9 = match(V9, c("protein_coding","miRNA", "lncRNA","tRNA","scaRNA", "snRNA", "snoRNA", "misc_RNA", "piRNA","ribozyme", "Mt_rRNA")),
      rank_V10 = match(V10, c("ncRNA", "three_prime_utr", "five_prime_utr", "start_codon", "stop_codon", "exon", "transcript", "gene"))
    ) %>%
    filter_all(all_vars(!is.na(.))) %>%
    slice_min(rank_V9) %>%
    slice_min(rank_V10) %>%
    select(-rank_V9, -rank_V10) %>%
    ungroup() %>%
    mutate(length=V3-V2) %>%
    filter(length >= 15)

#check fragment size distribution for threshold 15 reads
	  ggplot(result_df, aes(x=log(length)+1)) + 
	    geom_density(fill="blue", alpha=0.5) +
	    ggtitle("Density Plot of Fragment Sizes") +
	    xlab("Fragment Size") +
	    ylab("Density") +
	    theme_classic() 
		#WT_CM_chimeric_15reads_length_distribution
		
write.table(result_df, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/WT_CM_chimericCLIP_all_15reads_annotation.bed", quote = F, col.names = F, row.names = F, sep = "\t")



############################# AGO data
#intersect with AGO KO CM to find overlaps
bedtools intersect -a AGO1234_KO_CM_smRNA_merged_filtered_sorted.bed -b ./GSE198251_RAW/GSM5942090_Expt1_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942089_Expt1_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942093_Expt2_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942094_Expt2_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942095_Expt3_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942096_Expt3_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed ./GSE198251_RAW/GSM5942096_Expt3_293T_NoGelChimeCLIP_total_rep2.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz ./GSE198251_RAW/GSM5942097_Expt4_293T_NoGelChimeCLIP_total_rep1.sorted.eclip.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.sorted.exclude-regions.bed.gz -s > AGO1234_CM_chimericCLIP_all.bed #114436 -f 0.5

WT_CM_chimericCLIP_all #114436
AGO1234_CM_chimericCLIP_all.bed #112899

#annotation genes, biotypes
bedtools intersect -s -a AGO1234_CM_chimericCLIP_all.bed \
                   -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/hg38_mrna_exons_no_miRNA_no_snRNA.bed -wa -wb \
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
    print $1, $2, $3, $4, $5, $6, $7, gene_name, gene_type, $15  ;
}' > AGO1234_CM_chimericCLIP_all_annotation.txt

library(dplyr)

#import annotated file as df
df <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/AGO1234_CM_chimericCLIP_all_annotation.txt", header=FALSE)

#remove redundant entries
df$V8[grepl("tRNA", df$V9)] <- df$V9[grepl("tRNA", df$V9)]
df$V9[grepl("tRNA", df$V9)] <- "tRNA"
df$V10[grepl("tRNA", df$V9)] <- "ncRNA"
df$V10[grepl("snRNA", df$V9)] <- "ncRNA"
df$V10[grepl("snoRNA", df$V9)] <- "ncRNA"
df$V10[grepl("miRNA", df$V9)] <- "ncRNA"
df$V10[grepl("scaRNA", df$V9)] <- "ncRNA"
df$V10[grepl("piRNA", df$V9)] <- "ncRNA"
df$V10[grepl("misc_RNA", df$V9)] <- "ncRNA"
#df$V10[grepl("lncRNA", df$V9)] <- "ncRNA"
df$V10[grepl("Mt_rRNA", df$V9)] <- "ncRNA"
df$V10[grepl("vault_RNA", df$V9)] <- "ncRNA"

#remove redundant entries
df <- mutate(df, V8 = gsub("hsa-mir", "hsa-miR", V8))

unique_df <- unique(df)
unique_df2 <- unique_df[!(unique_df$V9 == "protein_coding_CDS_not_defined" | unique_df$V9 == "nonsense_mediated_decay" | unique_df$V9 == "retained_intron" | unique_df$V9 == "non_stop_decay" | unique_df$V9 == "" | unique_df$V9 == "artifact"), , drop = FALSE]

#get the best annotation for each peak
#selected using the following priority: CDS exon > 3′ UTR > 5′ UTR > protein-coding gene intron > noncod- ing RNA exon > noncoding RNA intron > intergenic.
   result_df <- unique_df2 %>%
    arrange(V1, V2, V3) %>%
    group_by(V1, V2, V3, V4, V5, V6, V7) %>%
    mutate(
      rank_V9 = match(V9, c("protein_coding","miRNA", "lncRNA","tRNA","scaRNA", "snRNA", "snoRNA", "misc_RNA", "piRNA","ribozyme", "Mt_rRNA")),
      rank_V10 = match(V10, c("ncRNA", "three_prime_utr", "five_prime_utr", "start_codon", "stop_codon", "exon", "transcript", "gene"))
    ) %>%
    filter_all(all_vars(!is.na(.))) %>%
    slice_min(rank_V9) %>%
    slice_min(rank_V10) %>%
    select(-rank_V9, -rank_V10) %>%
    ungroup() %>%
    mutate(length=V3-V2) %>%
    filter(length >= 15)


ggplot(result_df, aes(x=log(length)+1)) +
        geom_density(fill="blue", alpha=0.5) +
        ggtitle("Density Plot of Fragment Sizes") +
        xlab("Fragment Size") +
        ylab("Density") +
        theme_classic()
        #AGO1234_CM_chimeric_15reads_length_distribution
        
write.table(result_df, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/AGO1234_CM_chimericCLIP_all_15reads_annotation.bed", quote = F, col.names = F, row.names = F, sep = "\t")



#######look for top20 most abundant miRNA in WT CM (also ago dependent), then generate with chimeric data a bed file with top20 mirna coordinates and associated mRNA

#select AGO dependent miRNAs
library(dplyr)
library(tidyr)

# data done with annotations and not coordinates
wt <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/WT_CM_chimericCLIP_all_15reads_annotation.bed", header=FALSE)
ago <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/AGO1234_CM_chimericCLIP_all_15reads_annotation.bed", header=FALSE)

#look for differential miRNAs between wt and ago <-  totally lost (63)
diff_miRNA <- anti_join(wt, ago, by = "V8") %>%
  filter(V9 == "miRNA") %>%
  mutate(Status = "total lost")
#diff_mRNA <- anti_join(wt, ago, by = "V8") %>% filter(V9 == "protein_coding")

#look for differential miRNAs between wt and ago <-  >80% lost coverage
ago_miRNA <- ago %>% filter(V9 == "miRNA") #%>% select()

# Join the WT and AGO datasets on miRNA names
combined_data <- wt %>%
  filter(V9 == "miRNA") %>%
  inner_join(ago_miRNA, by = "V8", suffix = c("_wt", "_ago"))

# Calculate the percentage decrease in coverage and filter miRNAs with more than 80% loss (144 miRNAs)
miRNAs_with_loss <- combined_data %>%
  mutate(percentage_loss = (V7_wt - V7_ago) / V7_wt * 100) %>%
  filter(percentage_loss > 80) %>%
  select(V8, V7_wt, V7_ago, percentage_loss ) %>%
  distinct() %>%
  mutate(Status = "more 80% lost")

miRNA_CM_lost <- rbind(miRNAs_with_loss[,c("V8","Status")], diff_miRNA[,c("V8","Status")])
write.table(miRNA_CM_lost, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/lost_miRNAs_AGO.txt", quote = F, col.names = F, row.names = F, sep = "\t")


## create chimera file for selected miRNA-mRNA couples

final_mirnaTarget <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/top10_miRNA_chimera_EXTRA.bed", header=FALSE)
final_mirnaTarget <-  final_mirnaTarget %>% rename("V4" = "miRNA", "V5" = "V8")
mirnas <- read.delim("~/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/GSE198251_RAW/hsa-miR-26a-5p_annotation.bed", header=FALSE)
mirnas1 <- mirnas %>% filter(V9 == "protein_coding")  %>% mutate(miRNA = "hsa-miR-26a-5p") %>% select(V1, V2, V3, miRNA, V8, V6) %>% distinct()
final_mirnaTarget <- rbind(final_mirnaTarget, mirnas1)

write.table(final_mirnaTarget, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/top10_miRNA_chimera_EXTRA.bed", quote = F, col.names = F, row.names = F, sep = "\t")



#double check new annotation from bogdan – keep only high confidence mRNA
bedtools intersect -s -a top10_miRNA_chimera_EXTRA.bed \
                   -b hg38_mrna_exons_no_miRNA_no_snRNA.bed -wa -wb | awk -F'\t' 'BEGIN {OFS="\t"} {
    gene_name = "";

    # Extract transcript type and gene region from the annotation columns
    split($NF, fields, ";");
    for (i = 1; i <= length(fields); i++) {
        if (fields[i] ~ /gene_name/) {
            gsub(/.*gene_name "/, "", fields[i]);
            gsub(/".*/, "", fields[i]);
            gene_name = fields[i];
        }
    }

    # Output the original line with transcript type and gene region
    print $1, $2, $3, $4, $5, $6, gene_name ;
}' | uniq > top10_miRNA_chimera_mRNA_EXTRA.bed

### annotations are good!


#in R
#get target and miRNA coverages, use 3 replicates to make statistics for correlation
top10_mirna_cm <- read.delim("~/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/WT_CM_smRNA_AGObatch_top10_chimera_miRNA_EXTRA.bed", header=FALSE)
colnames(top10_mirna_cm) <- c("chr", "start", "end", "WT_batch", "coverage", "strand", "miRNA", "target")

top10_mirna_cell <- read.delim("~/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/WT_Cell_smRNA_AGObatch_top10_chimera_miRNA_EXTRA.bed", header=FALSE)
colnames(top10_mirna_cell) <- c("chr", "start", "end", "WT_batch", "coverage", "strand", "miRNA", "target")

#replicates
wt_cm1 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/WT_CM_smRNA_AGObatch_merged_filtered_reliable_annotation.bed", header=FALSE) #L0423A9_WT_CM_smRNA_AGObatch_filtered_sorted_annotation.bed, WT_CM_smRNA_AGObatch_merged_filtered_reliable_annotation.bed
wt_cm2 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/WT_CM_smRNA_RBPbatch_merged_filtered_reliable_annotation.bed", header=FALSE) #L0423A41_WT_CM, WT_CM_smRNA_RBPbatch_merged_filtered_reliable_annotation.bed
wt_cm3 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/WT_CM_smRNA_GAbatch_merged_filtered_reliable_annotation.bed", header=FALSE) #L0922B3_WT_CM, WT_CM_smRNA_GAbatch_merged_filtered_reliable_annotation.bed
colnames(wt_cm1) <- c("chr", "start", "end", "WT_batch", ".", "strand","coverage", "miRNA", "biotype", "feature")
colnames(wt_cm2) <- c("chr", "start", "end", "WT_batch", ".", "strand","coverage",  "miRNA", "biotype", "feature")
colnames(wt_cm3) <- c("chr", "start", "end", "WT_batch", ".", "strand","coverage",  "miRNA", "biotype", "feature")
wt_cm1 <- wt_cm1 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)
wt_cm2 <- wt_cm2 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)
wt_cm3 <- wt_cm3 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)

ago_cm1 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/AGO1234_CM_L0423A16_filtered_sorted_annotation.bed", header=FALSE)
ago_cm2 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/AGO1234_CM_L0423A32_filtered_sorted_annotation.bed", header=FALSE)
ago_cm3 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/data/AGO1234_CM_L0423A48_filtered_sorted_annotation.bed", header=FALSE)
colnames(ago_cm1) <- c("chr", "start", "end", "batch", ".", "strand","coverage", "miRNA", "biotype", "feature")
colnames(ago_cm2) <- c("chr", "start", "end", "batch", ".", "strand","coverage",  "miRNA", "biotype", "feature")
colnames(ago_cm3) <- c("chr", "start", "end", "batch", ".", "strand","coverage",  "miRNA", "biotype", "feature")
ago_cm1 <- ago_cm1 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)
ago_cm2 <- ago_cm2 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)
ago_cm3 <- ago_cm3 %>% mutate(coverage = coverage/(sum(coverage)* 1e6)) #%>%  select(-length)


#rename each target fragment (analysis based on fragments and not full transcipts)
top10_mirna_cm1 <- top10_mirna_cm %>%
  mutate(coord_group = paste(target, chr, start, end, sep="_")) %>%
  group_by(target) %>%  # Group only by target to ensure numbering is unique within each target
  #mutate(fragment_number = dense_rank(coord_group)) %>%
  mutate(target = paste(target, dense_rank(coord_group), sep="_")) %>%
  ungroup()  %>%
  select(-WT_batch,-coord_group,-coverage) %>% distinct()


############## CM: assign coverage mirna and target for each replicate
gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(wt_cm1[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal", ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm1 <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) #%>% group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm1 <- wt_cm1 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-WT_batch,-.,-biotype, -feature)
miRNA_cm1 <- miRNA_cm1[-c(18,13,3,27,9,14,26),] #6,12,13,14
rep1 <- mRNA_cm1 %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm1[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep1")

gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(wt_cm2[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal", ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm2 <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) #%>% group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm2 <- wt_cm2 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-WT_batch,-.,-biotype, -feature)
miRNA_cm2 <- miRNA_cm2[-c(18,13,3,27,9,14,26),] #6,12,13
rep2 <- mRNA_cm2 %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm2[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep2")

gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(wt_cm3[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal",  ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm3 <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) #%>% group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm3 <- wt_cm3 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-WT_batch,-.,-biotype, -feature)
miRNA_cm3 <- miRNA_cm3[-c(18,13,3,27,9,14,26),]
rep3 <- mRNA_cm3 %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm3[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep3")

#get coverage replicates for each mirna-target CM AGO
gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(ago_cm1[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal", ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm1_ago <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) %>%  group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm1_ago <- ago_cm1 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-batch,-.,-biotype, -feature)
miRNA_cm1_ago <- miRNA_cm1_ago[-c(18,13,3,27,9,14,26),]
rep1_ago <- mRNA_cm1_ago %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm1_ago[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep1")

gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(ago_cm2[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal", ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm2_ago <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) %>%  group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm2_ago <- ago_cm2 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-batch,-.,-biotype, -feature)
miRNA_cm2_ago <- miRNA_cm2_ago[-c(18,13,3,27,9,14,26),]
rep2_ago <- mRNA_cm2_ago %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm2_ago[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep2")

gr1 <- makeGRangesFromDataFrame(top10_mirna_cm1, keep.extra.columns = TRUE)
gr2 <- makeGRangesFromDataFrame(ago_cm3[, -c(4,5, 8:11)], keep.extra.columns = TRUE)
overlaps <- findOverlaps(gr1, gr2, maxgap = 20, minoverlap = 10, type= "equal", ignore.strand = FALSE)
df_gr1 <- as.data.frame(gr1[queryHits(overlaps)])
df_gr2 <- as.data.frame(gr2[subjectHits(overlaps)])
mRNA_cm3_ago <- cbind(df_gr1, df_gr2)[,c(8:10,12,6,7,13, 11)] %>% filter(width < 51) %>% select(-width) %>%  group_by(seqnames, start, end, strand, miRNA) %>% mutate(coverage = sum(coverage)) %>%  arrange(target) %>% filter(row_number() == 1) %>%  ungroup()  # Keep only the top entry per group

miRNA_cm3_ago <- ago_cm3 %>% filter(miRNA %in% top10_mirna_cm1$miRNA) %>% select(-batch,-.,-biotype, -feature)
miRNA_cm3_ago <- miRNA_cm3_ago[-c(18,13,3,27,9,14,26),]
rep3_ago <- mRNA_cm3_ago %>% select(target, miRNA, coverage) %>% distinct() %>% left_join(miRNA_cm3_ago[,c("miRNA", "coverage")], by = "miRNA", suffix = c("_target", "_miRNA")) %>% drop_na() %>% mutate(rep = "rep3")

#only get bindings present in all the replicates (or at least 2!!!!) –-> needed for statistics
total <- rep1 %>%
  full_join(rep2, by = c("target", "miRNA")) %>%
  full_join(rep3, by = c("target", "miRNA")) %>%
  mutate(non_na_count = rowSums(!is.na(select(., contains("coverage_target"))))) %>%
  filter(non_na_count >= 3) %>%
  mutate(miRNA_target = paste(miRNA, target, sep = "_")) %>%
  select(-non_na_count) %>% distinct() # Optional: remove the count column if it's no longer needed

total_ago <- rep1_ago %>%
  full_join(rep2_ago, by = c("target", "miRNA")) %>%
  full_join(rep3_ago, by = c("target", "miRNA")) %>%
  mutate(non_na_count = rowSums(!is.na(select(., contains("coverage_target"))))) %>%
  filter(non_na_count >= 3) %>%
  mutate(miRNA_target = paste(miRNA, target, sep = "_")) %>%
  select(-non_na_count) %>% distinct() # Optional: remove the count column if it's no longer needed

#scatterplot
total_scatterplot <- rbind(rep1,rep2,rep3)
filtered_scatterplot <- total_scatterplot %>% mutate(miRNA_target = paste(miRNA, target, sep = "_")) %>%  filter(miRNA_target %in% total$miRNA_target)
total_scatterplot_ago <- rbind(rep1_ago,rep2_ago,rep3_ago)
filtered_scatterplot_ago <- total_scatterplot_ago %>% mutate(miRNA_target = paste(miRNA, target, sep = "_")) %>%  filter(miRNA_target %in% total_ago$miRNA_target)


#create new df to collect correlations for each couple
results_df_cm <- tibble(miRNA_target = character(),
                        miRNA = factor(),
                        cor_value_cm = double())
                        #cor_p_value_cm = double())

results_df_cm_ago <- tibble(miRNA_target = character(),
                            miRNA = factor(),
                            cor_value_cm_ago = double())
                        #cor_p_value_cm_ago = double())

# Iterate over each miRNA_target and replicate and calculate correlation <- do it for each WT-KO comparison you need
for (Target in unique(filtered_scatterplot$miRNA_target)) {
  # Filter data for WT and AGO1234 for the current target and replicate
  filtered_data <- filtered_scatterplot %>%
    filter(miRNA_target == Target)
  
  if (nrow(filtered_data) == 0) next
  
  # Compute changes
  cor_value = cor(filtered_data$coverage_miRNA, filtered_data$coverage_target, use = "complete.obs", method = "pearson")
  #cor_p_value = cor.test(filtered_data$coverage_miRNA, filtered_data$coverage_target, method = "pearson")$p.value
  
  # Append to results dataframe
  results_df_cm <- results_df_cm %>%
    add_row(miRNA_target = Target,
            miRNA = unique(filtered_data$miRNA),
            cor_value_cm = if_else(length(cor_value) > 0, cor_value, NA_real_))
            #cor_p_value_cm = if_else(length(cor_p_value) > 0, cor_p_value, NA_real_))
}

for (Target in unique(filtered_scatterplot_ago$miRNA_target)) {
  # Filter data for WT and AGO1234 for the current target and replicate
  filtered_data <- filtered_scatterplot_ago %>%
    filter(miRNA_target == Target)
  
  if (nrow(filtered_data) == 0) next
  
  # Compute changes
  cor_value = cor(filtered_data$coverage_miRNA, filtered_data$coverage_target, use = "complete.obs", method = "pearson")
  #cor_p_value = cor.test(filtered_data$coverage_miRNA, filtered_data$coverage_target, method = "pearson")$p.value
  
  # Append to results dataframe
  results_df_cm_ago <- results_df_cm_ago %>%
    add_row(miRNA_target = Target,
            miRNA = unique(filtered_data$miRNA),
            cor_value_cm_ago = if_else(length(cor_value) > 0, cor_value, NA_real_))
            #cor_p_value_cm_ago = if_else(length(cor_p_value) > 0, cor_p_value, NA_real_))
}

results_df_cm$miRNA <- as.factor(results_df_cm$miRNA)
levels(results_df_cm$miRNA)
results_df_cm$miRNA <- droplevels(results_df_cm$miRNA)
mean_cor_values <- results_df_cm %>%
  group_by(miRNA) %>%
  summarize(mean_cor_value = mean(cor_value_cm, na.rm = TRUE))
results_df_cm <- results_df_cm %>%
  left_join(mean_cor_values, by = "miRNA")



write.table(results_df_cell, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/correlation-results_cell_WT.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(results_df_cm, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/correlation-results_cm_WT.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(results_df_cm_ago, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/correlation-results_cm_ago.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(results_df_cell_ago, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/correlation-results_cell_ago.txt", quote = F, col.names = T, row.names = F, sep = "\t")


#for permutations in fig 4E
input <- rbind(miRNA_cm1,miRNA_cm2, miRNA_cm3,miRNA_cm1,miRNA_cm2, miRNA_cm3, miRNA_cm1,miRNA_cm2, miRNA_cm3, miRNA_cm1,miRNA_cm2, miRNA_cm3,miRNA_cm1,miRNA_cm2, miRNA_cm3, miRNA_cm1,miRNA_cm2, miRNA_cm3) %>% mutate(coverage_miRNA = coverage) %>% select(-coverage)

write.table(input, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/input_permutations.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#for plot in fig 4F
result1 <- total %>%
  group_by(miRNA) %>%
  distinct()  %>%
  summarise(unique_targets = n_distinct(target)) %>% mutate(subcategory = "WT")

result2 <- total_ago %>%
  group_by(miRNA) %>%
  distinct()  %>%
  summarise(unique_targets = n_distinct(target)) %>% mutate(subcategory = "AGO mutant")

result3 <- rbind(result1, result2)

write.table(result3, "/Users/alessandra/polybox/sequencing/RNA-seq_analysis/miRNA_targets_analysis/chimeric_eCLIP/No_targets_top20_dependencyAGO.txt", quote = F, col.names = T, row.names = F, sep = "\t")



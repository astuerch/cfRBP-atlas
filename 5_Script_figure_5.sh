
######### FIGURE 5B: cfRBP footprint characteristics

CHECK IF FIGURE FILE HAS DONE STEP 167 OF MERGE FRAGMENTS 5NT BEFORE DOING FINAL FOOTPRINT FILE WITH ALL RBPS, BECAUSE NUMBERS DONT MATCH !!!!!


#count footprints
#count genes

#plots done in stats_footprints.prism


###fragment size panel – violin plot
awk 'BEGIN{OFS=FS="\t"} {print $3 - $2, $8}' footprint_allRBPs_master_fragments_sorted.bed > ./3_Footprint_characterization/fragment_size/all_RBPs_master_length.txt

#R
reads <- read.table("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/3_Footprint_characterization/fragment_size/all_RBPs_master_length.txt", quote="\"", comment.char="")
desired_order <- c("RBP_1", "RBP_2", "RBP_3", "RBP_4", "RBP_5", "RBP_6")

# Convert reads$V2 to a factor with the specified levels
reads$V2 <- factor(reads$V2, levels = desired_order)

ggplot(reads, aes(x = V2, y = V1)) +
  geom_violin( trim = T, fill = "grey") +
  #geom_jitter(shape=16, position=position_jitter(0.2))
  geom_boxplot(width=0.1) +
  #ylim(0,50)+
 # geom_jitter(width = 0.1, alpha = 0.6, color = "grey") +  # All data points
  #stat_summary(fun = "median", geom = "point", shape = 20, size = 3, fill = "red") +  # Mean point
  labs( x = "", y = "Fragment length [nt]") +
  theme_classic()







######### FIGURE 5D: cfRBP footprint intersection with plasma dataset
#no of unique peaks in CM
awk -F'\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,$6}'  WT_CM_smRNA_master_fragments.bed | uniq | wc -l #241296
#no of unique footprint peaks
awk -F'\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,$6}' footprint_allRBPs_master_fragments_sorted.bed | uniq | wc -l #81602


#prepare plasma samples (from Detector-seq)
#reliable peaks
cd /Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/
bed_file1="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/plasma1_sorted_filtered.bed"
bed_file2="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/plasma2_sorted_filtered.bed"
bed_file3="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/plasma3_sorted_filtered.bed"
bed_file4="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/plasma4_sorted_filtered.bed"
bed_file5="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/plasma5_sorted_filtered.bed"

input_bed="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/Plasma_detector_merged_filtered.bed"
output_bed="/Volumes/AS_2023/4_other_DB/Detector-seq_processed_files/Plasma_detector_merged_filtered_reliable.bed"

bedtools multiinter -i $input_bed $bed_file1 $bed_file2 $bed_file3 $bed_file4 $bed_file5 -names input bed_file1 bed_file2 bed_file3 bed_file4 bed_file5 -s | awk '($4 >= 6)' > reliable_peaks1.bed

bedtools intersect -u -sorted -a $input_bed -b reliable_peaks1.bed  > $output_bed


#####intersect biofluid with cm and rbp footprints
bedtools intersect -u -s -sorted -a Plasma_detector_merged_filtered_reliable.bed -b WT_CM_smRNA_master_fragments.bed > PlasmaDetector_WTcm_intersection.bed

bedtools intersect -u -s -sorted -a Plasma_detector_merged_filtered_reliable.bed -b footprint_allRBPs_master_fragments_sorted.bed > PlasmaDetector_RBPfootprint_intersection.bed


#calculate percentage of overlapping peaks for different filters
wc -l Plasma_detector_merged_filtered_reliable.bed #236780
wc -l WT_CM_smRNA_master_fragments.bed #241296

wc -l PlasmaDetector_WTcm_intersection.bed #17882
wc -l PlasmaDetector_RBPfootprint_intersection.bed #14377

#constructed venn diagram



######### FIGURE 5E: cfRBP footprint upsetPlot with plasma dataset
#this figure is based on gene annotations, not on fragments
#annotate files
bedtools intersect -s -f 1 -a WT_CM_smRNA_master_fragments.bed \
                   -b gencode.v45.annotation_utr_tRNAs_miRNAs.gtf -wa -wb \
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
    print $1, $2, $3, $8, $5, $6, $7, gene_name, gene_type, $10 ;
}' > WT_CM_smRNA_master_fragments_annotation.txt

Rscript footprint_annotation.R WT_CM_smRNA_master_fragments_annotation.txt WT_CM_smRNA_master_fragments_annotation.txt > WT_CM_smRNA_master_fragments_annotation.log 2>&1



bedtools intersect -s -f 1 -a ./2_Footprint_biofluids/Plasma_detector_merged_filtered_reliable.bed \
                   -b gencode.v45.annotation_utr_tRNAs_miRNAs.gtf -wa -wb \
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
}' > Plasma_detector_merged_filtered_reliable_annotation.txt


Rscript footprint_annotation.R Plasma_detector_merged_filtered_reliable_annotation.txt Plasma_detector_merged_filtered_reliable_annotation.txt > Plasma_detector_merged_filtered_reliable_annotation.log 2>&1


#put together unique genes for each source to use for upset plot – data_upsetplot.txt
data1 <- read.delim("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_master_fragments_annotation.txt", header=FALSE)
data2 <- read.delim("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/footprint_allRBPs_master_fragments_sorted_annotation.txt", header=FALSE) #to do for each
data3 <- read.delim("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/Plasma_detector_merged_filtered_reliable_annotation.txt", header=FALSE) #to do for each

genes1 <- data %>% select(V8) %>% unique %>% mutate(Gene = V8, source = "HEK_CM") %>% select(-V8)
genes2 <- data %>% select(V8) %>% unique %>% mutate(Gene = V8, source = "RBP_footprint") %>% select(-V8)
genes3 <- data %>% select(V8) %>% unique %>% mutate(Gene = V8, source = "Plasma") %>% select(-V8)

data1 <- rbind(genes1,genes2,genes3)
colnames(data1) <- c("Gene", "source")

write.table(data1, "~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/data_upsetplot.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#prepare upsetplot
library(UpSetR)
library(dplyr)

data <- read.delim("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/data_upsetplot.txt")

# Grouping genes by KO and aggregating them into lists
grouped <- data %>%
  group_by(source) %>%
  summarize(Gene = list(unique(Gene)))

# Create a list of gene sets for each KO
gene_sets <- grouped$Gene

# Create an empty list to store the genes for each KO
gene_list <- list()

# Iterate through each gene set and create a list of genes for each KO
for (i in 1:length(gene_sets)) {
  gene_list[[grouped$source[i]]] <- gene_sets[[i]]
}

# Generate UpSet plot
upset(fromList(gene_list), sets.bar.color = "grey60",
      nsets = 50, order.by = c("freq"), nintersects = 50, cutoff = 10, point.size = 6, line.size = 1.2,
      mainbar.y.label = "Gene Intersections", sets.x.label = "cfRNA source",
      text.scale = c(3, 2.5, 2.5, 2, 3, 2.5))






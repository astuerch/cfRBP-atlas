##### PROCESSING OF RBP IP SEQUENCING
# run sequencing data with excerpt pipeline adapted for Nextflex sRNA v4

cd /cluster/work/mansuy/mateescu_lab/alessandra/exceRpt_smallRNA
module load stack/2024-06 jdk/8u141-b15 bowtie2/2.5.1-u2j3omo fastx-toolkit/0.0.14 samtools/1.17 fastqc/0.12.1 sra-tools/3.0.3 r/4.3.2

export PATH=$PATH:/cluster/work/mansuy/mateescu_lab/alessandra/exceRpt_smallRNA/STAR-2.7.0a/bin/Linux_x86_64


mkdir /cluster/scratch/astuerch/output_RBP_IP

for d in /cluster/work/mansuy/mateescu_lab/alessandra/fastq_files/Libraries_IP_RBP_PNK_1/*
do
  sbatch --time=10:00:00 -n 20 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'make -f /cluster/work/mansuy/mateescu_lab/alessandra/exceRpt_smallRNA/exceRpt_smallRNA_NextFlex_sRNA_v4_12bp \
       INPUT_FILE_PATH=$d \
       OUTPUT_DIR=/cluster/scratch/astuerch/output_RBP_IP' "
done


#prepare diagnostic file
cd /cluster/work/mansuy/mateescu_lab/alessandra/exceRpt_smallRNA
module load stack/2024-06 r/4.4.0

mkdir /cluster/work/mansuy/mateescu_lab/alessandra/excerpt_diagnostic_RBP_IP_1/

R #start interactive session
source("mergePipelineRuns_functions_new.R") #added data.dir pathway directly in the script
library(plyr)
library(reshape2)
library(ggplot2)
processSamplesInDir(data.dir, output.dir) #change variables manually in the script



#QC libraries
#PCA analysis in R
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load the RNA-seq data
table1 <- read.delim("/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/exceRpt_all_ReadCounts.txt")
table1[is.na(table1)] <- 0
table1 <- round(table1[,-c(1)])
#table1 <- table1[,-46]


# Create a DESeq2 object
colData_batch1 <- read.delim("/Users/alessandra/polybox/sequencing/RBP_IP_1/colData_batch_IP_1.txt") ##ORDER OF SAMPLES HAS TO BE THE SAME AS IN THE TABLE DATA!!!!!
colData_batch1$source <- paste(colData_batch1$Sample, colData_batch1$genotype, sep="_")
colData_batch1 <- subset(colData_batch1, genotype == 'WT')
columns_to_keep_cm <- c(colData_batch1$Sample)

batch_cm <- dplyr::select(table1, all_of(columns_to_keep_cm))
#colData_batch1 <- colData_batch1[-46,]
keep2 <- rowMeans(batch_cm) >= 5
data <- batch_cm[keep2,]


dds <- DESeqDataSetFromMatrix(batch_cm, colData_batch1, design = ~ treatment)

# Normalize the data
dds <- DESeq(dds)

# Perform variance stabilization (PCA for genes)
vsd <- varianceStabilizingTransformation(dds)

# Transpose the variance-stabilized data matrix (PCA for samples)
vsd_t <- t(assay(vsd))

# Run PCA analysis
pca <- prcomp(vsd_t)

# Convert the PCA results to a data frame
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = colData_batch1$treatment, Origin = colData_batch1$source)
variance_explained <- round(summary(pca)$importance[2, c("PC1", "PC2")] * 100, 2)

# Create a PCA plot using ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(Sample), label = Origin)) +
  geom_point() +
  geom_text_repel(size = 3, nudge_y = 1, arrow = arrow(length = unit(0.1, "inches"), type = "closed"), min.segment.length = 0.5) +
  labs(x = "PC1", y = "PC2", title = paste("PCA Plot (PC1:", variance_explained[1], "%, PC2:", variance_explained[2], "%)")) +
  theme_bw()

#PCA_smRNA_RBP_IP
#5x7


########## process data to find high confidence protected RNA fragments
cd /cluster/work/mansuy/mateescu_lab/alessandra/BAM_files/Bam_files_IP_target

module load stack/2024-06 samtools/1.17 bedtools2/2.31.0

# extract bam files from excerpt output
for d in /cluster/scratch/astuerch/output_RBP_IP/*/
do
cd $d
cp *endogenousAlignments_genome_Aligned.out.bam /cluster/scratch/astuerch/output_RBP_IP/bam_files/
done


#combine replicates
sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools merge RBP_2_IP_CM_smRNA_merged.bam L0724A18_endogenousAlignments_genome_Aligned.out.bam L0724A19_endogenousAlignments_genome_Aligned.out.bam L0724A20_endogenousAlignments_genome_Aligned.out.bam'" #samtools merge with bam files only concatenates the input files without changing their coordinates

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools merge RBP_1_IP_CM_smRNA_merged.bam L0724A21_endogenousAlignments_genome_Aligned.out.bam L0724A22_endogenousAlignments_genome_Aligned.out.bam L0724A23_endogenousAlignments_genome_Aligned.out.bam'" #samtools merge with bam files only concatenates the input files without changing their coordinates

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools merge IgG_IP_CM_smRNA_merged.bam L0724A15_endogenousAlignments_genome_Aligned.out.bam L0724A16_endogenousAlignments_genome_Aligned.out.bam L0724A17_endogenousAlignments_genome_Aligned.out.bam'" #samtools merge with bam files only concatenates the input files without changing their coordinates

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools merge WT_CM_input_smRNA_merged.bam L0724A10_endogenousAlignments_genome_Aligned.out.bam L0724A11_endogenousAlignments_genome_Aligned.out.bam L0724A12_endogenousAlignments_genome_Aligned.out.bam'" #samtools merge with bam files only concatenates the input files without changing their coordinates


#sort merged files
sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools sort RBP_1_IP_CM_smRNA_merged.bam > RBP_1_IP_CM_smRNA_merged_sorted.bam'"

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools sort RBP_2_IP_CM_smRNA_merged.bam > RBP_2_IP_CM_smRNA_merged_sorted.bam'"

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools sort IgG_IP_CM_smRNA_merged.bam > IgG_IP_CM_smRNA_merged_sorted.bam'"

sbatch --time=1:00:00 -n 4 --mem-per-cpu=2048 --tmp=40000  --wrap="bash -c 'samtools sort WT_CM_input_smRNA_merged.bam > WT_CM_input_smRNA_merged_sorted.bam'"


## clean IP bam files from IgG reads
#get bed files
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/read_length_fragmentomics

bedtools bamtobed -i /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/IGV_plots/RBP_1_IP_CM_smRNA_merged_sorted.bam > RBP_1_IP_CM_smRNA_merged_sorted_readsID.bed
bedtools bamtobed -i /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/IGV_plots/RBP_2_IP_CM_smRNA_merged_sorted.bam > RBP_2_IP_CM_smRNA_merged_sorted_readsID.bed
bedtools bamtobed -i /Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/IGV_plots/IgG_IP_CM_smRNA_merged_sorted.bam > IgG_IP_CM_smRNA_merged_sorted_readsID.bed

##intersect to find common reads between IP and IgG
cd /cluster/scratch/astuerch/
module load stack/2024-06 samtools/1.17 bedtools2/2.31.0

#sort files
for file in *sorted_readsID.bed; do
        NAME=$(basename $file .bed)
        sbatch --time=4:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'sort -k1,1 -k2,2n $file -o ${NAME}_2.bed'"
done

#files too big, need to be split
split -l 10000000 RBP_2_IP_CM_smRNA_merged_sorted_readsID_2.bed RBP_2_chunk_
split -l 10000000 RBP_1_IP_CM_smRNA_merged_sorted_readsID_2.bed RBP_1_chunk_
split -l 10000000 IgG_IP_CM_smRNA_merged_sorted_readsID_2.bed IgG_chunk_

#for file in RBP_2_chunk_*; do
#    split -l 1000000 "$file" "${file}_second_"
#done

#intersection IP with IgG
for chunk in RBP_2_chunk_*; do
    NAME=$(basename $chunk)
    sbatch --time=24:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'bedtools intersect -f 0.9 -u -s -sorted -a $chunk -b IgG_chunk_* > intersect_$chunk'"
done

#keep only read ID information for intersecting reads
cut -f4 intersect_RBP_1_chunk_* | uniq > intersecting_reads_RBP_1.bed
cut -f4 intersect_RBP_2_chunk_* | uniq > intersecting_reads_RBP_2.bed

#remove intersecting reads
awk '{print $4}' RBP_1_IP_CM_smRNA_merged_sorted_readsID_2.bed > all_IDs_RBP_1.txt
sbatch --time=4:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'grep -Fvx -f intersecting_reads_RBP_1.bed all_IDs_RBP_1.txt > keep_reads_RBP_1.txt'"

sbatch --time=4:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'samtools view -h -N keep_reads_RBP_1.txt RBP_1_IP_CM_smRNA_merged_sorted.bam -o filtered_RBP_1_IP_CM_smRNA_merged_sorted.bam'"


awk '{print $4}' RBP_2_IP_CM_smRNA_merged_sorted_readsID_2.bed > all_IDs_RBP_2.txt
sbatch --time=4:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'grep -Fvx -f intersecting_reads_RBP_2.bed all_IDs_RBP_2.txt > keep_reads_RBP_2.txt'"

sbatch --time=4:00:00 -n 20 --mem-per-cpu=4096 --tmp=40000  --wrap="bash -c 'samtools view -h -N keep_reads_RBP_2.txt RBP_2_IP_CM_smRNA_merged_sorted.bam -o filtered_RBP_2_IP_CM_smRNA_merged_sorted.bam'"


#create bed files
bedtools bamtobed -i filtered_RBP_1_IP_CM_smRNA_merged_sorted.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""RBP_1_IP_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) print i, coverage[i]}' | sortBed -i -  > filtered_RBP_1_IP_CM_smRNA_merged_sorted.bed

bedtools bamtobed -i filtered_RBP_2_IP_CM_smRNA_merged_sorted.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""RBP_2_IP_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) print i, coverage[i]}' | sortBed -i -  > filtered_RBP_2_IP_CM_smRNA_merged_sorted.bed

bedtools bamtobed -i WT_CM_input_smRNA_merged_sorted.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""Input_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) print i, coverage[i]}' | sortBed -i -  > WT_CM_input_smRNA_merged_sorted.bed

#master final data
#select fragments present in all replicates
for d in /Volumes/AS_2023/1_exRBP_data/excerpt_output/RBP_IP_1/bam_files/L0724A12*
do
NAME=$(basename $d _endogenousAlignments_genome_Aligned_sorted.bam)
bedtools bamtobed -i $d | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""IP_CM_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 3) print i, coverage[i]}'  | sortBed -i -  > /Volumes/AS_2023/1_exRBP_data/excerpt_output/RBP_IP_1/${NAME}.bed
done

bed_file1="/Volumes/AS_2023/1_exRBP_data/excerpt_output/RBP_IP_1/L0724A10.bed"
bed_file2="/Volumes/AS_2023/1_exRBP_data/excerpt_output/RBP_IP_1/L0724A11.bed"
bed_file3="/Volumes/AS_2023/1_exRBP_data/excerpt_output/RBP_IP_1/L0724A12.bed"


# Define the input BED file (replace with your actual file path)
input_bed="/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/WT_CM_input_smRNA_merged_sorted.bed"

# Define the output file
output_bed="/Users/alessandra/polybox/sequencing/RNA-seq_analysis/RBP_IP_analysis/3-RBP_IP/WT_CM_input_smRNA_merged_sorted_MASTER_a.bed"

# Perform overlap filtering using BEDTools
bedtools multiinter -i $input_bed $bed_file1 $bed_file2 $bed_file3 -s | awk '($4 >= 4)' > reliable_peaks.bed

bedtools intersect -u -sorted -a $input_bed -b reliable_peaks.bed  > $output_bed



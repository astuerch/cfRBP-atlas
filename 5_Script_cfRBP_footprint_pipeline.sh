# Prepare data before analysis 
######################################################################################
#1. get genome aligned files from excerpt output
######################################################################################
#collect all the endogenousAlignments_genome_Aligned.out files!!!!!
#sort BAM files
cd /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis

mkdir ./WT
for d in ./WT/*.bam
do
    NAME=$(basename $d .bam)
    samtools sort $d > ./WT/${NAME}_sorted.bam
done

######################################################################################
#2. generated WT BED files for the pipeline (do for each batch!!!!!)
######################################################################################
#merge bam replicates (CM small, Cell small, Cell long)
samtools merge WT_CM_smRNA_AGObatch_merged.bam ./WT/*_CM_sorted.bam
samtools merge WT_Cell_smRNA_AGObatch_merged.bam ./WT/*small_Cell_sorted.bam
samtools merge WT_Cell_longRNA_AGObatch_merged.bam ./WT/*long_Cell_sorted.bam

### simplify files by keeping only unique entry coordinates (and not reads name), and adding the coverage.
bedtools bamtobed -i WT_Cell_smRNA_AGObatch_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""WT_AGObatch_Cell_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 4) print i, coverage[i]}'  | sortBed -i -  > WT_Cell_smRNA_AGObatch_merged_filtered.bed

bedtools bamtobed -i WT_CM_smRNA_AGObatch_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""WT_AGObatch_CM_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 4) print i, coverage[i]}'  | sortBed -i -  > WT_CM_smRNA_AGObatch_merged_filtered.bed

bedtools bamtobed -i WT_Cell_longRNA_AGObatch_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""WT_AGObatch_Cell_longRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) print i, coverage[i]}' | sortBed -i -  > WT_Cell_longRNA_AGObatch_merged_filtered.bed


###prepare WT Cell data
#cluster reads with span of 10bp
bedtools cluster -i WT_Cell_longRNA_AGObatch_merged.bed -s -d 10 > WT_Cell_longRNA_AGObatch_merged_filtered_cluster.txt
bedtools cluster -i WT_Cell_smRNA_AGObatch_merged.bed -s -d 10 > WT_Cell_smRNA_AGObatch_merged_filtered_cluster.txt

#clean file: it combines together overlapping regions, filters for more than 3 reads
Rscript filtered_selected_peaks.R WT_Cell_longRNA_AGObatch_merged_filtered_cluster.txt WT_Cell_longRNA_AGObatch_merged_filtered.bed > WT_Cell_longRNA_AGObatch_merged_filtered.log 2>&1
Rscript filtered_selected_peaks.R WT_Cell_smRNA_AGObatch_merged_filtered_cluster.txt WT_Cell_smRNA_AGObatch_merged_filtered.bed > WT_Cell_smRNA_AGObatch_merged_filtered.log 2>&1

#sort files
sortBed -i WT_Cell_longRNA_AGObatch_merged_filtered.bed > WT_Cell_longRNA_AGObatch_merged_filtered_sorted.bed
sortBed -i WT_Cell_smRNA_AGObatch_merged_filtered.bed > WT_Cell_smRNA_AGObatch_merged_filtered_sorted.bed

rm WT_Cell_longRNA_AGObatch_merged_filtered_cluster.txt
rm WT_Cell_longRNA_AGObatch_merged_filtered.bed
rm WT_Cell_smRNA_AGObatch_merged_filtered_cluster.txt
rm WT_Cell_smRNA_AGObatch_merged_filtered.bed

#normalize coverage cell data using RPM - needed later to filter for cfRNA fragments that are not affected into cells
totalReadsWT1=$(awk '{sum+=$7} END {print sum}' WT_Cell_longRNA_AGObatch_merged_filtered_sorted.bed)
totalReadsWT2=$(awk '{sum+=$7} END {print sum}' WT_Cell_smRNA_AGObatch_merged_filtered_sorted.bed)

awk -v totalReadsWT="$totalReadsWT1" 'BEGIN{OFS=FS="\t"}{rpmWT = ($7 / totalReadsWT) * 1000000; print $1, $2, $3, $4, $5, $6, rpmWT}' WT_Cell_longRNA_AGObatch_merged_filtered_sorted.bed > WT_Cell_longRNA_AGObatch_merged_filtered_sorted_rpm.bed
awk -v totalReadsWT="$totalReadsWT2" 'BEGIN{OFS=FS="\t"}{rpmWT = ($7 / totalReadsWT) * 1000000; print $1, $2, $3, $4, $5, $6, rpmWT}' WT_Cell_smRNA_AGObatch_merged_filtered_sorted.bed > WT_Cell_smRNA_AGObatch_merged_filtered_sorted_rpm.bed



###prepare WT CM smRNA files –> here we filter for reliable peaks (2 out of 3 replicates) since those are our final source for cfRBP-RNA footprints

#reliable peaks
#module load stack/2024-06 samtools/1.17 bedtools2/2.31.0
bed_file1="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A25_WT_CM_smRNA_AGObatch.bed"
bed_file2="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A26_WT_CM_smRNA_AGObatch.bed"
bed_file3="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A27_WT_CM_smRNA_AGObatch.bed"

# Define the input BED file (replace with your actual file path)
input_bed="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_AGObatch_merged_filtered.bed"

# Define the output file
output_bed="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_AGObatch_merged_filtered_reliable.bed"

# Perform overlap filtering using BEDTools –> keep only fragments are present in 2 out 3 replicates
bedtools multiinter -i $input_bed $bed_file1 $bed_file2 $bed_file3 -min 0.5 -names input bed_file1 bed_file2 bed_file3 -s | awk '($4 >= 3)' > reliable_peaks.bed

bedtools intersect -u -sorted -a $input_bed -b reliable_peaks.bed  > $output_bed




######################################################################################
#3. prepare RBP KO data (for each RBP independently)
######################################################################################
#start pipeline for each RBP (to change name RBP in pipeline)
cd /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis/
mkdir /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis/RBP_1/01.CM_differential_peaks
mkdir /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis/RBP_1/02.CM_differential_peaks_protected
mkdir /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis/RBP_1/03.fasta_selected_peaks


#sort bam files
#mkdir ./RBP_1
for d in ./RBP_1/*.bam
do
    NAME=$(basename $d .bam)
    samtools sort $d > ./RBP_1/${NAME}_sorted.bam
done


#generate bam merge files with selected replicates, samtools merge with bam files only concatenates the input files without changing their coordinates
samtools merge ./RBP_1/RBP_1_KO_CM_smRNA_merged.bam ./RBP_1/*_CM_sorted.bam
samtools merge ./RBP_1/RBP_1_KO_Cell_smRNA_merged.bam ./RBP_1/*small_Cell_sorted.bam
samtools merge ./RBP_1/RBP_1_KO_Cell_longRNA_merged.bam ./RBP_1/*long_Cell_sorted.bam

### simplify files by keeping only unique entry coordinates (and not reads name), and adding the coverage.
# for smRNA data, a filtering step for more than 4 and 9 reads for cell and CM respectively was added (we are less stringent with cell data because later we discard any RBP-dependent fragment in CM that has some defect in cells)
bedtools bamtobed -i ./RBP_1/RBP_1_KO_Cell_smRNA_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""RBP_1_KO_Cell_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 4) print i, coverage[i]}'  | sortBed -i -  > ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed

bedtools bamtobed -i ./RBP_1/RBP_1_KO_CM_smRNA_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""RBP_1_KO_CM_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 9) print i, coverage[i]}'  | sortBed -i -  > ./RBP_1/RBP_1_KO_CM_smRNA_merged_filtered.bed

bedtools bamtobed -i ./RBP_1/RBP_1_KO_Cell_longRNA_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""RBP_1_KO_Cell_longRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) print i, coverage[i]}' | sortBed -i -  > ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.bed


rm ./RBP_1/RBP_1_KO_Cell_smRNA_merged.bed
rm ./RBP_1/RBP_1_KO_CM_smRNA_merged.bed
rm ./RBP_1/RBP_1_KO_Cell_longRNA_merged.bam


###prepare WT Cell data
#cluster reads with span of 10bp
bedtools cluster -i ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.bed -s -d 10 > ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_cluster.txt
bedtools cluster -i ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed -s -d 10 > ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_cluster.txt

#clean file: it combines together overlapping regions, filters for more than 3 reads
Rscript filtered_selected_peaks.R ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_cluster.txt ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.bed > ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.log 2>&1
Rscript filtered_selected_peaks.R ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_cluster.txt ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed > ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.log 2>&1

#sort files
sortBed -i ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.bed > ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_sorted.bed
sortBed -i ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed > ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_sorted.bed

rm ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_cluster.txt
rm ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered.bed
rm ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_cluster.txt
rm ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed


#normalize coverage cell data using RPM - needed later to filter for cfRNA fragments that are not affected into cells
totalReadsKO1=$(awk '{sum+=$7} END {print sum}' ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_sorted.bed)
totalReadsKO2=$(awk '{sum+=$7} END {print sum}' ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed)

awk -v totalReadsWT="$totalReadsKO1" 'BEGIN{OFS=FS="\t"}{rpmWT = ($7 / totalReadsWT) * 1000000; print $1, $2, $3, $4, $5, $6, rpmWT}' ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_sorted.bed > ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_sorted_rpm.bed
awk -v totalReadsWT="$totalReadsKO2" 'BEGIN{OFS=FS="\t"}{rpmWT = ($7 / totalReadsWT) * 1000000; print $1, $2, $3, $4, $5, $6, rpmWT}' ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered.bed > ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_rpm.bed




######################################################################################
#4. find cfRNA-RBP footprints using RBP KO data (at this stage still batch dependent)
######################################################################################

###find RBP dependent fragments in CM - intersect wt and KO peaks to find disappearing peaks in KO
bedtools intersect -v -s  -sorted -bed -a WT_CM_smRNA_AGObatch_merged_filtered_reliable.bed -b ./RBP_1/RBP_1_KO_CM_smRNA_merged_filtered.bed > ./RBP_1/01.CM_differential_peaks/diff_WT_RBP_1_CM_smRNA.bed

wc -l ./RBP_1/01.CM_differential_peaks/diff_WT_RBP_1_CM_smRNA.bed

###find peaks in cells which biogenesis is not related to the studied RBP
#get peaks with same expression in cell longRNA and smallRNA -> keep overlap filtering stricter because it's cell data, it has to mantain at least 1/4 of the expression in cells (otherwise filtered out as "biogenesis" defect)
bedtools intersect -wb -s -f 0.5 -r -sorted -a WT_Cell_longRNA_AGObatch_merged_filtered_sorted_rpm.bed -b ./RBP_1/RBP_1_KO_Cell_longRNA_merged_filtered_sorted_rpm.bed | awk 'BEGIN{OFS=FS="\t"} $7/4 > $14' | cut -f1-6 | sortBed -i - | uniq  > ./RBP_1/02.CM_differential_peaks_protected/diff_exp_WT_RBP_1_Cell_longRNA.bed

bedtools intersect -wb -s -f 0.5 -r -sorted -a WT_Cell_smRNA_AGObatch_merged_filtered_sorted_rpm.bed -b ./RBP_1/RBP_1_KO_Cell_smRNA_merged_filtered_rpm.bed | awk 'BEGIN{OFS=FS="\t"} $7/4 > $14' | cut -f1-6 | sortBed -i - | uniq > ./RBP_1/02.CM_differential_peaks_protected/diff_exp_WT_RBP_1_Cell_smallRNA.bed

#get final protected peaks for exRBPs with filtering for target expression in cells –> peaks not related to biogenesis but only to export/sort/protection of exRNA
bedtools intersect -v -s -sorted -a ./RBP_1/01.CM_differential_peaks/diff_WT_RBP_1_CM_smRNA.bed -b ./RBP_1/02.CM_differential_peaks_protected/diff_exp_WT_RBP_1_Cell_longRNA.bed  ./RBP_1/02.CM_differential_peaks_protected/diff_exp_WT_RBP_1_Cell_smallRNA.bed > ./RBP_1/02.CM_differential_peaks_protected/RBP_1_CM_protected_peaks.bed



######################################################################################
#5. create a MASTER WT CM file (batch independent) and intersect with cfRNA-RBP footprints to have a final reliable cfRBP-RNA database
######################################################################################

###create master WT CM dataset with only peaks present in all replicates
#sort all bam files for wt cm replicates and merge in an unique file
samtools sort ./WT/L1223A27_endogenousAlignments_genome_Aligned.out.bam > /Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/BAM_files_IGV/L1223A27_WT_CM_sorted.bam

samtools merge WT_CM_smRNA_all_replicates_merged.bam ./BAM_files_IGV/*WT_CM_sorted.bam

# simplify files by keeping only unique entry coordinates (and not reads name) and adding the coverage
bedtools bamtobed -i WT_CM_smRNA_all_replicates_merged.bam | awk 'BEGIN{OFS="\t"} {key=$1"\t"$2"\t"$3"\t""WT_CM_smRNA""\t"".""\t"$6} {coverage[key]+=1} END{for (i in coverage) if(coverage[i] > 4) print i, coverage[i]}'  | sortBed -i -  > WT_CM_smRNA_all_replicates_merged.bed

#select fragments present in all replicates
bed_file1="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0423A9_WT_CM_smRNA_AGObatch.bed"
bed_file2="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0423A41_WT_CM_smRNA_AGObatch.bed"
bed_file3="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A25_WT_CM_smRNA_miRNAbatch.bed"
bed_file4="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A26_WT_CM_smRNA_miRNAbatch.bed"
bed_file5="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723A27_WT_CM_smRNA_miRNAbatch.bed"
bed_file6="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723B24_WT_CM_smRNA_AGObatch.bed"
bed_file7="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0723B29_WT_CM_smRNA_AGObatch.bed"
bed_file8="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L0922B3_WT_CM_smRNA_AGObatch.bed"
bed_file9="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L1223A25_WT_CM_smRNA_GAbatch.bed"
bed_file10="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_rep/L1223A26_WT_CM_smRNA_GAbatch.bed"

# Define the input and output BED file (replace with your actual file path)
input_bed="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_all_replicates_merged.bed"
output_bed="/Volumes/AS_2023/1_exRBP_data/analysis/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_master_fragments.bed"

# Perform overlap filtering using BEDTools
bedtools multiinter -i $input_bed $bed_file1 $bed_file2 $bed_file3 $bed_file4 $bed_file5 $bed_file6 $bed_file7 $bed_file8 $bed_file9 $bed_file10 -names input bed_file1 bed_file2 bed_file3 bed_file4 bed_file5 bed_file6 bed_file7 bed_file8 bed_file9 bed_file10 -s | awk '($4 >= 11)' > reliable_peaks.bed

bedtools intersect -u -sorted -a $input_bed -b reliable_peaks.bed  > $output_bed



###collect footprints from all the rbps tested – protected peaks without annotation
for d in /Volumes/My\ Book/analysis/RNA-seq_analysis/IGV_CLIP_analysis/*
do
NAME=$(basename $d)
bedtools merge -i $d/02.CM_differential_peaks_protected/*_CM_protected_peaks.bed -s -c 4,6,7 -o distinct,distinct,sum > $d/02.CM_differential_peaks_protected/${NAME}_CM_protected_peaks_merged.bed
awk -v RBP="$NAME" 'BEGIN{OFS=FS="\t"}{print $0, RBP}'  $d/02.CM_differential_peaks_protected/${NAME}_CM_protected_peaks_merged.bed >> footprint_allRBPs_merged.bed
done

sortBed -i footprint_allRBPs_merged.bed > footprint_allRBPs_sorted.bed

#intersect with MASTER WT CM
bedtools intersect -wa -wb -f 1 -s -sorted -bed -a WT_CM_smRNA_master_fragments.bed -b footprint_allRBPs_sorted.bed | awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $15, $6, $7}' | uniq >  footprint_allRBPs_master_fragments_sorted.bed

#in R
footprint_allRBPs_master_fragments_sorted1 <- read.delim("~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/footprint_allRBPs_master_fragments_sorted.bed", header=FALSE)
final <- unique(footprint_allRBPs_master_fragments_sorted1)

write.table(final, "~/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/footprint_allRBPs_master_fragments_sorted.bed", quote = F, col.names = F, row.names = F, sep = "\t")
#THIS IS THE FINAL FOOTPRINT FILE: coordinates (from MASTER WT CM, with read ID) of cfRBP-RNA footprints

#add annotations
bedtools intersect -s -f 1 -a footprint_allRBPs_master_fragments_sorted.bed \
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
    print $1, $2, $3, $4, $5, $6, $7, gene_name, gene_type, $10 ; }' > footprint_allRBPs_master_fragments_sorted_annotation2.txt

Rscript footprint_annotation.R footprint_allRBPs_master_fragments_sorted_annotation2.txt footprint_allRBPs_master_fragments_sorted_annotation2.txt > footprint_allRBPs_master_fragments_sorted_annotation2.log 2>&1



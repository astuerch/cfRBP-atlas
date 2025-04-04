########## PREPROCESSING OF RAW DATA FOR ANALYSIS

# UPOP data is only available as bam files, not called peaks

#1) merge bam replicates
samtools merge UPOP_merge.bam *sorted.bam #samtools merge with bam files only concatenates the input files without changing their coordinates

#sort bam file
samtools sort -o UPOP_merged_sorted.bam UPOP_merge.bam 
#index
samtools index UPOP_merged_sorted.bam


#2) interesect merged bam file with called peaks to retrieve peak coverage
bedtools genomecov -ibam UPOP_merged_sorted.bam -g /Users/alessandra/homer/data/genomes/hg38/genome.fa  -bg > UPOP_merged_sorted.bedgraph #transform from bam to bedgraph

bedtools intersect -a UPOP_merged_sorted.bedgraph -b UPOP_peaks_sorted.bed  > UPOP_merged_sorted_peaks.bedgraph #intersect

#combine coverage for same peak
bedtools merge -d 1 -c 4 -o sum -i UPOP_merged_sorted_peaks_sorted.bedgraph > UPOP_merged_sorted_union_peaks.bed

#3) get read length
awk '{print $3 - $2}' UPOP_merged_sorted_union_peaks.bed  | sort | uniq -c  > peak_lengths_summary_upop.txt


#ePRINT: possible to download peak coverage

#1) combine coverage of reads under same peak
bedtools cluster -i ep293TsiNeg_merged_peaks_sorted.bed -s -d 10 > ep293TsiNeg_merged_peaks_sorted_cluster.txt

#2) clean file in R
Rscript filtered_selected_peaks.R WT_Cell_longRNA_GAbatch_merged_filtered_cluster.txt ep293TsiNeg_merged_peaks_sorted_filtered.bed > ep293TsiNeg_merged_peaks_sorted_filtered.log 2>&1


sortBed -i ep293TsiNeg_merged_peaks_sorted_filtered.bed > ep293TsiNeg_merged_peaks_sorted_filtered_sorted.bed

rm ep293TsiNeg_merged_peaks_sorted_cluster.txt
rm ep293TsiNeg_merged_peaks_sorted_filtered.bed

#3) get read length
awk '{print $3 - $2}' /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig2/ep293TsiNeg_merged_peaks_sorted_filtered_sorted.bed  | sort | uniq -c | awk '$2 >= 10 && $2 <= 80' > /Users/alessandra/polybox/Shared/cfRBP_atlas_paper/data/from_fig2/peak_lengths_summary_eprint.txt


#4) combine wt cell, wt cm, upop and eprint read lengths in one file: peak_lengths_summary.txt


#### FIGURE 2B-E: common fragments between WT CM and UPOP/ePRINT datasets

#find overlaps peaks UPOP (Fig 2B)
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/input_CLIP_comparison/UPOP-seq

#prepare UPOP data
Rscript /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/merge_fragments_footprint.R UPOP_peaks_sorted.bed UPOP_peaks_sorted_merged.bed
sortBed -i UPOP_peaks_sorted_merged.bed > UPOP_peaks_sorted_merged_sorted.bed
rm UPOP_peaks_sorted_merged.bed

#intersect data
bedtools intersect -v -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed -b UPOP_peaks_sorted_merged_sorted.bed | wc -l #unique_to_A.bed 200549

bedtools intersect -v -s -a UPOP_peaks_sorted_merged_sorted.bed -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed | wc -l #unique_to_B.bed 293277
bedtools intersect -u -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed -b UPOP_peaks_sorted_merged_sorted.bed | wc -l #intersection 30156 (13%)


####only for mRNA (Fig 2D)
# annotate peaks
bedtools intersect -s -f 1 -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed \
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
}' > /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_annotation.txt

#clean file
Rscript /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/footprint_annotation.R /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_annotation.txt /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_annotation.txt > /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_annotation.log 2>&1

#filter for mrna reads
awk '$9 == "protein_coding"' /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_annotation.txt > /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed

#intersect
bedtools intersect -v -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed -b UPOP_peaks_sorted_merged_sorted.bed | wc -l #unique_to_A.bed 58725 (tot 75856)

bedtools intersect -v -s -a UPOP_peaks_sorted_merged_sorted.bed -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed | wc -l #unique_to_B.bed 297145
bedtools intersect -u -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed -b UPOP_peaks_sorted_merged_sorted.bed | wc -l #intersection 17131 (22.6%)



#find overlaps peaks ePRINT (Fig 2C)
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/input_CLIP_comparison/ePRINT/

for d in *Aligned.out.bam
do
    NAME=$(basename $d .bam)
    samtools sort $d > ${NAME}_sorted.bam
done

samtools merge ep293TsiNeg_merged.bam *sorted.bam #taken 3 WT libraries (L1023A49, A50, A73) to be comparable with longRNA seq from RBP KOs since here coverage could be important
bedtools bamtobed -i ep293TsiNeg_merged.bam | awk 'BEGIN{OFS="\t"} {print $1"\t"$2"\t"$3"\t""ep293TsiNeg""\t"".""\t"$6}' | sort -k1,1 -k2,2n | uniq > ep293TsiNeg_merged.bed

#convert coordinates bedgraph peaks from hg19 to hg38

bedtools intersect -a ep293TsiNeg_merged.bed -b GSM7186916_ep_293T_siNeg_R1_hg38.bed GSM7186917_ep_293T_siNeg_R2_hg38.bed > ep293TsiNeg_merged_peaks.bed

sort -k1,1 -k2,2n ep293TsiNeg_merged_peaks.bed | uniq > ep293TsiNeg_merged_peaks_sorted.bed

#clean data
cd /Users/alessandra/polybox/sequencing/RNA-seq_analysis/input_CLIP_comparison/ePRINT

awk '{print $0 "\t0"}' ep293TsiNeg_merged_peaks_sorted1.bed > ep293TsiNeg_merged_peaks_sorted.bed

Rscript /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/merge_fragments_footprint.R ep293TsiNeg_merged_peaks_sorted.bed ep293TsiNeg_merged_peaks_sorted_merged.bed

sortBed -i ep293TsiNeg_merged_peaks_sorted_merged.bed > ep293TsiNeg_merged_peaks_merged_sorted.bed
rm ep293TsiNeg_merged_peaks_sorted_merged.bed

#intersect with cm data
bedtools intersect -v -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed -b ep293TsiNeg_merged_peaks_merged_sorted.bed | wc -l #unique_to_A.bed 184335

bedtools intersect -v -s -a ep293TsiNeg_merged_peaks_merged_sorted.bed -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed | wc -l #unique_to_B.bed 9999148
bedtools intersect -u -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted.bed -b ep293TsiNeg_merged_peaks_merged_sorted.bed | wc -l #intersection 46370 (20%)


####only for mRNA (Fig 2E)
bedtools intersect -v -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed -b ep293TsiNeg_merged_peaks_merged_sorted.bed | wc -l #unique_to_A.bed 48989 (tot 75856)

bedtools intersect -v -s -a ep293TsiNeg_merged_peaks_merged_sorted.bed -b /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed | wc -l #unique_to_B.bed 10029033
bedtools intersect -u -s -a /Users/alessandra/polybox/sequencing/RNA-seq_analysis/IGV_CLIP_analysis/WT_CM_smRNA_allbatches_merged_filtered_reliable_sorted_mRNA.bed -b ep293TsiNeg_merged_peaks_merged_sorted.bed | wc -l #intersection 26867 (35.4%)



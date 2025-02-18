#!/bin/sh

#  read-lenght_script.sh
#  
#
#  Created by Alessandra Stürchler on 09.01.25.
#  
#biotype read length for WT cell and cm samples
module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.9


#1) filter bam files for unique mapping reads and extract annotation and read length
for d in /cluster/work/demello/Project1-exRBP/RBP_KO_Seq/excerpt_output_WT_12bp/*/
do
NAME=$(basename $d)
sbatch --time=4:00:00 -n 20 --mem-per-cpu=2048 --tmp=40000  --wrap="samtools view $d/endogenousAlignments_genomeMapped_transcriptome_Aligned.out.sorted.bam | grep -E '(NH:i:1|^@)' | awk '{print \$3, length(\$10)}' > ${NAME}_unique_mapped.txt"
done

#2) combine cell and cm replicates, divide by biotypes
cat CM* > WT_CM_read_length.txt
wc -l WT_CM_read_length.txt # 34775265

cat Cell* > WT_Cell_read_length.txt
wc -l WT_Cell_read_length.txt # 38719915

grep -E 'protein_coding' WT_Cell_read_length.txt > WT_Cell_mRNA_read_length.txt
grep -E 'tRNA' WT_Cell_read_length.txt > WT_Cell_tRNA_read_length.txt
grep -E 'misc_RNA' WT_Cell_read_length.txt > WT_Cell_miscRNA_read_length.txt
grep -E 'snoRNA' WT_Cell_read_length.txt > WT_Cell_snoRNA_read_length.txt
grep -E 'snRNA' WT_Cell_read_length.txt > WT_Cell_snRNA_read_length.txt
grep -E 'piRNA' WT_Cell_read_length.txt > WT_Cell_piRNA_read_length.txt
grep -E 'miRNA' WT_Cell_read_length.txt > WT_Cell_miRNA_read_length.txt
grep -vE 'protein_coding|tRNA|misc_RNA|snoRNA|snRNA|piRNA|miRNA' WT_Cell_read_length.txt > WT_Cell_other_read_length.txt

grep -E 'protein_coding' WT_CM_read_length.txt > WT_CM_mRNA_read_length.txt
grep -E 'tRNA' WT_CM_read_length.txt > WT_CM_tRNA_read_length.txt
grep -E 'misc_RNA' WT_CM_read_length.txt > WT_CM_miscRNA_read_length.txt
grep -E 'snoRNA' WT_CM_read_length.txt > WT_CM_snoRNA_read_length.txt
grep -E 'snRNA' WT_CM_read_length.txt > WT_CM_snRNA_read_length.txt
grep -E 'piRNA' WT_CM_read_length.txt > WT_CM_piRNA_read_length.txt
grep -E 'miRNA' WT_CM_read_length.txt > WT_CM_miRNA_read_length.txt
grep -vE 'protein_coding|tRNA|misc_RNA|snoRNA|snRNA|piRNA|miRNA' WT_CM_read_length.txt > WT_CM_other_read_length.txt


#3) count unique read lengths for each biotype
for d in WT_CM*.txt
do
NAME=$(basename $d .txt)
awk '{print $2}' $d | sort | uniq -c > ${NAME}_summary.txt
done


#!/bin/sh

#
#  
#
#  Created by Alessandra Stürchler on 03.11.24.
#  

#import annotated file as df

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 2) {
  stop("Usage: Rscript merge_fragments_footprint.R <input_file> <output_file>")
}

# Assign arguments to variables
input_file <- args[1]
output_file <- args[2]

library(dplyr)

df <- read.delim(input_file, header=FALSE)

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
unique_df2 <- unique_df2 %>%
  mutate(across(where(~ methods::is(.x, "Rle")), as.vector))

#get the best annotation for each peak
#selected using the following priority: CDS exon > 3′ UTR > 5′ UTR > protein-coding gene intron > noncod- ing RNA exon > noncoding RNA intron > intergenic.
result_df <- unique_df2 %>%
  arrange(V1, V2, V3) %>%
  group_by(V1, V2, V3, V4, V5, V6, V7) %>%
  mutate(
    biotype = case_when(
      V9 %in% c("miRNA", "tRNA", "snRNA", "snoRNA", "protein_coding", "misc_RNA", "piRNA") ~ as.character(V9),
      TRUE ~ "others"  # Assign "others" to any biotype not listed above
    ),
    rank_V9 = match(biotype, c("miRNA", "tRNA", "snRNA", "snoRNA", "protein_coding", "misc_RNA", "piRNA", "others")),
    rank_V10 = match(V10, c("ncRNA", "three_prime_utr", "five_prime_utr", "start_codon", "stop_codon", "exon", "transcript", "CDS","gene"))
  ) %>%
  arrange(rank_V9, rank_V10) %>%
  filter_all(all_vars(!is.na(.))) %>%
  slice(1) %>%  # Take the first row after sorting by both ranks
  select(-rank_V9, -rank_V10) %>%  ungroup()

write.table(result_df, output_file, quote = F, col.names = F, row.names = F, sep = "\t")
            
            

#tot fragments footprint for each RBP
rbp_footprint <- result_df %>%
  group_by(V5) %>%
  summarise(unique_peaks = n_distinct(V1,V2,V3))

#biotype fragments footprint for each RBP
result <- result_df %>%
  group_by(V5,biotype) %>%
  summarise(unique_peaks = n_distinct(V1,V2,V3))

#tot fragments
df1 <- unique(result_df[, c(1,2,3,5,8, 10,11)]) %>% select(-V1,-V2,-V3)
colnames(df1) <- c("cfRBP","gene", "feature", "biotype")

result1 <- df1 %>%
  group_by(biotype) %>%
  summarise(unique_genes = n_distinct(gene))

#tot number of features
result2 <- df1 %>%
  count(feature, name = "feature_frequency")
  

cat("--- Number of fragments per RBP:\n")
print(rbp_footprint)

cat("--- Number of fragments biotype per RBP:\n")
print(result)


cat("--- total number of annotated fragment per biotype", sum(result$unique_peaks), "\n")
cat("--- total number of genes ", sum(result1$unique_genes), "\n")

cat("--- the different features:\n")
print(result2)


cat("Work complete. Output written to", output_file, "\n")


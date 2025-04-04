# summarizing selected peaks in clusters
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
library(tidyr)

data3 <- read.delim(input_file, header=FALSE)

colnames(data3) <- c("chr", "start", "end", "sample", ".", "strand", "coverage", "cluster")

result7 <- data3 %>%
  group_by(cluster) %>%
  summarise(start = min(start),
            end = max(end),
            coverage = sum(coverage),
            across(everything(), ~ first(.))) %>% # Retain all other columns in their original order
  ungroup() %>%
  filter(coverage >= 3)


result7 <- result7[, names(data3)]

write.table(result7[,-8], output_file, quote = F, col.names = F, row.names = F, sep = "\t")

cat("Work complete. Output written to", output_file, "\n")

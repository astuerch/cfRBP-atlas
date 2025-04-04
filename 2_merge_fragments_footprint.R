
library(data.table)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 2) {
  stop("Usage: Rscript merge_fragments_footprint.R <input_file> <output_file>")
}

# Assign arguments to variables
input_file <- args[1]
output_file <- args[2]


# Load the data
dt <- fread(input_file, header = FALSE)
setnames(dt, c("chrom", "start", "end", "name", "dot", "strand", "score"))

# Data Validation

# a. Check for missing 'start' or 'end'
missing_start <- dt[is.na(start)]
missing_end <- dt[is.na(end)]

if (nrow(missing_start) > 0) {
  cat("Rows with missing 'start' positions:\n")
  print(missing_start)
}

if (nrow(missing_end) > 0) {
  cat("Rows with missing 'end' positions:\n")
  print(missing_end)
}

# Remove rows with missing 'start' or 'end'
dt_clean <- dt[!is.na(start) & !is.na(end)]

# b. Check for 'start' > 'end'
invalid_coords <- dt_clean[start > end]

if (nrow(invalid_coords) > 0) {
  cat("Rows with 'start' > 'end':\n")
  print(invalid_coords)
  
  # Remove invalid rows
  dt_clean <- dt_clean[start <= end]
}

# c. Check for non-integer 'start' or 'end'
non_integer_start <- dt_clean[!grepl("^\\d+$", start)]
non_integer_end <- dt_clean[!grepl("^\\d+$", end)]

if (nrow(non_integer_start) > 0) {
  cat("Rows with non-integer 'start' positions:\n")
  print(non_integer_start)
}

if (nrow(non_integer_end) > 0) {
  cat("Rows with non-integer 'end' positions:\n")
  print(non_integer_end)
  
  # Remove non-integer rows
  dt_clean <- dt_clean[grepl("^\\d+$", start) & grepl("^\\d+$", end)]
}

# Convert 'start' and 'end' to integers
dt_clean[, start := as.integer(start)]
dt_clean[, end := as.integer(end)]


# Sort the data by chromosome, strand, start, and end
setorder(dt_clean, chrom, strand, start, end)

# Define the are_similar function
are_similar <- function(start1, end1, start2, end2) {
  # Ensure all inputs are numeric and not NA
  if (any(is.na(c(start1, end1, start2, end2)))) {
    return(FALSE)
  }
  
  # Calculate proximity
  proximity <- (abs(start1 - start2) <= 5) && (abs(end1 - end2) <= 5)
  
  return(proximity)
}

are_similar_vectorized <- Vectorize(are_similar)


# Assign group IDs using data.table's vectorized operations
dt_clean[, similar_to_prev := are_similar_vectorized(
  shift(start, type = "lag"),
  shift(end, type = "lag"),
  start,
  end
), by = .(chrom, strand)]

# Handle NA in 'similar_to_prev' by setting them to FALSE
dt_clean[is.na(similar_to_prev), similar_to_prev := FALSE]

# Assign group IDs: increment group when 'similar_to_prev' is FALSE
dt_clean[, group := cumsum(!similar_to_prev), by = .(chrom, strand)]


# Merge the fragments based on groups
final_dt <- dt_clean[, .(
  chrom = first(chrom),
  start = min(start),
  end = max(end),
  name = first(name),      # Modify if needed
  dot = first(dot),
  strand = first(strand),
  score = sum(score)       # Modify aggregation as needed
), by = .(chrom, strand, group)]


# Write the merged data to a file
fwrite(final_dt[,4:10], file = output_file, sep = "\t", col.names = FALSE)

cat("Merging complete. Output written to", output_file, "\n")


  

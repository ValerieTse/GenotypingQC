#!/usr/bin/env Rscript

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript call_rate.r <smiss_file> <vmiss_file>")
}

smiss_file <- args[1]
vmiss_file <- args[2]

# ----- Read and Process .smiss File -----
smiss <- read.table(smiss_file, header = FALSE)
colnames(smiss) <- c("FID", "IID", "PHENO1", "MISSING_CT", "OBS_CT", "F_MISS")
smiss$call_rate <- 1 - smiss$F_MISS

# Print the maximum missing rate for samples
max_missing_sample <- max(smiss$F_MISS)
cat("Maximum missing rate for samples:", max_missing_sample, "\n")

# ----- Read and Process .vmiss File -----
vmiss <- read.table(vmiss_file, header = FALSE)
colnames(vmiss) <- c("CHROM", "ID", "MISSING_CT", "OBS_CT", "F_MISS")
vmiss$call_rate <- 1 - vmiss$F_MISS

# Print the maximum missing rate for each chromosome
max_missing_by_chrom <- aggregate(F_MISS ~ CHROM, data = vmiss, FUN = max)
cat("Maximum missing rate for each chromosome:\n")
print(max_missing_by_chrom)


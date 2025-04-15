#!/usr/bin/env Rscript

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript het_rates.r <het_file> <prefix> <Folder>")
}
het_file <- args[1]
prefix <- args[2]
Folder <- args[3]

# Read the PLINK .het file. Adjust the filename as needed.
het_df <- read.table(het_file, header = TRUE)
colnames(het_df) <- c("FID", "IID", "O_HOM", "E_HOM", "N_NM", "F")

# Compute the heterozygosity rate for each sample.
het_df$het_rate <- (het_df$N_NM - het_df$O_HOM) / het_df$N_NM

result <- het_df[, c("FID", "IID", "het_rate")]

write.table(result, paste0(Folder, "/DATA/", prefix,"_het.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)

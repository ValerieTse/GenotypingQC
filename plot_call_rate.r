#!/usr/bin/env Rscript

# ----- Load Required Library -----
library(ggplot2)

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_call_rate.r <smiss_file> <vmiss_file> <prefix> <Folder>")
}

smiss_file <- args[1]
vmiss_file <- args[2]
prefix <- args[3]
Folder <- args[4]

# ----- Read and Process .smiss File -----
smiss <- read.table(smiss_file, header = FALSE)
colnames(smiss) <- c("FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS")
smiss$call_rate <- 1 - smiss$F_MISS

# Print the maximum missing rate for samples
max_missing_sample <- max(smiss$F_MISS)
cat("Maximum missing rate for samples:", max_missing_sample, "\n")

# Create histogram of sample call rate using ggplot2
hist_smiss <- ggplot(smiss, aes(x = call_rate)) +
  geom_histogram(bins = 30, fill = "olivedrab3", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Sample Call Rate",
       x = "Sample Call Rate",
       y = "Frequency") +
  xlim(0.7, 1) +
  theme_bw()

# Save the sample call rate histogram with prefix
ggsave(paste0(Folder,"/OUTPUTS/", prefix, "_hist_sample_call_rate.png"), plot = hist_smiss, width = 6, height = 4, dpi = 300)

# ----- Read and Process .vmiss File -----
vmiss <- read.table(vmiss_file, header = FALSE)
colnames(vmiss) <- c("CHROM", "ID", "MISSING_CT", "OBS_CT", "F_MISS")
vmiss$call_rate <- 1 - vmiss$F_MISS

# Print the maximum missing rate for each chromosome
max_missing_by_chrom <- aggregate(F_MISS ~ CHROM, data = vmiss, FUN = max)
cat("Maximum missing rate for each chromosome:\n")
print(max_missing_by_chrom)

# Create histogram of variant call rate using ggplot2
hist_vmiss <- ggplot(vmiss, aes(x = call_rate)) +
  geom_histogram(bins = 30, fill = "goldenrod1", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Variant Call Rate",
       x = "Variant Call Rate",
       y = "Frequency") +
  xlim(0.7, 1) +
  theme_bw()

# Save the variant call rate histogram with prefix
ggsave(paste0(Folder, "/OUTPUTS/",prefix, "_hist_variant_call_rate.png"), plot = hist_vmiss, width = 6, height = 4, dpi = 300)

# ----- Save Boxplot of Sample Call Rate (Base R Plot) -----
png(paste0(Folder, "/OUTPUTS/", prefix, "_boxplot_sample_call_rate.png"), width = 600, height = 400)
boxplot(smiss$call_rate, main = "Sample Call Rate")
dev.off()

# ----- Save Boxplot of Variant Call Rate by Chromosome (Base R Plot) -----
png(paste0(Folder, "/OUTPUTS/", prefix, "_boxplot_variant_call_rate_by_chrom.png"), width = 600, height = 400)
par(mar = c(7, 4, 4, 2) + 0.1)
boxplot(call_rate ~ CHROM, data = vmiss,
        main = "SNP Call Rate by Chromosome",
        xlab = "Chromosome",
        ylab = "SNP Call Rate",
        col = "lightblue",
        xaxt = "n")

desired_order <- c(as.character(1:22), "X", "Y", "XY", "MT")
axis(1, at = 1:length(desired_order),
     labels = desired_order,
     las = 2, cex.axis = 0.8)
dev.off()

#!/usr/bin/env Rscript

# Load required packages
suppressMessages(library(ggplot2))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript plot_sexcheck.r <sexcheck_file_path> <output_prefix> <Folder>")
}

sexcheck_file <- args[1]
prefix <- args[2]
Folder <- args[3]

# Read the file; adjust the file path if needed.
sex_data <- read.table(sexcheck_file <- args[1], header = TRUE)

hist_plot <- ggplot(sex_data, aes(x = F)) +
                geom_histogram(binwidth = 0.05, fill = "gray", color = "black") +
                geom_vline(xintercept = c(0.2, 0.8), color = "red", linetype = "dashed") +
                geom_vline(xintercept = c(0.4, 0.6), color = "blue", linetype = "dashed") +
                labs(title = "Histogram of F-statistics from Sex Check",
                    x = "F-statistic",
                    y = "Count") +
                theme_bw()

ggsave(paste0(Folder, "/OUTPUTS/", prefix, "_sexcheck_histogram.png"), hist_plot, width = 8, height = 5)

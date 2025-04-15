#!/usr/bin/env Rscript

# Load required packages
suppressMessages(library(dplyr))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript prepare_pop.r <fam_file_path> <pop_file> <output_prefix> <Folder>")
}

fam_file <- args[1]
pop_file <- args[2]
prefix <- args[3]
Folder <- args[4]

df <- read.table(pop_file, header = TRUE, sep = ",")
df$SPop <- ifelse(df$SPop == "Sample","-",df$SPop)

df2 <- read.table(fam_file)
colnames(df2) <- c("FID","IID","V3","V4","V5","V6")

df <- merge(df2, df, by="IID", all.x = TRUE)
df <- df %>% select(SPop) 

write.table(df, paste0(Folder, "/DATA/", prefix, ".pop"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

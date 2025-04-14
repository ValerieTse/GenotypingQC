#!/usr/bin/env Rscript

# Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript counts.r <plink_prefix>")
}

plink <- args[1]

#  count number of individuals in the study
FAM <-fread(paste0(plink, ".fam"), header=FALSE)
cat("Number of individuals:", dim(FAM)[1], "\n")

# Map Information on number of SNPS count number of genotyped SNPs
map <- fread(paste0(plink, ".bim"), header=FALSE)
cat("Number of SNPs:", dim(map)[1], "\n")


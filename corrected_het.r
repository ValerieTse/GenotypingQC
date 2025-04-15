#!/usr/bin/env Rscript
library(dplyr)

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript corrected_het.r <het_file> <pc_file> <output_prefix> <Folder>")
}
het_file <- args[1]
pc_file <- args[2]
prefix <- args[3]
Folder <- args[4]

het_df <- read.table(het_file, header = TRUE)
pc_file <- data.table::fread(pc_file) 
df <- merge(het_df, pc_file, by = "IID", all.x = TRUE)

## 2. Fit a linear model with all linear, quadratic, 
##    and pairwise interaction terms for PC1–PC4.
##    In R, (PC1 + PC2 + PC3 + PC4)^2 will include:
##      - Linear terms: PC1 + PC2 + PC3 + PC4
##      - Quadratic terms: I(PC1^2), I(PC2^2), I(PC3^2), I(PC4^2)
##      - Pairwise interactions: PC1:PC2, PC1:PC3, PC1:PC4, PC2:PC3, PC2:PC4, PC3:PC4
##
fit <- lm(het_rate ~ (PC1 + PC2 + PC3 + PC4)^2, data = df)

##
## 3. Check the results
##
summary(fit)

##
## 4. Extract the fitted ancestry-corrected heterozygosity
##    The fitted values are effectively h_0 + β(x) for each individual.
##
df$het_corrected <- fit$fitted.values



df %>% select(IID, het_corrected) %>% write.csv(paste0(Folder, "/DATA/",prefix,"_corrected_het.csv"), row.names = FALSE)

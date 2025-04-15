#!/usr/bin/env Rscript

# ----- Load Required Library -----
library(dplyr)

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript sexcheck_chrY.r <smiss_file> <sex_info> <prefix> <Folder>")
}

smiss_file <- args[1]
sex_info <- args[2]
prefix <- args[3]
Folder <- args[4]

# ----- Read and Process .smiss File -----
smiss <- read.table(smiss_file, header = FALSE)
colnames(smiss) <- c("FID", "IID", "PHENO1", "MISSING_CT", "OBS_CT", "F_MISS")
smiss <- smiss %>% select(IID, F_MISS)

sex_info <- read.table(sex_info, header = FALSE)
colnames(sex_info) <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
sex_info <- sex_info %>% select(IID, PEDSEX, SNPSEX, F)

df <- merge(smiss, sex_info , by = "IID", all = TRUE) %>%
      mutate(STATUS = ifelse((PEDSEX == 2 & is.na(F_MISS) == FALSE), 'PROBLEM','NORMAL'))

write.table(df, file = paste0(Folder,"/DATA/", prefix,"_sexcheck_Ymiss.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



#!/usr/bin/env Rscript

# ----- Load Required Library -----
library(dplyr)
library(data.table)

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript step1_update_info_extract.r <FAM_file> <updated_info_file>")
}

fam_file <- args[1]
updated_info_file <- args[2]

# ----- Loading data -----
update.info <- fread(updated_info_file, header = TRUE,  fill = TRUE)
old.fam <- fread(fam_file)
colnames(old.fam) <- c("FID", "IID", "PID", "MID", "Sex", "P")

# ----- Generate the reference data -----
update_ids1 <- old.fam %>% select(FID, IID)
update_ids2 <- update.info %>% select(sfid, spid)
update_ids <- merge(update_ids1, update_ids2, by.x = "IID", by.y = "spid", all.x = TRUE)
sum(is.na(update_ids$sfid))
update_ids$spid <- update_ids$IID
update_ids <- update_ids %>% select(FID, IID, sfid, spid)

# ----- Write to a file for PLINK -----
write.table(update_ids, file = file.path("/Release_Notes/update_id.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

pheno1 <- update.info %>% select(spid, asd) %>% mutate(asd = ifelse(asd == "True", 1, ifelse(asd == "False", 2, 0)))
pheno2 <- update_ids %>% select(sfid, spid)
pheno <- merge(pheno2, pheno1, by.x = "spid", by.y = "spid", all.x = TRUE)
sum(is.na(pheno$asd))
pheno <- pheno %>% select(sfid, spid, asd)

write.table(pheno, file = file.path("/Release_Notes/update_pheno.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

update_sex1 <- update.info %>% select(spid, sex) %>% mutate(sex = ifelse(sex == "Male", 1, ifelse(sex == "Female", 2, 0)))
update_sex2 <- update_ids %>% select(sfid, spid)
update_sex <- merge(update_sex2, update_sex1, by.x = "spid", by.y = "spid", all.x = TRUE)
sum(is.na(update_sex$sex))
update_sex <- update_sex %>% select(sfid, spid, sex)

write.table(update_sex, file = file.path("/Release_Notes/update_sex.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#!/usr/bin/env Rscript

# ----- Load Required Library -----
library(illuminaio)
library(dplyr)
library(purrr)
library(furrr)
library(data.table)
library(parallelly)  # For CPU limit detection
library(tidyr)
library(dplyr)
library(ggplot2)

# ----- Parse Command Line Arguments -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_call_rate.r <idat_dir> <array_manifest_file> <problematics_list> <sexcheck_file> <prefix> <Folder>")
}

idat_dir <- args[1]
manifest_path <- args[2]
problematics_list <- args[3]
sexcheck_file <- args[4]
prefix <- args[5]
Folder <- args[6]

problematics <- read.table(problematics_list)
sexcheck <- read.table(sexcheck_file, header = TRUE) %>% filter(IID %in% summary_df$SampleID) %>% select(IID, PEDSEX, SNPSEX, X_het = F)

# ---- 1. Set safe parallel workers ----
# Detect available CPU resources
available_cores <- availableCores(constraints = "multicore")  # Respects system limits
safe_workers <- max(1, floor(available_cores * 0.8))  # Use 80% of available cores

cat(sprintf(
  "System reports %d available cores. Using %d parallel workers.\n",
  available_cores,
  safe_workers
))

# Set parallel plan
plan(multisession, workers = safe_workers)

# Paths (UPDATE THESE)
# idat_dir <- "/project/cshu0237_1136/Data/SPARK/iWES_v2/genotypes/genotyping_array/idat/"
# manifest_path <- "/project/cshu0237_1136/Users/jingyux/SPARK_GWAS_QC/Release_Notes/GSA-24v2-0_A2.csv"
# output_dir <- "/project/cshu0237_1136/Users/jingyux/SPARK_GenoQC/DATA"
# problematics <- read.table("/project/cshu0237_1136/Users/jingyux/SPARK_GWAS_QC/Array_v2/DATA/step4_problematic_list.txt")

# Create output directory
# if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 1. Load manifest with proper typing ----
lines <- readLines(manifest_path)
if (length(lines) >= 7) {
  lines <- lines[-(1:7)]
}
idx <- grep("\\[Controls\\]", lines)
if(length(idx) > 0) {
  lines <- lines[seq_len(idx[1] - 1)]
}
manifest <- read.csv(textConnection(paste(lines, collapse="\n")),
               stringsAsFactors = FALSE, header = TRUE)
manifest <- manifest %>%
  transmute(
    IlmnID,
    AddressA_ID = as.character(AddressA_ID),  # Green channel
    AddressB_ID = as.character(AddressB_ID),  # Red channel
    Chr = as.character(Chr),
    MapInfo
  ) %>%
  filter(Chr %in% c(as.character(1:22), "X", "Y", "23", "24")) %>%  # Keep autosomes + XY
  mutate(Chr = case_when(  # Standardize chromosome notation
    Chr == "23" ~ "X",
    Chr == "24" ~ "Y",
    TRUE ~ Chr
  ))

manifest_long <- pivot_longer(manifest, cols=c("AddressA_ID", "AddressB_ID"), names_to = "col", values_to = "Address")

# ---- 2. Find all IDAT pairs recursively ----
idat_pairs <- list.files(idat_dir, pattern = "_\\w+\\.idat$", 
                        recursive = TRUE, full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    SampleID = gsub("_(Red|Grn)\\.idat$", "", basename(Path)),
    Channel = ifelse(grepl("_Red\\.idat$", Path), "Red", "Green")
  ) %>%
  group_by(SampleID) %>%
  filter(n() == 2) %>%
  summarise(
    Paths = list(setNames(Path, Channel)),
    .groups = "drop"
  )


idat_pairs_filtered <- idat_pairs %>%
  filter(SampleID %in% problematics$V2)

# ---- 3. Processing function with X/Y intensity calculation ----
process_sample <- function(sample_paths, manifest) {
  tryCatch({
    # Read IDAT files
    signals <- map_dfr(names(sample_paths), ~ {
      data <- readIDAT(sample_paths[.x])
      tibble(
        Address = names(data$Quants[, "Mean"]),
        Signal = data$Quants[, "Mean"],
        Channel = .x
      )
    })
    
    # Mapping
    mapped_data <- merge(signals, manifest_long, by = "Address", all.x = TRUE)
    mapped_data_check <- mapped_data %>% filter((Channel == "Green" & col == "AddressA_ID") |(Channel == "Red" & col == "AddressB_ID") )

    # Calculate chromosome-specific metrics
    chrom_stats <- mapped_data_check %>%
      filter(Chr %in% c("X", "Y")) %>%
      group_by(Chr) %>%
      summarise(
        Mean_Intensity = mean(Signal, na.rm = TRUE),
        Median_Intensity = median(Signal, na.rm = TRUE),
        N_Probes = n(),
        .groups = "drop"
      ) %>%
      mutate(SampleID = gsub("_(Red|Grn)\\.idat$", "", basename(sample_paths[1])))
    
    # Full probe data for X/Y
    xy_probes <- mapped_data_check %>%
      filter(Chr %in% c("X", "Y")) %>%
      mutate(SampleID = chrom_stats$SampleID[1])
    
    list(
      SampleID = unique(xy_probes$SampleID),
      Chromosome_Stats = chrom_stats,
      XY_Probes = xy_probes,
      Source_Path = dirname(sample_paths[1])
    )
  }, error = function(e) {
    message("Failed ", sample_paths[1], ": ", e$message)
    return(NULL)
  })
}

# ---- 4. Process all samples in parallel ----
results <- future_map(
  idat_pairs_filtered$Paths, 
  process_sample, 
  manifest = manifest,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
) %>% compact()

saveRDS(results,paste0(Folder, "/DATA/", prefix, "_XY_raw_results.rds"))

# ---- 5. Save results ----
# Summary CSV with X/Y intensities
summary_df <- map_dfr(results, ~ .x$Chromosome_Stats) %>%
  left_join(
    map_dfr(results, ~ tibble(SampleID = .x$SampleID, Subfolder = .x$Source_Path)),
    by = "SampleID"
  ) %>%
  select(SampleID, Subfolder, everything())

fwrite(summary_df, paste0(Folder, "/DATA/", prefix, "_XY_chromosome_summary.csv"))

summary_X <- summary_df %>% filter(Chr == "X") %>% select(IID = SampleID, ChrX=Chr, X_Mean_Intensity = Mean_Intensity, X_N_Probes = N_Probes)
summary_Y <- summary_df %>% filter(Chr == "Y") %>% select(IID = SampleID, ChrY=Chr, Y_Mean_Intensity = Mean_Intensity, Y_N_Probes = N_Probes)
summary_new <- merge(summary_X,summary_Y, by = "IID")
summary_new <- merge(summary_new, sexcheck, by = "IID")

# ---- 6. Plotting ----
p1 <- ggplot(summary_new, aes(x = X_Mean_Intensity, y = Y_Mean_Intensity, color = as.factor(SNPSEX))) +
          geom_point() +
          labs(x = "X intensity (30522 probes)", y = "Y intensity (5689 probes)", color = "Genetic Sex") +
          scale_color_manual(
               values = c("0" = "grey", "1" = "red", "2" = "blue"),
               labels = c("0" = "NA", "1" = "Male", "2" = "Female")) +
          theme_bw()

ggsave(paste0(Folder, "/OUTPUTS/", prefix, "_YX_intesity_SNPSEX.png"), plot = p1, width = 6, height = 5, dpi = 300)

p2 <- ggplot(summary_new, aes(x = X_Mean_Intensity, y = Y_Mean_Intensity, color = as.factor(PEDSEX))) +
          geom_point() +
          labs(x = "X intensity (30522 probes)", y = "Y intensity (5689 probes)", color = "Annotated Sex") +
          scale_color_manual(
               values = c("1" = "red", "2" = "blue"),
               labels = c("1" = "Male", "2" = "Female")) +
          theme_bw()

ggsave(paste0(Folder, "/OUTPUTS/", prefix, "_YX_intesity_PEDSEX.png"), plot = p2, width = 6, height = 5, dpi = 300)

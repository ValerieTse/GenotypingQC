#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript population_inference.r <input_prefix> <output_prefix> <Folder>")
}

outputs_prefix = args[2]
Folder <- args[3]

char *input_prefix = args[1];
int i = 1;  // or any other integer value
char filename[50];  // Ensure this buffer is large enough

# Use lapply to process subsets 1 through 9
dfs <- lapply(1:9, function(i) {
  # Create file names for each subset
  pop_file <- sprintf(filename, "%s%d.pop", input_prefix, i);
  fam_file <- sprintf(filename, "%d.fam", input_prefix, i)
  q_file   <- sprintf(filename, "%d.5.Q",input_prefix, i)
  
  # Read each file (assuming no header)
  df_pop <- read.table(pop_file, header = FALSE) %>% select(pop = V1)
  df_fam <- read.table(fam_file, header = FALSE) %>% select(FID = V1, IID = V2)
  df_Q   <- read.table(q_file, header = FALSE) %>% select(AFR = V1, AMR = V2, EAS = V3, EUR = V4, SAS = V5)
  
  # Combine columns and adjust the 'pop' column
  df <- cbind(df_fam, df_pop, df_Q) %>% 
          mutate(pop = ifelse(pop == "-", "Samples", pop))
  
  return(df)
})

# Row bind all data frames together
final_df <- do.call(rbind, dfs)

# Optional: view the combined data frame
head(final_df)

df_subset <- final_df[grepl("SP", final_df$IID), ]
nrow(df_subset) # 42942

df_ref <- final_df[!grepl("SP", final_df$IID), ]
df_ref <- df_ref %>% distinct(IID, .keep_all = TRUE)
nrow(df_ref)

df_long <- df_subset %>%
  pivot_longer(
    cols = c("AFR", "AMR", "EAS", "EUR", "SAS"),   
    names_to = "InferredPop",                     
    values_to = "Proportion"                     
  )

p <- ggplot(df_long, aes(x = IID, y = Proportion, fill = InferredPop)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ pop, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Samples", y = "Proportion", fill = "Inferred Population")

ggsave(paste0(Folder, "/OUTPUTS/",outputs_prefix, "_popinfer.png"), plot = p, dpi = 300, width = 12, height = 4)

df_long <- df_ref %>%
  pivot_longer(
    cols = c("AFR", "AMR", "EAS", "EUR", "SAS"),   
    names_to = "InferredPop",                     
    values_to = "Proportion"                     
  )

p <- ggplot(df_long, aes(x = IID, y = Proportion, fill = InferredPop)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ pop, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Samples", y = "Proportion", fill = "Inferred Population")

ggsave(paste0(Folder, "/OUTPUTS/",outputs_prefix, "_pop1000G.png"), plot = p, dpi = 300, width = 12, height = 4)

df_unique <- final_df %>% distinct(IID, .keep_all = TRUE)

df_long <- df_unique %>%
  pivot_longer(
    cols = c("AFR", "AMR", "EAS", "EUR", "SAS"),   
    names_to = "InferredPop",                     
    values_to = "Proportion"                     
  )

p <- ggplot(df_long, aes(x = IID, y = Proportion, fill = InferredPop)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ pop, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Samples", y = "Proportion", fill = "Inferred Population")

ggsave(paste0(Folder, "/OUTPUTS/",outputs_prefix, "_pop_us&1000G.png"), plot = p, dpi = 300, width = 12, height = 4)

pop_cols <- c("AFR", "AMR", "EAS", "EUR", "SAS")

df_subset <- df_subset %>%
  rowwise() %>%
  mutate(
    main_pop = c("AFR", "AMR", "EAS", "EUR", "SAS")[which.max(c(AFR, AMR, EAS, EUR, SAS))]
  ) %>%
  ungroup()

df_subset %>% filter(main_pop == "AFR") %>% select(FID, IID) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_AFR.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)
df_subset %>% filter(main_pop == "EUR") %>% select(FID, IID) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_EUR.txt"),quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)
df_subset %>% filter(main_pop == "AMR") %>% select(FID, IID) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_AMR.txt"),quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)
df_subset %>% filter(main_pop == "EAS") %>% select(FID, IID) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_EAS.txt"),quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)
df_subset %>% filter(main_pop == "SAS") %>% select(FID, IID) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_SAS.txt"),quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)
df_subset %>% select(FID, IID, main_pop) %>% write.table(paste0(Folder, "/DATA/", outputs_prefix,"_popinfer.txt"),quote = FALSE, sep = "\t", col.names = FALSE, row.names=FALSE)

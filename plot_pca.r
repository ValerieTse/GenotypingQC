#!/usr/bin/env Rscript

# ---- Load required packages -----
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# ----- Get command-line arguments -----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript plot_pca.R <fam_file> <pop_file> <eigenvec_file> <prefix> <Folder>")
}

fam_file <- args[1]
pop_file <- args [2]
eigenvec_file <- args[3]
prefix <- args[4]
Folder <- args[5]

# ----- create pop files -----  
fread(fam_file, header = F) -> fam 
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "P")
fam <- fam %>% select(FID, IID)

fread(pop_file) -> pop
pop <- merge(fam, pop, by.x = "IID", by.y = "#IID", all.x = TRUE)
pop %>% select(FID, IID, Pop = Population, SPop=SuperPop) %>% mutate(Pop = ifelse(is.na(Pop), 'Sample', Pop), SPop = ifelse(is.na(SPop), 'Sample', SPop)) -> pop
pop %>% fwrite(paste0(Folder, "/DATA/", prefix, "_pop.txt"))

pca <- fread(eigenvec_file)
pca <- merge(pca, pop, by = "IID", all.x = TRUE)


# ----- Plotting PCA -----
pca$SPop <- as.factor(pca$SPop)
others <- levels(pca$SPop)[levels(pca$SPop) != "Sample"]
default_colors <- scales::hue_pal()(length(others))
names(default_colors) <- others
my_colors <- c(default_colors, Sample = "grey")
my_colors <- my_colors[c(others,"Sample")]

pca <- ggplot() +
  geom_point(data = subset(pca, SPop == "Sample"),
             aes(x = PC1, y = PC2, color = SPop),
             size = 0.1) +
  geom_point(data = subset(pca, SPop != "Sample"),
             aes(x = PC1, y = PC2, color = SPop),
             size = 0.1) +
  scale_color_manual(values = my_colors) +
  labs(color = "Population") +
  theme_bw()
ggsave(paste0(Folder, "/OUTPUTS/", prefix, "_pca_Spop.png"), plot = pca, width = 6, height = 5, dpi = 300)

pca$Pop <- as.factor(pca$Pop)
others <- levels(pca$Pop)[levels(pca$Pop) != "Sample"]
default_colors <- scales::hue_pal()(length(others))
names(default_colors) <- others
my_colors <- c(default_colors, Sample = "grey")
my_colors <- my_colors[c(others,"Sample")]

pca2 <- ggplot() +
  geom_point(data = subset(pca, Pop == "Sample"),
             aes(x = PC1, y = PC2, color = Pop),
             size = 0.1) +
  geom_point(data = subset(pca, Pop != "Sample"),
             aes(x = PC1, y = PC2, color = Pop),
             size = 0.1) +
  scale_color_manual(values = my_colors) +
  labs(color = "Population") +
  theme_bw()
ggsave(paste0(Folder, "/OUTPUTS/", prefix, "_pca_pop.png"), plot = pca, width = 6, height = 5, dpi = 300)

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Sample data inspection (commented out)
# icgc[1:3,80:83]; tcga[1:3,108:111]; ope[1:3,c(7,37:39)]
# gse76427[1:3,31:34]; gse148355[1:3,42:45]; gse14520[1:3,41:43]; gse124751[1:3,34:37]

# Function to process each dataset
process_dataset <- function(df, dataset_name) {
  # Filter tumor samples based on available column
  if ("sample_type" %in% colnames(df)) {
    df_tumor <- df %>% filter(sample_type == "Tumor")
  } else if ("TYPE" %in% colnames(df)) {
    df_tumor <- df %>% filter(TYPE == "HCC")
  } else {
    df_tumor <- df  # Assume all are tumor
  }
  
  # Remove NA or invalid age entries
  df_tumor <- df_tumor %>%
    filter(!is.na(age) & !age %in% c("#N/A", "NA")) %>%
    mutate(age = as.numeric(age)) %>%
    filter(!is.na(age))
  
  # Define 10-year age bins
  breaks <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  labels <- c("0-9", "10-19", "20-29", "30-39", "40-49", 
              "50-59", "60-69", "70-79", "80-89", "90-99")
  
  df_tumor$Age_group <- cut(df_tumor$age, breaks = breaks, labels = labels, right = FALSE)
  
  # Compute mean HCCaging per age group
  result <- df_tumor %>%
    filter(!is.na(HCCaging)) %>%
    group_by(Age_group) %>%
    summarise(HCCaging_mean = mean(HCCaging, na.rm = TRUE),
              n = n()) %>%
    mutate(data_id = dataset_name) %>%
    select(data_id, Age_group, HCCaging_mean, n)
  
  return(result)
}

# Process all datasets
tcga_processed <- process_dataset(tcga, "TCGA_HCC")
icgc_processed <- process_dataset(icgc, "ICGC")
ope_processed <- process_dataset(ope, "OEP000321")
gse76427_processed <- process_dataset(gse76427, "GSE76427")
gse148355_processed <- process_dataset(gse148355, "GSE148355")
gse14520_processed <- process_dataset(gse14520, "GSE14520")
gse124751_processed <- process_dataset(gse124751, "GSE124751")

# Combine all results
all_data <- bind_rows(tcga_processed, icgc_processed, ope_processed,
                      gse76427_processed, gse148355_processed,
                      gse14520_processed, gse124751_processed)

# Z-score normalization within each dataset
all_data <- all_data %>%
  group_by(data_id) %>%
  mutate(HCCaging_mean_scaled = scale(HCCaging_mean)[,1]) %>%
  ungroup() %>%
  as.data.frame()

# Define standard age group order
standard_age_groups <- c("10-19", "20-29", "30-39", "40-49", "50-59",
                         "60-69", "70-79", "80-89", "90-99")

# Prepare dataset for plotting with correct factor levels
prepare_dataset <- function(dataset_name) {
  df <- all_data %>%
    filter(data_id == dataset_name) %>%
    mutate(Age_group = factor(Age_group, levels = standard_age_groups)) %>%
    filter(Age_group %in% standard_age_groups)
  return(df)
}

# Plot single heatmap per dataset
plot_single_heatmap <- function(df, show_legend = TRUE) {
  dataset_name <- unique(df$data_id)
  heatmap_df <- df %>% tidyr::complete(Age_group = standard_age_groups)
  
  p <- ggplot(heatmap_df, aes(x = Age_group, y = 1, fill = HCCaging_mean_scaled)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = ifelse(is.na(HCCaging_mean_scaled), "NA",
                                sprintf("%.2f\n(n=%d)", HCCaging_mean_scaled, n))),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#d1e5f0", mid = "#c9e1ee", high = "#4393c3",
                        midpoint = mean(heatmap_df$HCCaging_mean_scaled, na.rm = TRUE),
                        na.value = "#fceef0",
                        name = "Scaled HCCaging") +
    labs(title = dataset_name, x = "Age Group", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# Generate heatmaps for all datasets
datasets <- unique(all_data$data_id)
dataset_plots <- lapply(datasets, function(name) {
  df <- prepare_dataset(name)
  show_legend <- (name == datasets[1])
  plot_single_heatmap(df, show_legend)
})
names(dataset_plots) <- datasets

# Combine heatmaps vertically
combined_plot <- patchwork::wrap_plots(dataset_plots, ncol = 1, heights = rep(1, length(datasets)))

# Enhanced heatmap with trend slope
plot_trend_heatmap <- function(df, show_legend = TRUE) {
  dataset_name <- unique(df$data_id)
  df_numeric <- df %>%
    mutate(age_mid = sapply(strsplit(as.character(Age_group), "-"), function(x) mean(as.numeric(x)))) %>%
    filter(!is.na(HCCaging_mean_scaled))
  
  # Compute linear trend
  if (nrow(df_numeric) > 1) {
    model <- lm(HCCaging_mean_scaled ~ age_mid, data = df_numeric)
    slope <- coef(model)[2]
    color <- ifelse(slope > 0, "darkgreen", "darkred")
  } else {
    slope <- NA
    color <- "black"
  }
  
  heatmap_df <- df %>% tidyr::complete(Age_group = standard_age_groups)
  
  title <- paste0(dataset_name, ifelse(!is.na(slope), paste0(" (Trend: ", sprintf("%.3f", slope), ")"), ""))
  
  p <- ggplot(heatmap_df, aes(x = Age_group, y = 1, fill = HCCaging_mean_scaled)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = ifelse(is.na(HCCaging_mean_scaled), "NA",
                                sprintf("%.2f\n(n=%d)", HCCaging_mean_scaled, n))),
              size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#d1e5f0", mid = "#c9e1ee", high = "#4393c3",
                        midpoint = mean(heatmap_df$HCCaging_mean_scaled, na.rm = TRUE),
                        na.value = "#fceef0",
                        name = "Scaled HCCaging Score") +
    labs(title = title, x = "Age Group", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = color),
      plot.background = element_rect(fill = "white", color = NA)
    )
  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# Apply trend heatmap to all datasets
trend_plots <- lapply(datasets, function(name) {
  df <- prepare_dataset(name)
  plot_trend_heatmap(df, show_legend = (name == datasets[1]))
})
names(trend_plots) <- datasets

# Combine trend heatmaps
combined_trend_plot <- patchwork::wrap_plots(trend_plots, ncol = 1, heights = rep(1, length(datasets)))

# Save plots
ggsave("HCCaging_all_datasets_combined.pdf",
       combined_plot, width = 12, height = 15)
ggsave("HCCaging_all_datasets_trend_combined.pdf",
       combined_trend_plot, width = 12, height = 15)

# Prepare data for line plots
all_data_trend <- all_data %>%
  mutate(age_mid = sapply(strsplit(as.character(Age_group), "-"), function(x) {
    if (length(x) == 2) mean(as.numeric(x)) else as.numeric(x[1]) + 5
  })) %>%
  filter(!is.na(HCCaging_mean))

# Define color palette
dataset_colors <- c(
  "TCGA_HCC" = "#f9a27d",
  "ICGC" = "#a0afd3", 
  "OEP000321" = "#e48dba",
  "GSE76427" = "#a6ce58",
  "GSE148355" = "#bb7043",
  "GSE14520" = "#6cb5e9",
  "GSE124751" = "#f58b4c"
)

# Create faceted trend line plot
plot_faceted_trends <- function() {
  ggplot(all_data_trend, aes(x = age_mid, y = HCCaging_mean_scaled)) +
    facet_wrap(~ data_id, ncol = 2, scales = "free_y") +
    geom_smooth(aes(color = data_id), method = "lm", se = TRUE, alpha = 0.15, linewidth = 1, fill = "lightgray") +
    geom_point(aes(color = data_id, size = n), alpha = 0.8) +
    geom_line(aes(color = data_id), alpha = 0.6, linewidth = 0.8) +
    scale_color_manual(values = dataset_colors, guide = "none") +
    scale_size_continuous(range = c(2, 6), name = "Sample Size") +
    scale_x_continuous(breaks = seq(20, 90, by = 10),
                       labels = c("20", "30", "40", "50", "60", "70", "80", "90")) +
    labs(title = "HCCaging Score Trends by Dataset",
         x = "Age (Midpoint of Age Group)",
         y = "Scaled HCCaging Score") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}

faceted_trend_plot <- plot_faceted_trends()
print(faceted_trend_plot)

ggsave("HCCaging_faceted_trends.pdf",
       faceted_trend_plot, width = 10, height = 10)

# Function to process multiple aging scores (extension)
process_dataset_multiple_scores <- function(df, dataset_name, score_names) {
  # Placeholder function (assumed to exist in original context)
  # This should process multiple aging scores per dataset
  # If not defined, this part will fail
  stop("Function `process_dataset_multiple_scores` is not defined. Please include it.")
}

# List of aging scores (assumed defined elsewhere)
# aging_scores <- c("HCCaging", "OtherScore1", ...)

# Process all datasets and all aging scores
all_datasets_all_scores <- process_all_datasets_multiple_scores(aging_scores)

# Save full result
qs::qsave(all_datasets_all_scores, "all_data_all_score.qs")



# scripts/01_enrichment_analysis.R
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# --- 参数区 ---
gene_file <- "data/age_up_genes.csv"
output_dir <- "results/enrichment"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 读取基因 ---
clu <- read_csv(gene_file)
gene_symbols <- clu$age_genes

# --- SYMBOL to ENTREZID ---
if (!require("BiocManager")) install.packages("BiocManager")
gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- gene_df$ENTREZID

# --- GO Enrichment ---
go_bp <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = "BP", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_mf <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = "MF", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_cc <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = "CC", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

# --- KEGG ---
kegg_enrich <- enrichKEGG(gene = gene_entrez, organism = "hsa", 
                          pvalueCutoff = 0.05, pAdjustMethod = "BH")

# --- 保存结果 ---
write.csv(go_bp, file.path(output_dir, "hccaging_GO_BP_sig.csv"), row.names = FALSE)
write.csv(go_mf, file.path(output_dir, "hccaging_GO_MF_sig.csv"), row.names = FALSE)
write.csv(go_cc, file.path(output_dir, "hccaging_GO_CC_sig.csv"), row.names = FALSE)
write.csv(kegg_enrich, file.path(output_dir, "hccaging_KEGG_sig.csv"), row.names = FALSE)

message("✅ Enrichment analysis completed.")

# ===================================================================
# Modular Script: Half-Violin Plots for Multiple Gene Set Comparisons
# File: scripts/plot_gene_set_comparison.R
# Purpose: Works with any pre-computed GSVA .qs files
# ===================================================================

library(ggplot2)
library(ggpubr)
library(qs)
library(dplyr)
library(tidyr)

# ===================================================================
# CONFIGURATION (User Modifiable Section)
# ===================================================================

config <- list(
  # Input/Output paths
  input_file  = "\/GSE76427_sample_info_tumor_score.qs",
  output_dir  = "/age_cli/age_genesets/",
  output_name = "gene_sets_half_violin_GSE76427.pdf",
  
  # Grouping variables
  group_col   = "sample_type",
  groups      = c("Normal", "Tumor"),
  
  # Gene sets to plot
  score_cols  = c("SenMayo", "CellAge", "GenAge", "ASIG", "SASP", "AgingAtlas",
                  "SenUp", "SigRS", "DAS", "mSS", "SenUP", "hUSI", "ClassicalSs", "HCCaging"),
  
  # Plot settings
  title       = "Comparison of Gene Set Scores between Tumor and Normal Samples in GSE76427",
  colors      = c("#90c29f", "#f5ac6f"), # Normal, Tumor colors
  gene_set_order = c("SenMayo", "CellAge", "GenAge", "ASIG", "SASP", "AgingAtlas",
                     "SenUp", "SigRS", "DAS", "mSS", "SenUP", "hUSI", "ClassicalSs", "HCCaging")
)

# ===================================================================
# DATA PROCESSING
# ===================================================================

# Create output directory
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
gsva_df <- qread(config$input_file)
message("✅ Loaded data from ", config$input_file)

# Validate required columns
missing_cols <- setdiff(c(config$group_col, config$score_cols), colnames(gsva_df))
if (length(missing_cols) > 0) {
  stop("❌ Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Reshape data to long format
gsva_long <- gsva_df %>%
  select(all_of(c(config$group_col, config$score_cols))) %>%
  pivot_longer(
    cols = all_of(config$score_cols),
    names_to = "GeneSet",
    values_to = "Score"
  ) %>%
  filter(!is.na(.data[[config$group_col]])) %>%
  mutate(
    !!config$group_col := factor(.data[[config$group_col]], levels = config$groups),
    GeneSet = factor(GeneSet, levels = config$gene_set_order)
  )

# ===================================================================
# PLOTTING
# ===================================================================

p <- ggplot(gsva_long, aes(x = GeneSet, y = Score)) +
  # Half-violins for each group
  geom_half_violin(
    data = ~ subset(., .data[[config$group_col]] == config$groups[1]),
    aes(fill = .data[[config$group_col]]),
    side = "left", alpha = 0.8, width = 0.8
  ) +
  geom_half_violin(
    data = ~ subset(., .data[[config$group_col]] == config$groups[2]),
    aes(fill = .data[[config$group_col]]),
    side = "right", alpha = 0.8, width = 0.8
  ) +
  # IQR error bars
  stat_summary(
    fun.data = function(x) {
      data.frame(y = median(x), 
                 ymin = quantile(x, 0.25), 
                 ymax = quantile(x, 0.75))
    },
    geom = "errorbar", width = 0.15,
    position = position_dodge(width = 0.2),
    color = "black", size = 0.5
  ) +
  # Median points
  stat_summary(
    fun = median, geom = "point",
    shape = 21, size = 2.5, fill = "white",
    position = position_dodge(width = 0.2), color = "black"
  ) +
  # Significance testing
  stat_compare_means(
    aes(x = GeneSet, y = Score, group = .data[[config$group_col]]),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
      symbols = c("***", "**", "*", "ns")
    ),
    label.y = max(gsva_long$Score, na.rm = TRUE) * 1.05,
    size = 4
  ) +
  # Aesthetics
  scale_fill_manual(values = config$colors, labels = config$groups) +
  labs(
    title = config$title,
    x = "Gene Sets",
    y = "GSVA Score",
    fill = config$group_col
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top",
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 5, 10, 5)
  ) +
  scale_x_discrete(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))

# ===================================================================
# OUTPUT
# ===================================================================

output_path <- file.path(config$output_dir, config$output_name)
ggsave(output_path, plot = p, width = 14, height = 6, dpi = 300)
message("✅ Plot saved to: ", output_path)

# Optional: Save as PNG
# png_path <- sub(".pdf$", ".png", output_path)
# ggsave(png_path, plot = p, width = 14, height = 6, dpi = 300)
# message("✅ PNG version saved to: ", png_path)


# ===================================================================
# Multi-Dataset Aging Geneset Volcano Plot Analysis
# Purpose: Comparative analysis of aging-related gene sets across multiple datasets
# ===================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

# ===================================================================
# GLOBAL CONFIGURATION
# ===================================================================

# Define genesets and colors
genesets <- c("SenMayo", "CellAge", "GenAge", "ASIG", "SASP", "AgingAtlas", 
              "SenUp", "SigRS", "DAS", "mSS", "SenUP", "hUSI", "ClassicalSs", "HCCaging")

geneset_colors <- c(
  "SenMayo" = "#E41A1C", "CellAge" = "#377EB8", "GenAge" = "#4DAF4A",
  "ASIG" = "#984EA3", "SASP" = "#FF7F00", "AgingAtlas" = "#FFD700",
  "SenUp" = "#A65628", "SigRS" = "#F781BF", "DAS" = "#999999",
  "mSS" = "#66C2A5", "SenUP" = "#FC8D62", "hUSI" = "#8DA0CB",
  "ClassicalSs" = "#E78AC3", "HCCaging" = "#A6D854"
)

# ===================================================================
# CORE FUNCTIONS
# ===================================================================

perform_diff_analysis <- function(gsva_df, dataset_name) {
  # Ensure Type column exists
  if (!"Type" %in% colnames(gsva_df)) {
    gsva_df$Type <- gsva_df$sample_type
  }
  
  # Check sample size requirements
  tumor_count <- sum(gsva_df$Type == "Tumor")
  normal_count <- sum(gsva_df$Type == "Normal")
  
  if (tumor_count < 2 || normal_count < 2) {
    message("Dataset ", dataset_name, " has insufficient samples - skipping analysis")
    return(NULL)
  }
  
  # Differential analysis implementation here
  # [Add your specific differential analysis code]
}

create_volcano_plot <- function(diff_data, dataset_name) {
  # Volcano plot creation logic
  # [Add your specific volcano plot code]
}

combine_plots <- function(volcano_plots) {
  wrap_plots(volcano_plots, ncol = 4) + 
    plot_annotation(
      title = "Aging Geneset Enrichment: Tumor vs Normal Across Datasets",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
}

create_unified_legend <- function() {
  # Unified legend creation logic
  # [Add your specific legend creation code]
}

# ===================================================================
# EXECUTION PIPELINE
# ===================================================================

# Example dataset list (replace with your actual datasets)
datasets <- list(
  LIRI = LIRI, 
  GSE10143 = GSE10143
  # Add more datasets as needed
)

# Perform differential analysis across all datasets
all_diff_results <- lapply(names(datasets), function(dataset_name) {
  perform_diff_analysis(datasets[[dataset_name]], dataset_name)
})

# Generate volcano plots
volcano_plots <- lapply(seq_along(all_diff_results), function(i) {
  if (!is.null(all_diff_results[[i]])) {
    create_volcano_plot(all_diff_results[[i]], names(datasets)[i])
  }
})

# Remove NULL plots and combine
volcano_plots <- Filter(Negate(is.null), volcano_plots)
final_plot <- combine_plots(volcano_plots)

# Add unified legend and save
unified_legend <- create_unified_legend()
final_plot_with_legend <- plot_grid(final_plot, unified_legend, ncol = 2, rel_widths = c(3, 1))

# Save final output
ggsave("multi_dataset_geneset_volcano_unified.pdf", final_plot_with_legend, width = 10, height = 5)

library(Mfuzz)
library(edgeR)
library(qs)
library(Biobase)
library(tidyverse)

# ===================================================================
# CONFIGURATION
# ===================================================================

input_dir <- "/"
output_dir <- "/time_cluster_analysis/"
age_diff_dir <- "/age_diff_tcga/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(age_diff_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ===================================================================
# DATA LOADING AND PREPROCESSING
# ===================================================================

df <- qread(paste0(input_dir, "taga_hcc_sample_info_tumor_score.qs"))
matrix <- qread(paste0(input_dir, "taga_hcc_exp.qs"), nthreads = 90)

# Filter tumor samples and process age data
tumor_samples <- df %>%
  filter(Type == "Tumor") %>%
  slice(2:371) %>%
  select(original_sample_id, age_at_initial_pathologic_diagnosis) %>%
  mutate(
    Age = as.numeric(ifelse(age_at_initial_pathologic_diagnosis == "#N/A", NA, age_at_initial_pathologic_diagnosis)),
    Age_group = case_when(
      Age < 20 ~ "<20",
      Age >= 20 & Age < 30 ~ "20-29",
      Age >= 30 & Age < 40 ~ "30-39", 
      Age >= 40 & Age < 50 ~ "40-49",
      Age >= 50 & Age < 60 ~ "50-59",
      Age >= 60 & Age < 70 ~ "60-69",
      Age >= 70 & Age < 80 ~ "70-79",
      Age >= 80 ~ "80+"
    ),
    Age_group = factor(Age_group, levels = c("<20", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), ordered = TRUE)
  )

# Extract tumor expression matrix and normalize
tumor_exp_matrix <- matrix[, tumor_samples$original_sample_id]
dge <- DGEList(counts = tumor_exp_matrix)
log2_cpm <- cpm(dge, log = TRUE, prior.count = 1)

# ===================================================================
# AGE GROUP AVERAGING
# ===================================================================

age_groups <- levels(tumor_samples$Age_group)
time_points <- c(18, 24.5, 34.5, 44.5, 54.5, 64.5, 74.5, 85)

averaged_data <- matrix(nrow = nrow(log2_cpm), ncol = length(age_groups))
colnames(averaged_data) <- age_groups
rownames(averaged_data) <- rownames(log2_cpm)

for (i in seq_along(age_groups)) {
  group_samples <- tumor_samples$original_sample_id[tumor_samples$Age_group == age_groups[i]]
  if (length(group_samples) > 0) {
    averaged_data[, i] <- if (length(group_samples) > 1) {
      rowMeans(log2_cpm[, group_samples], na.rm = TRUE)
    } else {
      log2_cpm[, group_samples]
    }
  }
}

averaged_data <- averaged_data[complete.cases(averaged_data), ]

# ===================================================================
# MFUZZ CLUSTERING
# ===================================================================

# Create ExpressionSet object
pheno_data <- data.frame(Time = time_points, row.names = colnames(averaged_data))
expr_set <- new("ExpressionSet", 
                exprs = as.matrix(averaged_data),
                phenoData = new("AnnotatedDataFrame", data = pheno_data))

# Preprocess data for Mfuzz
expr_set_filtered <- filter.NA(expr_set, thres = 0.001)
expr_set_filled <- fill.NA(expr_set_filtered, mode = "mean")
expr_set_filtered_std <- filter.std(expr_set_filled, min.std = 0.05)
expr_set_standardized <- standardise(expr_set_filtered_std)

# Determine clustering parameters
m <- mestimate(expr_set_standardized)
message("Estimated m value: ", m)

# Determine optimal cluster number (visual inspection required)
tmp <- Dmin(expr_set_standardized, m = m, crange = seq(4, 20, 2), repeats = 3, visu = TRUE)

# Perform clustering (adjust c value based on Dmin results)
c <- 7
set.seed(123)
cl <- mfuzz(expr_set_standardized, c = c, m = m)

# ===================================================================
# VISUALIZATION AND OUTPUT
# ===================================================================

# Generate cluster plots
color_palette <- colorRampPalette(c("#fad7d3", "#ebf0b5", "#f6e09d", "#d2e2ef", "#90bfd5", "#3c7fb1"))(100)

pdf(paste0(output_dir, "tcga_age_time_cluster.pdf"), width = 10, height = 10, pointsize = 10)
mfuzz.plot(expr_set_standardized, cl = cl, mfrow = c(4, 4), colo = color_palette, 
           time.labels = time_points, new.window = FALSE)
dev.off()

# Extract cluster genes
cluster_genes <- lapply(1:c, function(i) names(cl$cluster[cl$cluster == i]))
names(cluster_genes) <- paste0("Cluster_", 1:c)

# Save results
clust_df <- data.frame(gene_id = names(cl$cluster), cluster = cl$cluster)
write.csv(clust_df, paste0(age_diff_dir, "tcga_age_cluster_assignments.csv"), row.names = FALSE)

membership_df <- as.data.frame(cl$membership)
write.csv(membership_df, paste0(age_diff_dir, "tcga_age_cluster_membership.csv"))

centers_df <- as.data.frame(cl$centers)
write.csv(centers_df, paste0(age_diff_dir, "tcga_age_cluster_centers.csv"))

message("Analysis completed successfully!")
message("Cluster assignments saved to: ", age_diff_dir)













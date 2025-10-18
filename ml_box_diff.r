rm(list=ls())
library(qs)
library(Matrix)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(homologene)
library(ConsensusClusterPlus)
library(limma)
library(future.apply)
options(future.globals.maxSize = 1.6 * 1024^4)


conda activate seurat4_plot
####################################################------------------------------------------------Bar plot
# Prepare data
library(tidyr)
library(dplyr)
library(ggplot2)


df <- read.csv("/Tumor_Normal_gbm_8samples.csv")
rownames(df) <- df$X
df$X <- NULL

dataset_colors <- c(
  "TCGA" = "#f9a27d",
  "LIRI" = "#a0afd3", 
  "OEP000321" = "#e48dba",
  "GSE76427" = "#a6ce58",
  "GSE10143" = "#bb7043",
  "GSE14520" = "#6cb5e9", 
  "GSE25097" = "#f58b4c",
  "GSE63898" = "#f9c27d" 
)
# Custom gene set colors (cyclic usage)
custom_colors <- c(
  "#a8b8d9",  # Desaturated cool blue
  "#7fd1e6",
  "#c6e4f0",
  "#93c2a9",
  "#c0dcb3",
  "#e0edd6",
  "#fef2c0",
  "#ffd9cc",
  "#f7c0a8",
  "#f7a592"   # Desaturated warm red
)

# Add gene set names as column
df_long <- df %>%
  rownames_to_column("GeneSet") %>%
  pivot_longer(cols = -GeneSet, names_to = "Dataset", values_to = "Value")

# Calculate statistics for each gene set
df_summary <- df_long %>%
  group_by(GeneSet) %>%
  summarise(
    MeanValue = mean(Value, na.rm = TRUE),
    SDValue = sd(Value, na.rm = TRUE),
    SEValue = SDValue / sqrt(n())
  )

# Merge back to original data
df_plot <- df_long %>%
  left_join(df_summary, by = "GeneSet")

# Set gene set order (by mean from high to low)
gene_set_order <- df_summary %>%
  arrange(desc(MeanValue)) %>%
  pull(GeneSet)

df_plot$GeneSet <- factor(df_plot$GeneSet, levels = gene_set_order)

# Assign colors to each gene set
gene_set_colors <- setNames(
  custom_colors[1:length(gene_set_order) %% length(custom_colors) + 1],
  gene_set_order
)

# Create plot
p <- ggplot(df_plot, aes(x = GeneSet)) +
  # Bar plot (mean) - different colors for each gene set
  geom_bar(aes(y = MeanValue, fill = GeneSet), 
           stat = "identity", 
           alpha = 0.8, 
           width = 0.8) +  # Set width to 0.8 for tighter bars
  
  # Error bars (standard deviation)
  geom_errorbar(aes(ymin = MeanValue - SDValue, ymax = MeanValue + SDValue),
                width = 0.3,  # Error bar width
                size = 0.8,   # Error bar thickness
                color = "black",
                alpha = 0.7) +
  
  # Points for four datasets
  geom_point(aes(y = Value, color = Dataset), 
             position = position_dodge(width = 0.5),  # Slightly spread out points
             size = 2.5, 
             alpha = 0.9) +
  
  # Set colors
  scale_fill_manual(values = gene_set_colors) +
  scale_color_manual(values = dataset_colors) +
  # Axes and titles
  labs(title = "Gene Set Performance Across Datasets",
       x = "Gene Set",
       y = "AUC",
       fill = "Gene Set",
       color = "Dataset") +
  
  # Theme settings
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
     legend.position = "right",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Ensure background is pure white
    plot.background = element_blank(),   # Ensure plot area background is pure white
    axis.line = element_line(color = "black", size = 0.5)  # Add XY axis lines
  ) +
  
  # Y-axis range adjustment to leave space for labels and error bars
  ylim(0, max(df_plot$Value, na.rm = TRUE) * 1.0)

# Display plot
print(p)
ggsave("/Tumor_Normal_gbm_new.pdf", p, width = 10, height = 7)



################################---------------------------------------
df <- read.csv("/Tumor_Normal_rf_8samples.csv")
rownames(df) <- df$X
df$X <- NULL

dataset_colors <- c(
  "TCGA" = "#f9a27d",
  "LIRI" = "#a0afd3", 
  "OEP000321" = "#e48dba",
  "GSE76427" = "#a6ce58",
  "GSE10143" = "#bb7043",
  "GSE14520" = "#6cb5e9", 
  "GSE25097" = "#f58b4c",
  "GSE63898" = "#f9c27d" 
)
# Custom gene set colors (cyclic usage)
custom_colors <- c(
  "#a8b8d9",  # Desaturated cool blue
  "#7fd1e6",
  "#c6e4f0",
  "#93c2a9",
  "#c0dcb3",
  "#e0edd6",
  "#fef2c0",
  "#ffd9cc",
  "#f7c0a8",
  "#f7a592"   # Desaturated warm red
)

# Add gene set names as column
df_long <- df %>%
  rownames_to_column("GeneSet") %>%
  pivot_longer(cols = -GeneSet, names_to = "Dataset", values_to = "Value")

# Calculate statistics for each gene set
df_summary <- df_long %>%
  group_by(GeneSet) %>%
  summarise(
    MeanValue = mean(Value, na.rm = TRUE),
    SDValue = sd(Value, na.rm = TRUE),
    SEValue = SDValue / sqrt(n())
  )

# Merge back to original data
df_plot <- df_long %>%
  left_join(df_summary, by = "GeneSet")

# Set gene set order (by mean from high to low)
gene_set_order <- df_summary %>%
  arrange(desc(MeanValue)) %>%
  pull(GeneSet)

df_plot$GeneSet <- factor(df_plot$GeneSet, levels = gene_set_order)

# Assign colors to each gene set
gene_set_colors <- setNames(
  custom_colors[1:length(gene_set_order) %% length(custom_colors) + 1],
  gene_set_order
)

# Create plot
p <- ggplot(df_plot, aes(x = GeneSet)) +
  # Bar plot (mean) - different colors for each gene set
  geom_bar(aes(y = MeanValue, fill = GeneSet), 
           stat = "identity", 
           alpha = 0.8, 
           width = 0.8) +  # Set width to 0.8 for tighter bars
  
  # Error bars (standard deviation)
  geom_errorbar(aes(ymin = MeanValue - SDValue, ymax = MeanValue + SDValue),
                width = 0.3,  # Error bar width
                size = 0.8,   # Error bar thickness
                color = "black",
                alpha = 0.7) +
  
  # Points for four datasets
  geom_point(aes(y = Value, color = Dataset), 
             position = position_dodge(width = 0.5),  # Slightly spread out points
             size = 2.5, 
             alpha = 0.9) +
  
  # Set colors
  scale_fill_manual(values = gene_set_colors) +
  scale_color_manual(values = dataset_colors) +
  # Axes and titles
  labs(title = "Gene Set Performance Across Datasets",
       x = "Gene Set",
       y = "AUC",
       fill = "Gene Set",
       color = "Dataset") +
  
  # Theme settings
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
     legend.position = "right",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Ensure background is pure white
    plot.background = element_blank(),   # Ensure plot area background is pure white
    axis.line = element_line(color = "black", size = 0.5)  # Add XY axis lines
  ) +
  
  # Y-axis range adjustment to leave space for labels and error bars
  ylim(0, max(df_plot$Value, na.rm = TRUE) * 0.99)

# Display plot
print(p)
ggsave("/Tumor_Normal_randomforest_new.pdf", p, width = 10, height = 7)
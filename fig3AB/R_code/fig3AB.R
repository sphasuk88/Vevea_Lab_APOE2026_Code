# fig3AB_heatmap_and_clusterplots.R
# Purpose:
#   Generate the Fig3A/3B heatmap and downstream cluster summary plots
#   from a proteomics CSV file.
#
# Expected folder layout:
#   fig3AB/
#   ├── input/
#   │   └── demo_fig3A.csv
#   ├── output/
#   └── R_code/
#
# Reviewer instructions:
#   1. Run this script in R or RStudio.
#   2. When prompted, manually select the input CSV file from the input folder.
#   3. The script will automatically determine the project folder and save all
#      output files to the output folder.

# -------------------------------
# Install/load required packages
# -------------------------------
required_packages <- c("pheatmap", "ggplot2", "dplyr", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)

# -------------------------------
# Select input file manually
# -------------------------------
input_file <- file.choose()
input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)

# -------------------------------
# Derive project and output paths
# -------------------------------
project_dir <- dirname(dirname(input_file))
output_dir <- file.path(project_dir, "output")

# Create output folder if missing
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------
# Print resolved paths
# -------------------------------
cat("Input file: ", input_file, "\n", sep = "")
cat("Project directory: ", project_dir, "\n", sep = "")
cat("Output directory: ", output_dir, "\n\n", sep = "")

# -------------------------------
# Load CSV file
# -------------------------------
all_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

cat("Input data dimensions: ", nrow(all_data), " rows x ", ncol(all_data), " columns\n\n", sep = "")

# -------------------------------
# Check required columns
# -------------------------------
if (!("FDR" %in% colnames(all_data))) {
  stop("Required column 'FDR' not found in the input file.")
}
if (!("GN" %in% colnames(all_data))) {
  stop("Required column 'GN' not found in the input file.")
}

# -------------------------------
# Filter the data for adjusted p-value less than 0.05
# -------------------------------
all_data_adjPvalue <- all_data[all_data$FDR < 0.05, ]

cat("Rows with FDR < 0.05: ", nrow(all_data_adjPvalue), "\n\n", sep = "")

# -------------------------------
# Build heatmap matrix
# Assuming columns 4:51 are the intended numeric expression columns
# -------------------------------
data_for_heatmap <- as.matrix(all_data_adjPvalue[, 4:51])
mode(data_for_heatmap) <- "numeric"

rownames(data_for_heatmap) <- all_data_adjPvalue$GN

# -------------------------------
# Extract sex information from column names
# -------------------------------
sex_info <- gsub(".*_([MF])_.*", "\\1", colnames(data_for_heatmap))

col_annotation <- data.frame(
  Sex = factor(sex_info, levels = c("M", "F")),
  stringsAsFactors = FALSE
)
rownames(col_annotation) <- colnames(data_for_heatmap)

# Reorder the matrix based on sex information
order_indices <- order(col_annotation$Sex)
data_for_heatmap <- data_for_heatmap[, order_indices, drop = FALSE]
col_annotation <- col_annotation[order_indices, , drop = FALSE]

# Define a color palette for the sex annotation
sex_colors <- c(M = "darkgreen", F = "orange")

# -------------------------------
# Generate heatmap
# -------------------------------
heatmap_result <- pheatmap(
  data_for_heatmap,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  cutree_rows = 4,
  color = colorRampPalette(c("blue", "white", "red"))(255),
  main = "Clustered Heatmap of Gene Expression by Protein",
  fontsize_row = 10,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = col_annotation,
  annotation_colors = list(Sex = sex_colors)
)

# Save heatmap
heatmap_pdf <- file.path(output_dir, "heatmap.pdf")
pdf(heatmap_pdf, width = 13, height = 7)
print(heatmap_result)
dev.off()

cat("Saved heatmap PDF: ", normalizePath(heatmap_pdf, winslash = "/"), "\n", sep = "")

# -------------------------------
# Extract clusters
# -------------------------------
clusters <- cutree(heatmap_result$tree_row, k = 4)

cluster_map <- data.frame(
  GN = rownames(data_for_heatmap),
  Cluster = clusters,
  row.names = NULL
)

# Join to full table
all_data <- all_data %>% left_join(cluster_map, by = "GN")

# -------------------------------
# Create DEPs-only file with cluster labels
# -------------------------------
deps30_with_clusters <- all_data_adjPvalue %>%
  left_join(cluster_map, by = "GN") %>%
  arrange(Cluster, GN)

deps_file <- file.path(output_dir, "clusters.csv")
write.csv(deps30_with_clusters, deps_file, row.names = FALSE)

cat("Saved DEP cluster table: ", normalizePath(deps_file, winslash = "/"), "\n", sep = "")

# -------------------------------
# Identify numeric sample columns for scaling
# -------------------------------
numeric_columns <- grep("b[1-3]_s[0-9]+_[FM]_[0-9]+", names(deps30_with_clusters), value = TRUE)

if (length(numeric_columns) == 0) {
  stop("No numeric sample columns matching 'b[1-3]_s[0-9]+_[FM]_[0-9]+' were found.")
}

# -------------------------------
# Perform row-wise z-score normalization
# -------------------------------
deps30_with_clusters_scaled <- deps30_with_clusters %>%
  select(all_of(numeric_columns)) %>%
  t() %>%
  scale(center = TRUE, scale = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  `colnames<-`(numeric_columns) %>%
  bind_cols(deps30_with_clusters %>% select(-all_of(numeric_columns)), .)

# -------------------------------
# Reshape data from wide to long format
# -------------------------------
long_data <- pivot_longer(
  deps30_with_clusters_scaled,
  cols = numeric_columns,
  names_to = "sample_label",
  values_to = "expression"
)

# Extract age and gender from sample labels
long_data <- long_data %>%
  mutate(
    Age = as.numeric(sub(".*_([0-9]+)$", "\\1", sample_label)),
    Gender = sub(".*_([MF])_.*", "\\1", sample_label)
  )

# -------------------------------
# Calculate median and SEM by cluster, age, and gender
# -------------------------------
results_scaled <- long_data %>%
  group_by(Cluster, Gender, Age) %>%
  summarise(
    Median = median(expression, na.rm = TRUE),
    SEM = sd(expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

results_file <- file.path(output_dir, "cluster_gender_age_medianSEM.csv")
write.csv(results_scaled, results_file, row.names = FALSE)

cat("Saved median/SEM summary: ", normalizePath(results_file, winslash = "/"), "\n", sep = "")

# -------------------------------
# Plot cluster summaries
# -------------------------------
num_clusters <- 3
gender_colors <- c("M" = "darkgreen", "F" = "orange")
age_breaks <- sort(unique(results_scaled$Age))

for (cluster_id in 1:num_clusters) {
  cluster_data <- results_scaled %>%
    dplyr::filter(Cluster == cluster_id) %>%
    dplyr::mutate(Gender = factor(Gender, levels = c("M", "F")))
  
  cluster_plot <- ggplot(cluster_data, aes(x = Age, y = Median, color = Gender, group = Gender)) +
    geom_ribbon(
      aes(ymin = Median - SEM, ymax = Median + SEM, fill = Gender),
      alpha = 0.4, color = NA
    ) +
    geom_line() +
    geom_point(size = 1) +
    scale_x_continuous(breaks = age_breaks, labels = age_breaks) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_color_manual(values = gender_colors) +
    scale_fill_manual(values = gender_colors) +
    labs(
      title = paste("Cluster", cluster_id),
      x = "Age (months)",
      y = "Normalized Median Expression"
    ) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none")
  
  cluster_pdf <- file.path(output_dir, sprintf("Cluster_%s_medianSEM_plot.pdf", cluster_id))
  ggsave(cluster_pdf, plot = cluster_plot, width = 3, height = 3, units = "in")
  
  cat("Saved cluster plot: ", normalizePath(cluster_pdf, winslash = "/"), "\n", sep = "")
}

cat("\nfig3AB workflow completed successfully.\n")
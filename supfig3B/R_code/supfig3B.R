# supfig3B_heatmap_clean.R
# Purpose:
#   Generate a Pearson sample-correlation heatmap from a proteomics CSV file.
#
# Reviewer instructions:
#   1. Run this script in R or RStudio.
#   2. When prompted, manually select the input CSV file.
#   3. The script will automatically determine the project folder and save output
#      files to the output folder.
#
# Assumptions:
#   - The first 3 columns are metadata columns.
#   - Sample columns follow the format:
#       b1_s13_M_3
#       b2_s18_IR_IR
#       b3_s20_Ctr_Ctr
#   - Extra non-sample columns such as representative_isoform will be ignored.

# -------------------------------
# Install/load required packages
# -------------------------------
required_packages <- c("tidyverse", "ggplot2", "ComplexHeatmap", "circlize")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}

library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

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

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------------
# Print resolved paths
# -------------------------------
cat("Input file: ", input_file, "\n", sep = "")
cat("Project directory: ", project_dir, "\n", sep = "")
cat("Output directory: ", output_dir, "\n\n", sep = "")

# -------------------------------
# Load input data
# -------------------------------
data <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input data dimensions: ", nrow(data), " rows x ", ncol(data), " columns\n\n", sep = "")

# -------------------------------
# 1) Expression matrix (sample columns only)
# -------------------------------
# Keep only columns that match expected sample naming patterns
sample_pattern <- "^b[123]_s[0-9]+_(M|F|Ctr|IR)_(3|6|12|18|24|30|Ctr|IR)$"
sample_cols <- grep(sample_pattern, colnames(data), value = TRUE)

if (length(sample_cols) == 0) {
  stop(
    paste0(
      "No sample columns detected.\n",
      "Expected columns like: b1_s13_M_3, b2_s18_IR_IR, b3_s20_Ctr_Ctr"
    )
  )
}

cat("Number of detected sample columns: ", length(sample_cols), "\n", sep = "")
cat("Sample columns used for heatmap:\n")
cat(paste0("  - ", sample_cols), sep = "\n")
cat("\n\n")

expr_df <- data[, sample_cols, drop = FALSE]
expr_df[] <- lapply(expr_df, as.numeric)
expr <- as.matrix(expr_df)

# Set row names using accession column if available
if ("ProteinAccession" %in% colnames(data)) {
  rownames(expr) <- make.unique(as.character(data$ProteinAccession))
} else if ("Protein.Accession.." %in% colnames(data)) {
  rownames(expr) <- make.unique(as.character(data$Protein.Accession..))
} else {
  rownames(expr) <- paste0("row_", seq_len(nrow(expr)))
}

cat("Expression matrix dimensions: ", nrow(expr), " proteins x ", ncol(expr), " samples\n\n", sep = "")

# -------------------------------
# 2) Sample metadata from column names
# -------------------------------
parts <- strsplit(colnames(expr), "_")

sample_info <- tibble(
  SampleID = colnames(expr),
  Batch = sapply(parts, `[`, 1),
  Sample = sapply(parts, `[`, 2),
  Group  = sapply(parts, `[`, 3),
  AgeChr = sapply(parts, `[`, 4)
) %>%
  mutate(
    Group  = factor(Group, levels = c("M", "F", "Ctr", "IR")),
    AgeChr = factor(AgeChr, levels = c("3", "6", "12", "18", "24", "30", "Ctr", "IR"), ordered = TRUE),
    AgeNum = suppressWarnings(as.numeric(as.character(AgeChr))),
    Batch  = factor(Batch, levels = c("b1", "b2", "b3"))
  )

# Order columns by Group -> Age -> Batch -> Sample
ord <- with(sample_info, order(Group, AgeNum, Batch, Sample, na.last = TRUE))
sample_info <- sample_info[ord, , drop = FALSE]
expr <- expr[, sample_info$SampleID, drop = FALSE]

cat("Sample ordering complete.\n\n")

# -------------------------------
# 3) Pearson sample-correlation heatmap (NO clustering)
# -------------------------------
cor_mat <- cor(expr, method = "pearson", use = "pairwise.complete.obs")
cor_mat <- cor_mat[sample_info$SampleID, sample_info$SampleID]

# Annotation palettes
group_cols <- c(
  M   = "darkgreen",
  F   = "darkorange",
  Ctr = "#6D6875",
  IR  = "red"
)

age_cols <- c(
  "3"   = "#FDE725",
  "6"   = "#B5DE2B",
  "12"  = "#6CCE59",
  "18"  = "#35B779",
  "24"  = "#1F9E89",
  "30"  = "#26828E",
  "Ctr" = "#C0C0C0",
  "IR"  = "#8C8C8C"
)

batch_cols <- c(
  b1 = "#5A4635",
  b2 = "#E76F51",
  b3 = "#E9C46A"
)

ha <- HeatmapAnnotation(
  Group = sample_info$Group,
  Age   = sample_info$AgeChr,
  Batch = sample_info$Batch,
  col = list(
    Group = group_cols,
    Age   = age_cols,
    Batch = batch_cols
  )
)

# Heatmap color scale
lo <- 0.6
mid <- 0.8
hi <- 1.00
col_fun <- colorRamp2(c(lo, mid, hi), c("#3B4CC0", "white", "#B40426"))

ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = sample_info$SampleID,
  column_order = sample_info$SampleID,
  row_split = sample_info$Group,
  column_split = sample_info$Group,
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(at = c(lo, mid, hi)),
  rect_gp = gpar(col = "grey85", lwd = 0.6)
)

# -------------------------------
# 4) Save heatmap to PDF
# -------------------------------
pdf_file <- file.path(output_dir, "demo_supfig3B.pdf")

pdf(pdf_file, width = 12, height = 10)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("Saved PDF: ", normalizePath(pdf_file, winslash = "/"), "\n", sep = "")
cat("\nSample-correlation heatmap workflow completed successfully.\n")
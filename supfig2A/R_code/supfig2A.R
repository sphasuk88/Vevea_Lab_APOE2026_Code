# supfig2A_organelle_marker_heatmap.R
# Purpose:
#   Generate organelle marker subset tables and a raw log2 heatmap
#   for 3-month male/female samples.
#
# Expected folder layout:
#   supfig2A/
#   ├── input/
#   │   ├── demo_main_dataset.csv
#   │   └── organelle_proteins_GN.csv
#   ├── output/
#   └── R_code/
#
# Reviewer instructions:
#   1. Run this script in R or RStudio.
#   2. When prompted, manually select:
#      - the main proteomics CSV file
#      - the organelle marker CSV file
#   3. The script will automatically determine the project folder and save all
#      output files to the output folder.

# -------------------------------
# Install/load required packages
# -------------------------------
required_packages <- c("pheatmap")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(pheatmap)

# -------------------------------
# Select input files manually
# -------------------------------
cat("Select MAIN proteomics dataset...\n")
input_file_main <- file.choose()
input_file_main <- normalizePath(input_file_main, winslash = "/", mustWork = TRUE)

cat("Select ORGANELLE marker dataset...\n")
input_file_markers <- file.choose()
input_file_markers <- normalizePath(input_file_markers, winslash = "/", mustWork = TRUE)

# -------------------------------
# Derive project and output paths
# -------------------------------
project_dir <- dirname(dirname(input_file_main))
output_dir <- file.path(project_dir, "output")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------
# Print resolved paths
# -------------------------------
cat("Main input file: ", input_file_main, "\n", sep = "")
cat("Marker input file: ", input_file_markers, "\n", sep = "")
cat("Project directory: ", project_dir, "\n", sep = "")
cat("Output directory: ", output_dir, "\n\n", sep = "")

# -------------------------------
# Load input files
# -------------------------------
ori_data <- read.csv(input_file_main, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
organelle_markers <- read.csv(input_file_markers, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

cat("Main dataset dimensions: ", nrow(ori_data), " rows x ", ncol(ori_data), " columns\n", sep = "")
cat("Marker dataset dimensions: ", nrow(organelle_markers), " rows x ", ncol(organelle_markers), " columns\n\n", sep = "")

# -------------------------------
# Check required columns
# -------------------------------
if (!("GN" %in% colnames(ori_data))) {
  stop("Required column 'GN' not found in the main proteomics dataset.")
}
if (!("GN" %in% colnames(organelle_markers))) {
  stop("Required column 'GN' not found in the organelle marker dataset.")
}

# -------------------------------
# Quick checks
# -------------------------------
cat("Case-insensitive match table:\n")
print(table(toupper(organelle_markers$GN) %in% toupper(ori_data$GN)))
cat("\n")

unmatched_markers <- setdiff(toupper(organelle_markers$GN), toupper(ori_data$GN))
if (length(unmatched_markers) > 0) {
  cat("Markers not matched in main dataset:\n")
  print(unmatched_markers)
  cat("\n")
}

# -------------------------------
# Select 3-month M/F sample columns
# -------------------------------
m3f3_cols <- grep("(_[MF]_0?3$)", colnames(ori_data), value = TRUE)

if (length(m3f3_cols) == 0) {
  stop("No 3-month male/female sample columns found matching '(_[MF]_0?3$)'.")
}

cat("3-month male/female columns used:\n")
print(m3f3_cols)
cat("\n")

# -------------------------------
# Match marker GNs to main data (case-insensitive)
# Preserve marker order
# -------------------------------
idx <- match(toupper(trimws(organelle_markers$GN)), toupper(trimws(ori_data$GN)))

# -------------------------------
# Build subset table
# -------------------------------
subset_df <- data.frame(
  GN = organelle_markers$GN,
  ori_data[idx, m3f3_cols, drop = FALSE],
  check.names = FALSE
)

out_csv <- file.path(output_dir, "demo_organelle_markers_MF_3mo_subset.csv")
write.csv(subset_df, out_csv, row.names = FALSE)
cat("Saved subset table: ", normalizePath(out_csv, winslash = "/"), "\n", sep = "")

# -------------------------------
# Log2-transform subset table
# -------------------------------
log2_subset_df <- subset_df

num_cols <- sapply(log2_subset_df, is.numeric)
log2_subset_df[, num_cols][log2_subset_df[, num_cols] <= 0] <- NA
log2_subset_df[, num_cols] <- log2(log2_subset_df[, num_cols])

out_log2_csv <- file.path(output_dir, "demo_organelle_markers_MF_3mo_subset_log2.csv")
write.csv(log2_subset_df, out_log2_csv, row.names = FALSE)
cat("Saved log2 subset table: ", normalizePath(out_log2_csv, winslash = "/"), "\n", sep = "")

# -------------------------------
# Build plot table with optional row labels
# -------------------------------
plot_df <- data.frame(
  GN = organelle_markers$GN,
  ori_data[idx, m3f3_cols, drop = FALSE],
  check.names = FALSE
)

if ("ProteinDescription" %in% colnames(ori_data)) {
  prot_desc <- ori_data$ProteinDescription[idx]
  short_desc <- gsub("\\sOS=.*$", "", ifelse(is.na(prot_desc), "", prot_desc))
  row_labels <- ifelse(
    nchar(short_desc) > 0,
    paste0(plot_df$GN, " — ", short_desc),
    plot_df$GN
  )
} else {
  row_labels <- plot_df$GN
}

# -------------------------------
# Prepare matrix for log2 transform
# -------------------------------
expression_data <- plot_df[, m3f3_cols, drop = FALSE]

expr_log2 <- expression_data
is_num <- sapply(expr_log2, is.numeric)
expr_log2[, is_num][expr_log2[, is_num] <= 0] <- NA
expr_log2[, is_num] <- log2(expr_log2[, is_num])

# -------------------------------
# Also save scaled matrix for records
# -------------------------------
scale_cols_by_min_max <- function(x) {
  xmax <- suppressWarnings(max(x, na.rm = TRUE))
  xmin <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(xmax) || !is.finite(xmin) || (xmax - xmin) == 0) {
    return(rep(NA_real_, length(x)))
  } else {
    return((x - xmin) / (xmax - xmin))
  }
}

scaled_mat <- as.data.frame(apply(expr_log2, 2, scale_cols_by_min_max), check.names = FALSE)
rownames(scaled_mat) <- plot_df$GN

scaled_csv <- file.path(output_dir, "demo_organelle_markers_MF_3mo_scaled_log2_matrix.csv")
write.csv(
  cbind(GN = rownames(scaled_mat), scaled_mat),
  scaled_csv,
  row.names = FALSE
)
cat("Saved scaled matrix CSV: ", normalizePath(scaled_csv, winslash = "/"), "\n", sep = "")

# -------------------------------
# Column annotation and ordering
# -------------------------------
sex_info <- gsub(".*_([MF])_.*", "\\1", colnames(expr_log2))
col_annotation <- data.frame(
  Sex = factor(sex_info, levels = c("M", "F")),
  stringsAsFactors = FALSE
)
rownames(col_annotation) <- colnames(expr_log2)

order_idx <- order(col_annotation$Sex)
mat_for_plot <- expr_log2[, order_idx, drop = FALSE]
col_annotation <- col_annotation[order_idx, , drop = FALSE]

sex_colors <- c(M = "darkgreen", F = "orange")

# -------------------------------
# Save raw log2 heatmap only
# -------------------------------
rownames(mat_for_plot) <- plot_df$GN

heatmap_pdf_raw <- file.path(output_dir, "demo_organelle_markers_MF_3mo_heatmap_log2_raw.pdf")

pheatmap(
  as.matrix(mat_for_plot),
  labels_row = row_labels,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Organelle Markers — 3-Month M/F (log2 intensities)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  legend = TRUE,
  na_col = "grey",
  annotation_col = col_annotation,
  annotation_colors = list(Sex = sex_colors),
  filename = heatmap_pdf_raw,
  width = 6,
  height = 8
)

cat("Saved raw log2 heatmap PDF: ", normalizePath(heatmap_pdf_raw, winslash = "/"), "\n", sep = "")
cat("\nsupfig2A workflow completed successfully.\n")
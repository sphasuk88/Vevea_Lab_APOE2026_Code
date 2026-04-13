# IgG_subtraction_clean.R
# Purpose:
#   Reorder proteomics sample columns, calculate the average of control columns,
#   subtract the control average from non-control samples, and save processed files.
#
# Reviewer instructions:
#   1. Run this script in R or RStudio.
#   2. When prompted, manually select the input CSV file.
#   3. The script will automatically determine the project folder and save all
#      output files to the output_and_processing folder.

# -------------------------------
# Select input file manually
# -------------------------------
input_file <- file.choose()
input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)

# -------------------------------
# Derive project and output paths
# -------------------------------
project_dir <- dirname(dirname(input_file))
output_dir <- file.path(project_dir, "output_and_processing")

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
# Reorder columns: first 3 fixed, then M, F, IR, Ctr
# -------------------------------
all_cols <- colnames(data)
ages <- c(3, 6, 12, 18, 24, 30)

get_cols <- function(sex, age, cols) {
  grep(paste0("_", sex, "_", age, "$"), cols)
}

m_idx <- unlist(lapply(ages, function(a) get_cols("M", a, all_cols)))
f_idx <- unlist(lapply(ages, function(a) get_cols("F", a, all_cols)))
ir_idx <- grep("_IR_IR$", all_cols)
ctr_idx <- grep("_Ctr_Ctr$", all_cols)

fixed_idx <- 1:3
final_order <- c(fixed_idx, m_idx, f_idx, ir_idx, ctr_idx)

data_sorted <- data[, final_order, drop = FALSE]

# -------------------------------
# Save grouped/sorted data
# -------------------------------
out_sorted <- file.path(output_dir, "demo_data_filtering_grpsort.csv")
write.csv(data_sorted, out_sorted, row.names = FALSE)

cat("Saved: ", normalizePath(out_sorted, winslash = "/"), "\n", sep = "")

# ===============================
# CONTROL AVERAGE & SUBTRACTION
# ===============================
sc <- colnames(data_sorted)

m_idx2 <- unlist(lapply(ages, function(a) grep(paste0("_M_", a, "$"), sc)))
f_idx2 <- unlist(lapply(ages, function(a) grep(paste0("_F_", a, "$"), sc)))
ir_idx2 <- grep("_IR_IR$", sc)
ctr_idx2 <- grep("_Ctr_Ctr$", sc)

fixed_idx2 <- 1:3

if (length(ctr_idx2) == 0) {
  stop("No control columns found matching '*_Ctr_Ctr'.")
}

# -------------------------------
# 1) Compute control average
# -------------------------------
Ctr_Avg <- rowMeans(data_sorted[, ctr_idx2, drop = FALSE], na.rm = TRUE)

data_with_CtrAvg <- cbind(
  data_sorted[, fixed_idx2, drop = FALSE],
  Ctr_Avg = Ctr_Avg,
  data_sorted[, c(m_idx2, f_idx2, ir_idx2, ctr_idx2), drop = FALSE]
)

# -------------------------------
# 2) Subtract Ctr_Avg from NON-control samples (M/F/IR)
#    Controls are dropped from subtraction output
# -------------------------------
non_ctr_idx2 <- c(m_idx2, f_idx2, ir_idx2)

sub_mat <- data_sorted[, non_ctr_idx2, drop = FALSE]
sub_mat[] <- lapply(sub_mat, as.numeric)

sub_mat <- sweep(sub_mat, 1, Ctr_Avg, FUN = "-")

final_minusCtr <- cbind(
  data_sorted[, fixed_idx2, drop = FALSE],
  Ctr_Avg = Ctr_Avg,
  sub_mat
)

# -------------------------------
# Optional: clean ProteinAccession
# Example: sp|P97490|ADCY8_MOUSE -> P97490
# -------------------------------
if ("ProteinAccession" %in% colnames(final_minusCtr)) {
  final_minusCtr$ProteinAccession <- sub("^.*\\|(.*)\\|.*$", "\\1", final_minusCtr$ProteinAccession)
}

# -------------------------------
# Add IR < 0 indicator columns
# -------------------------------
ir_col_names <- colnames(sub_mat)[grepl("_IR_IR$", colnames(sub_mat))]
for (nm in ir_col_names) {
  final_minusCtr[[paste0(nm, " < 0")]] <- ifelse(final_minusCtr[[nm]] < 0, 1L, 0L)
}

# -------------------------------
# Add b123_Neg
# Count negative values across all non-IR animal samples (M + F only)
# -------------------------------
mf_cols <- colnames(sub_mat)[!grepl("_IR_IR$", colnames(sub_mat))]
final_minusCtr$b123_Neg <- rowSums(final_minusCtr[, mf_cols, drop = FALSE] < 0, na.rm = TRUE)

out_sub <- file.path(output_dir, "demo_data_filtering_grpsort_IgGsub.csv")
write.csv(final_minusCtr, out_sub, row.names = FALSE)

cat("Saved: ", normalizePath(out_sub, winslash = "/"), "\n", sep = "")
cat("\nIgG subtraction workflow completed successfully.\n")
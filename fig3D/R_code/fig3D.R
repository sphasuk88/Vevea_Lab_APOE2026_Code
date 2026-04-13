# fig3D_log2FC_scatter.R
# Purpose:
#   Generate the Fig3D scatter plot comparing male and female log2 fold changes
#   for the 30-month versus 3-month comparisons.
#
# Expected folder layout:
#   fig3D/
#   ├── input/
#   │   └── female-male_30mo vs 3mo.csv
#   ├── output/
#   └── R_code/
#
# Reviewer instructions:
#   1. Run this script in R or RStudio.
#   2. When prompted, manually select the input CSV file from the input folder.
#   3. The script will automatically determine the project folder and save the
#      output PDF to the output folder.

# -------------------------------
# Install/load required packages
# -------------------------------
required_packages <- c("dplyr", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(ggplot2)

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
# Load data
# -------------------------------
all_data <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input data dimensions: ", nrow(all_data), " rows x ", ncol(all_data), " columns\n\n", sep = "")

# -------------------------------
# Check required columns
# -------------------------------
required_cols <- c(
  "GN",
  "Log2Fold(M_30/M_3)",
  "Log2Fold(F_30/F_3)",
  "FDR(M_30/M_3)",
  "FDR(F_30/F_3)"
)

missing_cols <- setdiff(required_cols, colnames(all_data))
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The following required columns are missing from the input file:\n",
      paste(missing_cols, collapse = "\n")
    )
  )
}

# -------------------------------
# Prepare plotting data
# -------------------------------
plot_df <- all_data %>%
  select(
    GN,
    logFC_M = `Log2Fold(M_30/M_3)`,
    logFC_F = `Log2Fold(F_30/F_3)`,
    FDR_M   = `FDR(M_30/M_3)`,
    FDR_F   = `FDR(F_30/F_3)`
  ) %>%
  filter(
    !is.na(logFC_M), !is.na(logFC_F),
    !is.na(FDR_M), !is.na(FDR_F)
  ) %>%
  mutate(
    Sig_M = FDR_M < 0.05,
    Sig_F = FDR_F < 0.05,
    SigGroup = dplyr::case_when(
      Sig_M & Sig_F ~ "Both",
      Sig_M & !Sig_F ~ "Male",
      !Sig_M & Sig_F ~ "Female",
      TRUE ~ "Not_Significant"
    ),
    FDR_combined = pmin(FDR_M, FDR_F),
    logFDR = -log10(FDR_combined),
    SigGroup = factor(SigGroup, levels = c("Not_Significant", "Male", "Female", "Both"))
  )

cat("Plotting data dimensions: ", nrow(plot_df), " rows x ", ncol(plot_df), " columns\n\n", sep = "")

# -------------------------------
# Manual label positions
# -------------------------------
label_positions <- data.frame(
  GN = c("C1qa", "C1qb", "C1qc", "Zc3hav1", "Apobr", "Lgals3", "Apoe"),
  x  = c(1.2, 1.7, 1.5, 1.5, 0.5, 2.2, -0.5),
  y  = c(3.5, 3.2, 4.5, -1.5, 2.5, 2.7, 3.0),
  stringsAsFactors = FALSE
)

# Only draw segments for genes that exist in plot_df
lines_to_labels <- plot_df %>%
  inner_join(label_positions, by = "GN")

# -------------------------------
# Build plot
# -------------------------------
log2FC_plot <- ggplot(plot_df, aes(x = logFC_M, y = logFC_F)) +
  
  # Non-significant points first
  geom_point(
    data = dplyr::filter(plot_df, SigGroup == "Not_Significant"),
    aes(size = logFDR, color = SigGroup),
    alpha = 0.45
  ) +
  
  # Significant points on top
  geom_point(
    data = dplyr::filter(plot_df, SigGroup != "Not_Significant"),
    aes(size = logFDR, color = SigGroup),
    alpha = 0.85
  ) +
  
  # Connecting lines from points to labels
  geom_segment(
    data = lines_to_labels,
    aes(x = logFC_M, y = logFC_F, xend = x, yend = y),
    color = "gray30", linewidth = 0.3
  ) +
  
  # Gene labels
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = GN),
    hjust = 0,
    color = "black",
    size = 2.8
  ) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Colors and sizes
  scale_color_manual(
    name = "Significance (FDR < 0.05)",
    values = c(
      "Not_Significant" = rgb(169, 169, 169, maxColorValue = 255, alpha = 160),
      "Male"            = rgb(0, 180, 0, maxColorValue = 255, alpha = 200),
      "Female"          = rgb(255, 140, 0, maxColorValue = 255, alpha = 200),
      "Both"            = rgb(0, 100, 255, maxColorValue = 255, alpha = 200)
    )
  ) +
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(1.5, 6),
    breaks = c(1, 2, 3, 5)
  ) +
  
  # Labels and theme
  labs(
    x = "Log2 Fold Change (Male 30 vs 3 months)",
    y = "Log2 Fold Change (Female 30 vs 3 months)",
    title = "Male vs Female Log2 Fold Change (30 vs 3 months)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title   = element_text(size = 8, face = "bold"),
    axis.title   = element_text(size = 8),
    axis.text    = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 5),
    legend.spacing.y  = unit(0.1, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.position   = "right",
    legend.margin     = margin(l = 0, r = 0, t = 0, b = 0),
    legend.box.margin = margin(0, 0, 0, -10)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# -------------------------------
# Print plot
# -------------------------------
print(log2FC_plot)

# -------------------------------
# Save plot
# -------------------------------
pdf_file <- file.path(output_dir, "demo_Log2FC_30vs3_MaleFemale2.pdf")

ggsave(
  filename = pdf_file,
  plot     = log2FC_plot,
  width    = 6,
  height   = 3,
  units    = "in",
  dpi      = 300,
  device   = cairo_pdf,
  family   = "sans"
)

cat("Saved PDF: ", normalizePath(pdf_file, winslash = "/"), "\n", sep = "")
cat("\nFig3D scatter plot workflow completed successfully.\n")
# supfig2A

## Description
This folder contains an R script used to generate organelle marker subset tables and a raw log2 heatmap for 3-month male and female samples from a proteomics CSV file and an organelle marker list.
----------------------------------------------------------------
## Folder structure

```text
supfig2A/
├── input/
│   ├── demo_main_dataset.csv
│   └── organelle_proteins_GN.csv
├── output/
└── R_code/
    └── supfig2A.R
```

----------------------------------------------------------------

**R version:**

This script was tested with: R version 4.5.2

----------------------------------------------------------------

**Packages used in this script:**

pheatmap  

The script automatically checks for and installs missing packages when run.

----------------------------------------------------------------

**How to run the demo:**

1. Open RStudio.

2. Open the script:

3. R_code/supfig2A.R

4. Run the script.

5. When prompted, choose the main proteomics dataset (demo_supfig2A.csv).

6. Then choose the organelle marker dataset (organelle_proteins_GN.csv).

7. The script will automatically:

8. detect the project folder

9. create the output folder if needed

10. match organelle marker gene names to the main dataset

11. extract 3-month male and female sample columns

12. save the subset table

13. save the log2-transformed subset table

14. save the scaled log2 matrix used for records

15. generate the raw log2 heatmap PDF

16. save all output files in the output folder

17. print the selected input file paths, project directory, output directory, and input data dimensions

----------------------------------------------------------------

**Expected output:**

The script produces the following output files:

> demo_organelle_markers_MF_3mo_subset.csv

Subset table containing the selected organelle marker genes and 3-month male/female sample columns.

> demo_organelle_markers_MF_3mo_subset_log2.csv

Log2-transformed version of the subset table.

> demo_organelle_markers_MF_3mo_scaled_log2_matrix.csv

Column-wise min-max scaled matrix derived from the log2-transformed subset, saved for records and reproducibility.

> demo_organelle_markers_MF_3mo_heatmap_log2_raw.pdf

Heatmap of raw log2 intensities for organelle marker proteins in 3-month male and female samples.

All files are saved in the output folder.

----------------------------------------------------------------

**Input file requirements:**

The main proteomics CSV file must contain:

GN  

and sample-expression columns matching this naming pattern:

(_[MF]_0?3$)

The organelle marker CSV file must contain:

GN

The script uses these files to:
- perform case-insensitive matching of organelle marker gene names to the main dataset
- extract 3-month male and female sample columns
- generate subset tables
- apply safe log2 transformation
- create the organelle marker heatmap

----------------------------------------------------------------

**The script will also print:**

the selected input file paths

the project directory

the output directory

the main dataset dimensions

the marker dataset dimensions

the case-insensitive match summary

the list of unmatched markers, if any

the list of 3-month male/female columns used

----------------------------------------------------------------

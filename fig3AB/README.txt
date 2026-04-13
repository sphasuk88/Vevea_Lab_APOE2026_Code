# fig3AB

## Description
This folder contains an R script used to generate the clustered heatmap and downstream cluster summary plots for Fig3A/3B from a proteomics CSV file.
----------------------------------------------------------------
## Folder structure

```text
fig3AB/
├── input/
│   └── demo_fig3AB.csv
├── output/
└── R_code/
    └── fig3AB.R
```

----------------------------------------------------------------

**R version:**

This script was tested with: R version 4.5.2

----------------------------------------------------------------

**Packages used in this script:**

pheatmap  
ggplot2  
dplyr  
tidyr  

The script automatically checks for and installs missing packages when run.

----------------------------------------------------------------

**How to run the demo:**

1. Open RStudio.

2. Open the script:

3. R_code/fig3AB.R

4. Run the script.

5. When prompted, choose the file:

6. input/demo_fig3AB.csv

7. The script will automatically:

8. detect the project folder

9. create the output folder if needed

10. generate the clustered heatmap PDF

11. assign proteins to clusters

12. save the clustered DEP table

13. calculate median and SEM values by cluster, sex, and age

14. generate cluster summary plots

15. save all output files in the output folder

16. print the selected input file path, project directory, output directory, input data dimensions, and filtered row count

----------------------------------------------------------------

**Expected output:**

The script produces the following output files:

> heatmap.pdf

Clustered heatmap of proteins with FDR < 0.05.

> clusters.csv

DEP table with assigned cluster labels.

> cluster_gender_age_medianSEM.csv

Summary table containing median and SEM values by cluster, sex, and age.

> Cluster_1_medianSEM_plot.pdf  
> Cluster_2_medianSEM_plot.pdf  
> Cluster_3_medianSEM_plot.pdf
> Cluster_4_medianSEM_plot.pdf

Cluster trend plots showing normalized median expression across age for male and female samples.

All files are saved in the output folder.

----------------------------------------------------------------

**Input file requirements:**

The input CSV file must contain:

1. a column named:

GN  
FDR

2. numeric expression columns in positions 4:51 for heatmap generation

3. sample-expression columns matching this naming pattern for downstream scaling and summary analysis:

b[1-3]_s[0-9]+_[MF]_[0-9]+

The script uses these columns to:
- filter proteins with FDR < 0.05
- generate the clustered heatmap
- assign cluster labels
- perform row-wise z-score normalization
- summarize expression trends by cluster, sex, and age

----------------------------------------------------------------

**The script will also print:**

the selected input file path

the project directory

the output directory

the input data dimensions

the number of rows with FDR < 0.05

----------------------------------------------------------------

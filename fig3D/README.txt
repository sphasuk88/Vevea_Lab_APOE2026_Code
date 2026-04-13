# fig3D

## Description
This folder contains an R script used to generate the Fig3D scatter plot comparing male and female log2 fold changes for the 30-month versus 3-month proteomic comparisons.
----------------------------------------------------------------
## Folder structure

```text
fig3D/
├── input/
│   └── demio_fig3D.csv 
├── output/
└── R_code/
    └── fig3D.R
```

----------------------------------------------------------------

**R version:**

This script was tested with: R version 4.5.2

----------------------------------------------------------------

**Packages used in this script:**

dplyr  
ggplot2  

The script automatically checks for and installs missing packages when run.

----------------------------------------------------------------

**How to run the demo:**

1. Open RStudio.

2. Open the script:

3. R_code/fig3D.R

4. Run the script.

5. When prompted, choose the file:

6. input/demo_fig3D.csv

7. The script will automatically:

8. detect the project folder

9. create the output folder if needed

10. generate the scatter plot PDF

11. save the output file in the output folder

12. print the selected input file path, project directory, output directory, and input data dimensions

----------------------------------------------------------------

**Expected output:**

The script produces a scatter plot comparing male and female log2 fold changes for the 30-month versus 3-month comparisons.

> demo_fig3D.pdf

This file is saved in the output folder.

----------------------------------------------------------------

**Input file requirements:**

The input CSV file must contain the following columns:

GN  
Log2Fold(M_30/M_3)  
Log2Fold(F_30/F_3)  
FDR(M_30/M_3)  
FDR(F_30/F_3)  

The script uses these columns to:
- assign significance groups
- calculate point size using the combined FDR
- generate the male-versus-female log2 fold-change scatter plot

----------------------------------------------------------------

**The script will also print:**

the selected input file path

the project directory

the output directory

the input data dimensions

----------------------------------------------------------------

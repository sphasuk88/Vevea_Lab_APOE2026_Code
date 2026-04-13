# IgG_subtraction and IR filtering

## Description
This folder contains an R script used to reorder proteomics sample columns, calculate the average of control columns, subtract the control average from non-control samples, and generate processed output files for downstream analysis.
----------------------------------------------------------------
## Folder structure

```text
IgG_subtraction/
├── input/demo_data_filtering.csv
├── output_and_processing/
└── R_script/
    └── IgG_subtraction_clean.R
```

----------------------------------------------------------------

**R version:**

This script was tested with: R version 4.5.2

----------------------------------------------------------------

**Packages used in this script:**

No additional R packages are required beyond base R.

----------------------------------------------------------------

**How to run the demo:**

1. Open RStudio.

2. Open the script:

3. R_script/data_filtering.R

4. Run the script.

5. When prompted, choose the input CSV file (demo_data_filtering.csv) from the input folder.

6. The script will automatically:

7. detect the project folder

8. create the output_and_processing folder if needed

9. reorder sample columns by group

10. calculate the average of control columns

11. subtract the control average from non-control samples

12. add IR < 0 indicator columns

13. calculate the b123_Neg column

14. save all output files in the output_and_processing folder

15. print the selected input file path, project directory, output directory, and input data dimensions

----------------------------------------------------------------

**Expected output:**

The script produces the following output files:

> demo_data_filtering_grpsort.csv

Grouped and reordered version of the input dataset.

> demo_data_filtering_grpsort_IgGsub.csv

IgG-subtracted dataset containing the control-average column, subtracted non-control sample columns, IR < 0 indicator columns, and b123_Neg.

All files are saved in the output_and_processing folder.

----------------------------------------------------------------

**Input file requirements:**

The input CSV file must contain:

1. metadata columns in the first 3 positions

2. sample-expression columns matching the following naming patterns:

_M_3, _M_6, _M_12, _M_18, _M_24, _M_30  
_F_3, _F_6, _F_12, _F_18, _F_24, _F_30  
_IR_IR  
_Ctr_Ctr

The script uses these columns to:
- reorder samples by group
- calculate the average control signal
- subtract control values from non-control samples
- identify negative IR values
- count negative values across male and female sample columns

----------------------------------------------------------------

**The script will also print:**

the selected input file path

the project directory

the output directory

the input data dimensions

----------------------------------------------------------------

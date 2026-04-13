Protein Isoform Selection Workflow
=================================

Overview
--------
This workflow identifies a representative protein isoform for each gene name (GN)
from a tab-delimited protein quantification file. The workflow reads the input file,
computes helper metrics used for isoform ranking, selects one representative accession
per gene, and writes an Excel output file with a final annotation column named:

    representative_isoform

Rows labeled `1` in this column are the selected representative accession for that gene.
Rows labeled `0` are not selected.


Project structure
-----------------
The workflow expects the following folder structure:

    project_folder/
    ├── input/
    │   └── id_uni_prot_quan.txt
    ├── output/
    ├── script/
    │   ├── main.py
    │   ├── ProteinIsoformLib.py
    │   └── PhosphoSiteLib.py
    ├── cmd_Isoform.bat
    ├── Isoform_input.params
    └── README.txt


Software requirements
---------------------
- Python 3
- pandas
- numpy
- openpyxl

Install the required Python packages with:

    py -m pip install pandas numpy openpyxl


Input file
----------
Place the tab-delimited input file in:

    .\input\id_uni_prot_quan.txt


Parameter file
--------------
The parameter file used by the workflow is:

    Isoform_input.params

It should contain:

    #--- input for Protein Isoform generation --- #########

    #full file of id_uni_prot_quan.txt
    batch_quant_path = .\input\id_uni_prot_quan.txt


How to run
----------
Option 1: Run from Windows Explorer
    Double-click:

        cmd_Isoform.bat

Option 2: Run from Command Prompt
    From the project folder, run:

        py .\script\main.py Isoform_input.params


Output
------
The workflow generates:

    .\output\id_uni_prot_quan_W1.xlsx

The output file contains the original quantification columns and one additional column:

    representative_isoform


Workflow summary
----------------
The workflow performs the following steps:

1. Read the tab-delimited protein quantification file
2. Compute helper metrics for each row
3. Remove rows with missing GN before representative isoform selection
4. Group rows by GN
5. Select one representative protein accession per GN
6. Map the selected accession back to all rows
7. Create the final column `representative_isoform`
8. Write the annotated table to Excel


Helper metrics used for selection
---------------------------------
The script computes the following values for each row before grouping by GN:

1. nbatch
   Number of PSM columns with a non-zero value for that row.

2. average_PSMs
   Mean of all columns whose names contain the string `PSM`.

3. average_Intensity
   Mean of all columns whose names contain the string `sig`.

Additional helper annotations are also computed internally:
- ProtType
- KRT

These are used during processing but are not retained in the final Excel output.


Representative isoform selection rules
--------------------------------------
Representative isoforms are selected independently for each gene name (GN).

Decision order:

1. Single-accession case
   If a GN maps to only one protein accession, that accession is selected directly.

2. Highest nbatch
   If a GN maps to multiple accessions, the accession with the highest `nbatch`
   is preferred first.

3. Canonical accession priority among tied accessions
   If multiple accessions tie at the highest `nbatch`, the tied subset is evaluated
   further. Canonical accessions are prioritized, using the exact rule implemented
   in the script:
   - accession contains `sp|`
   - accession does not contain `-`

4. Multiple canonical accessions
   If more than one canonical accession remains:
   - choose the accession with the highest `average_PSMs`
   - if still tied, choose the accession with the highest `average_Intensity`

5. No canonical accession available
   If no canonical accession is present in the tied subset:
   - if only one isoform remains, select it
   - if multiple isoforms remain, prefer the lowest isoform number after `-`
   - if still tied, choose the accession with the highest `average_PSMs`
   - if still tied again, choose the accession with the highest `average_Intensity`


Definition of representative_isoform
------------------------------------
After one accession is selected for each GN, the selected accession is mapped back
to all rows with that GN.

The final column is assigned as follows:

- `1` if `Protein Accession #` equals the selected accession for that GN
- `0` otherwise


Important implementation notes
------------------------------
1. Rows with missing GN are excluded before representative isoform selection.

2. The final exported Excel file contains:
   - the original input columns
   - the final column `representative_isoform`

   It does not retain intermediate helper columns such as:
   - nbatch
   - average_PSMs
   - average_Intensity
   - ProtType
   - KRT
   - accession_uniq

3. Canonical accession priority is based on the exact accession-pattern rule used in
   the code, not on an external UniProt annotation lookup.

4. In the isoform-number comparison step, the implementation uses the first character
   after `-`. For example:

       P12345-2

   is interpreted as isoform 2.

   Multi-digit isoform numbers may therefore not be interpreted as full integers.


Suggested methods text
----------------------
Representative isoforms were selected per gene name (GN) from the protein
quantification table. For each GN, the protein accession with the highest `nbatch`
(number of non-zero PSM columns) was prioritized. Ties were resolved by prioritizing
canonical Swiss-Prot-style accessions, defined in the script as entries containing
`sp|` and not containing `-`. Among tied canonical accessions, the accession with the
highest average PSM was selected, followed by the highest average intensity when
necessary. If no canonical accession was available, isoforms were prioritized by the
lowest isoform number, with remaining ties resolved by highest average PSM and then
highest average intensity. Rows matching the selected accession for each GN were
annotated with `representative_isoform = 1`, and all other rows were annotated as `0`.


Troubleshooting
---------------
If the workflow does not run successfully:

1. Confirm that Python is installed:

       py --version

2. Confirm that required packages are installed:

       py -m pip install pandas numpy openpyxl

3. Confirm that the input file exists at:

       .\input\id_uni_prot_quan.txt

4. Confirm that `Isoform_input.params` contains the correct path:

       batch_quant_path = .\input\id_uni_prot_quan.txt

5. Try running manually from Command Prompt:

       py .\script\main.py Isoform_input.params

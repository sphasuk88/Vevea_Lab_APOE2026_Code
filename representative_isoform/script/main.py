import pandas as pd
import numpy as np
import os,sys
from os.path import dirname
from ProteinIsoformLib import *
from PhosphoSiteLib import *
import warnings

if len(sys.argv)<2:
    sys.exit("  Please use command (python main.py W1P1input.params) to re-run")
params_file = sys.argv[1]
params = Get_Params(params_file)

in1 = params["batch_quant_path"] #"./combined_norm_uni_prot.txt"

print ("  Reading the txt file {}".format(in1))

df_quant = pd.read_csv(in1, delimiter="\t", skiprows=return_rows_head(in1, "\t"))


result_path = os.getcwd()+"/output"
if not os.path.exists(result_path):
    os.makedirs(result_path)

out1 = result_path+"/id_uni_prot_quan_W1.xlsx"


quan_cols=list(df_quant.columns)
tail_cols = ["representative_isoform"]


print ("  Computing nbatch, average psms, and average intensity")

dfCount = df_quant.filter(like = "PSM", axis=1)
df_quant["nbatch"] = len(dfCount.columns) - (dfCount == 0).astype(int).sum(axis=1)
df_quant["average_PSMs"]=df_quant.apply(average_PSM, df=df_quant, axis=1)
df_quant["average_Intensity"]=df_quant.apply(average_intensity, df=df_quant, axis=1)

print ("  Computing ProtType, KRT")

df_quant["ProtType"] = df_quant.apply(lambda x: get_uniprot(x["Protein Accession #"]), axis=1)
df_quant["KRT"] = df_quant.apply(lambda x: get_KRT(str(x["GN"])), axis=1)

print ("  Cleaning empty gene name")

df_quant_clean=df_quant.dropna(subset=["GN"])

print ("  Grouping the column values based on unique GN")

df_quant_groupGN=df_quant_clean.groupby(["GN"])[["Protein Accession #","nbatch","average_PSMs","average_Intensity"]].agg(lambda x: x.tolist()).reset_index()

print ("  Selecting correct isoform for each row (Please be patient)")

df_quant_groupGN["accession"]=df_quant_groupGN.apply(select_final_valtest_final, axis=1)

print ("  Forming dictionary of correct isoform with unique protein accession")

representative_isoform_dict = dict(zip(df_quant_groupGN.GN, df_quant_groupGN.accession))

print ("  Computing accession_uniq")

df_quant["accession_uniq"] = df_quant["GN"].map(representative_isoform_dict)

print ("  Computing representative_isoform")

df_quant["representative_isoform"] = df_quant.apply(lambda x: "1" if x["Protein Accession #"] == x.accession_uniq else "0", axis=1)

print ("  Outputing the new excel file {}".format(out1))
df_final = df_quant[quan_cols+tail_cols]
df_final.to_excel(out1, index=None)

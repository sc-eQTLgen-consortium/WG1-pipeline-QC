#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--ind_coupling", required=False, default=None, type=str, help="")
parser.add_argument("--main_dir", required=True, type=str, help="")
parser.add_argument("--threshold", required=False, default=0.02, type=float, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import pandas as pd

################################################################################

print("Loading data")

poolsheet = pd.read_csv(args.poolsheet, sep="\t", dtype=str)

man_select_df = None
for file_extension in ["selfSM", "depthSM", "bestSM", "selfRG", "depthRG", "bestRG"]:
    print("  Combining *.{} files".format(file_extension))
    df_list = []
    for _, row in poolsheet.iterrows():
        inpath = os.path.join(args.main_dir, row["Pool"], "verifybamid", "genoCheck." + file_extension)
        if not os.path.exists(inpath):
            continue
        df_list.append(pd.read_csv(inpath, sep="\t", header=0, index_col=None))

    print("\tLoaded {:,} files".format(len(df_list)))
    if len(df_list) == 0:
        exit()

    df = pd.concat(df_list, axis=0)

    outfile = "genoCheck." + file_extension
    df.to_csv(os.path.join(args.out, "CombinedResults", outfile), sep="\t", header=True, index=False)
    print("  Saved {}".format(outfile))

    if file_extension == "bestSM":
        man_select_df = df[["#SEQ_ID", "CHIP_ID", "FREEMIX", "CHIPMIX"]].copy()
    
    del df

man_select_df["#SEQ_ID"] = man_select_df["#SEQ_ID"].astype(str)
if args.ind_coupling is not None:
    print("Loading individual coupling file")
    ind_coupling_df = pd.read_csv(args.ind_coupling, sep="\t", dtype=str)
    print("\tLoaded {:,} row".format(ind_coupling_df.shape[0]))

    ind_coupling_df.columns = ["#SEQ_ID", "ASSIGNMENT"]
    man_select_df = man_select_df.merge(ind_coupling_df, on="#SEQ_ID", how="inner")
else:
    man_select_df["ASSIGNMENT"] = man_select_df["NA"]

print("\tFilling in recommendations using threshold == {}".format(args.threshold))
man_select_df["NOTE"] = ""
man_select_df.loc[(man_select_df["CHIPMIX"] > args.threshold) | (man_select_df["FREEMIX"] > args.threshold), "NOTE"] = "CONTAMINATION"
man_select_df.loc[(man_select_df["CHIPMIX"] > (1 - args.threshold)) & (man_select_df["FREEMIX"] < args.threshold), "NOTE"] = "SAMPLE SWAP"
man_select_df.loc[man_select_df["CHIPMIX"] < args.threshold, "NOTE"] = "USE CHIP_ID"

man_select_df["MATCH"] = False
man_select_df.loc[man_select_df["ASSIGNMENT"] == man_select_df["CHIP_ID"], "MATCH"] = True

man_select_df["UPDATE/REMOVE/KEEP"] = ""
man_select_df.loc[(man_select_df["NOTE"] == "USE CHIP_ID") & man_select_df["MATCH"], "UPDATE/REMOVE/KEEP"] = "KEEP"
print("\tFound {:,} issues out of {:,} samples".format((man_select_df["UPDATE/REMOVE/KEEP"] == "").sum(), man_select_df.shape[0]))

outfile = "verifyBamID_manual_selection.tsv"
man_select_df.to_csv(os.path.join(args.out, "manual_selection", outfile), sep="\t", header=True, index=False)
print("  Saved {}".format(outfile))
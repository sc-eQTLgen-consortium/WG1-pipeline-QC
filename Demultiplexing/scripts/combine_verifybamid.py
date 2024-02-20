#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--ind_coupling", required=False, default=None, type=str, help="")
parser.add_argument("--main_dir", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import pandas as pd

################################################################################

print("Loading data")

poolsheet = pd.read_csv(args.poolsheet, sep="\t", dtype=str)
df_list = []
for _, row in poolsheet.iterrows():
    inpath = os.path.join(args.main_dir, row["Pool"], "verifybamid", "genoCheck.bestSM")
    if not os.path.exists(inpath):
        continue
    df_list.append(pd.read_csv(inpath, sep="\t", header=0, index_col=None))

print("\tLoaded {:,} files".format(len(df_list)))
if len(df_list) == 0:
    exit()

df = pd.concat(df_list, axis=0)
df["#SEQ_ID"] = df["#SEQ_ID"].astype(str)

if args.ind_coupling is not None:
    print("Loading individual coupling file")
    ind_coupling_df = pd.read_csv(args.ind_coupling, sep="\t", dtype=str)
    print("\tLoaded {:,} row".format(ind_coupling_df.shape[0]))

    ind_coupling_df.columns = ["#SEQ_ID", "ASSIGNMENT"]
    df = df.merge(ind_coupling_df, on="#SEQ_ID", how="inner")

    df["SAMPLESWAP"] = df["ASSIGNMENT"] != df["CHIP_ID"]
    df["UPDATE/REMOVE/KEEP"] = ""
    df.loc[~df["SAMPLESWAP"], "UPDATE/REMOVE/KEEP"] = "KEEP"
    print("\tFound {} sample swaps out of {} samples".format(df["SAMPLESWAP"].sum(), df.shape[0]))
else:
    df["ASSIGNMENT"] = df["NA"]
    df["UPDATE/REMOVE/KEEP"] = ""

outfile = "verifyBamID_manual_selection.tsv"
df.to_csv(os.path.join(args.out, outfile), sep="\t", header=True, index=False)
print("  Saved {}".format(outfile))
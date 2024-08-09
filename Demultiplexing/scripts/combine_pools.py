#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--main_dir", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import pandas as pd

################################################################################

print("Loading data")

poolsheet = pd.read_csv(args.poolsheet, sep="\t", dtype=str)

data = {}
outfiles = ["combined_results.tsv.gz", "combined_results_w_combined_assignments.tsv.gz", "Final_Assignments_demultiplexing_doublets.tsv.gz"]
for outfile in outfiles:
    print("  Processing {}".format(outfile))

    pool_df_list = []
    for _, row in poolsheet.iterrows():
        inpath = os.path.join(args.main_dir, row["Pool"], "CombinedResults", outfile)
        if not os.path.exists(inpath):
            continue
        pool_df = pd.read_csv(inpath, sep="\t", header=0, index_col=None)
        pool_df.insert(0, "Pool", row["Pool"])
        pool_df["Barcode"] = pool_df["Barcode"].str.split('-').str[0] + "_" + pool_df["Pool"]
        pool_df_list.append(pool_df)

    if len(pool_df_list) == 0:
        print("\tNo input files found")
        continue

    df = pd.concat(pool_df_list, axis=0)
    if df.shape[0] == 0:
        print("\tNo data in input files")
        continue
    df.to_csv(os.path.join(args.out, outfile), sep="\t", header=True, index=False, compression="gzip")
    print("\tSaved combined {} with shape: {}".format(outfile, df.shape))

outfiles = ["combined_results_summary.tsv.gz", "combined_results_demultiplexing_summary.tsv.gz"]
for outfile in outfiles:
    print("  Processing {}".format(outfile))

    pool_df_list = []
    methods = set()
    for _, row in poolsheet.iterrows():
        inpath = os.path.join(args.main_dir, row["Pool"], "CombinedResults", outfile)
        if not os.path.exists(inpath):
            continue
        pool_df = pd.read_csv(inpath, sep="\t", header=0, index_col=None)
        methods.update(pool_df.columns)
        pool_df = pool_df.rename(columns={"N": row["Pool"]})
        pool_df_list.append(pool_df)
    df = pd.concat(pool_df_list, axis=0).fillna(0)
    methods.remove("N")

    # Group overlapping values and sort by the total.
    df = df.groupby(list(methods)).sum().reset_index(drop=False)
    df["N"] = df.loc[:, [col for col in df.columns if col not in methods]].sum(axis=1)
    df.sort_values(by="N", ascending=False, inplace=True)

    if df.shape[0] == 0:
        print("\tNo data in input files")
        continue

    df.to_csv(os.path.join(args.out, outfile), sep="\t", header=True, index=False, compression="gzip")
    print("\tSaved combined {} with shape: {}".format(outfile, df.shape))

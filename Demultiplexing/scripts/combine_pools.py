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

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import pandas as pd
import gzip


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


print("Loading data")
poolsheet = pd.read_csv(args.poolsheet, sep="\t", dtype=str)
npools = poolsheet.shape[0]

outfiles = ["combined_results.tsv.gz", "combined_results_w_combined_assignments.tsv.gz", "Final_Assignments_demultiplexing_doublets.tsv.gz"]
for outfile in outfiles:
    print("  Processing {}".format(outfile))

    fhout = gzopen(os.path.join(args.out, outfile), 'w')
    i = 0
    nrows = 0
    header = None
    ncolumns = None
    for i, row in poolsheet.iterrows():
        if i == 0 or i % 5 == 0:
            print("\tParsed {:,} / {:,} pools".format(i + 1, npools), end='\r')

        # Check if the output exists.
        inpath = os.path.join(args.main_dir, row["Pool"], "CombinedResults", outfile)
        if not os.path.exists(inpath):
            print("\tError, could not find output file for pool '{}'".format(row["Pool"]))
            exit()

        # Check how many barcodes we should have.
        j = 0
        with gzopen(row["Barcodes"], 'r') as f:
            for j, _ in enumerate(f):
                pass
        f.close()
        nbarcodes_pool = j + 1

        # Load the combined results file.
        nrows_pool = 0
        fhin = gzopen(inpath, 'r')
        for j, line in enumerate(fhin):
            values = line.rstrip("\n").split("\t")

            # Make sure the number of columns is always identical.
            if ncolumns is None:
                ncolumns = len(values)
            if len(values) != ncolumns:
                print("\tError, unexpected number of columns in {}".format(inpath))
                exit()

            # Handle the header.
            if j == 0:
                # Save which column is located where.
                indices = {column: index for index, column in enumerate(values)}

                # Only write header for the first file.
                if header is None:
                    header = values
                    fhout.write("\t".join(["Pool"] + values) + "\n")
                if values != header:
                    print("Error, all files must have the same column order.")
                    exit()

                continue

            # Reformat the barcodes value by adding the pool name.
            values[indices["Barcode"]] = values[indices["Barcode"]].split("-")[0] + "_" + row["Pool"]

            # Save line.
            fhout.write("\t".join([row["Pool"]] + values) + "\n")
            nrows_pool += 1
        fhin.close()

        # Check that the number of rows written matches the number of barcodes.
        if nrows_pool != nbarcodes_pool:
            print("\tError, number of rows does not match number of barcodes for pool '{}'. Barcodes file has {:,} barcodes but results file has {:,}.".format(row["Pool"], nbarcodes_pool, nrows_pool - 1))
            exit()

        nrows += nrows_pool
    fhout.close()

    print("\tParsed {:,} / {:,} pools".format(i + 1, npools), end='\n')
    print("\tSaved combined {} with shape: ({}, {})".format(outfile, nrows + 1, ncolumns))

# This is not done line by line since it has some complex merging that is more difficult.
outfiles = ["combined_results_summary.tsv.gz", "combined_results_demultiplexing_summary.tsv.gz"]
for outfile in outfiles:
    print("  Processing {}".format(outfile))

    pool_df_list = []
    methods = set()
    i = 0
    for i, row in poolsheet.iterrows():
        if i == 0 or i % 5 == 0:
            print("\tParsed {:,} / {:,} pools".format(i + 1, npools), end='\r')

        inpath = os.path.join(args.main_dir, row["Pool"], "CombinedResults", outfile)
        if not os.path.exists(inpath):
            continue
        pool_df = pd.read_csv(inpath, sep="\t", header=0, index_col=None)
        methods.update(pool_df.columns)
        pool_df = pool_df.rename(columns={"N": row["Pool"]})
        pool_df_list.append(pool_df)
    print("\tParsed {:,} / {:,} pools".format(i + 1, npools), end='\n')

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

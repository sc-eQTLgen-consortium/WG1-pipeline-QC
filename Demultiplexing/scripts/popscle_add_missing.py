#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--best", required=True, type=str, help="")
parser.add_argument("--barcodes", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import pandas as pd

print("Loading data ...")
df = pd.read_csv(args.best, sep="\t", header=0, index_col=None)
print("\tLoaded {} with shape: {}".format(os.path.basename(args.best), df.shape))

barcodes = pd.read_csv(args.barcodes, sep="\t", header=None, index_col=None)
print("\tLoaded {} with shape: {}".format(os.path.basename(args.barcodes), barcodes.shape))

if df.shape[0] > barcodes.shape[0]:
    print("Error, best file cannot have more rows than the barcodes file.")
    exit()

print("Adding missing barcodes ...")

# Save the counts to verify nothing went wrong later.
counts = df["DROPLET.TYPE"].value_counts()

# Drop the INT_ID column since it does not match the real alphabetical sorting order.
df.drop(["INT_ID"], axis=1, inplace=True)

# Reformat the barcodes in the alphabetical order that demuxlet uses. Here we add
# a new INT_ID column since the other one has missing barcodes.
barcodes = barcodes.sort_values(by=0).reset_index(drop=True)
barcodes = barcodes.reset_index(drop=False).rename(columns={"index": "INT_ID", 0: "BARCODE"})

# Check there are no duplicate barcodes.
if len(df["BARCODE"].unique()) != df.shape[0]:
    print("Error, best file contains duplicate barcodes.")
    exit()
if len(barcodes["BARCODE"].unique()) != barcodes.shape[0]:
    print("Error, barcodes file contains duplicate barcodes.")
    exit()

# Merge them together and reorder.
pre_shape = df.shape[0]
df = barcodes.merge(df, how="left").sort_values(by="INT_ID")
print("\tAdded {:,} empty barcode(s)".format(df.shape[0] - pre_shape))
del pre_shape, barcodes

# Check the counts to verify nothing went wrong.
new_counts = df["DROPLET.TYPE"].value_counts()
if not counts.equals(new_counts):
    print("Error, something went wrong with merging the missing barcodes.")
    exit()

# Fill in the missing droplet types with AMB.
df["DROPLET.TYPE"] = df["DROPLET.TYPE"].fillna("AMB")

print("Saving output ...")
df.to_csv(args.out, sep="\t", header=True, index=False)
print("\tSaved {} with shape: {}".format(os.path.basename(args.out), df.shape))
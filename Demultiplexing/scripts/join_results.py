#!/usr/bin/env python
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="")
# parser.add_argument("--demuxlet", required = True, help = "")
# parser.add_argument("--souporcell", required = True, help = "")
parser.add_argument("--scrublet", required = True, help = "")
parser.add_argument("--scds", required = True, help = "")
parser.add_argument("--doubletdetection", required = True, help = "")
parser.add_argument("--output", required = True, help = "")
args = parser.parse_args()

# demuxlet = pd.read_csv(args.demuxlet, sep="\t", header=0, index_col=0)
# souporcell = pd.read_csv(args.souporcell, sep="\t", header=0, index_col=0)
scrublet = pd.read_csv(args.scrublet, sep="\t", header=0, index_col=0)
scds = pd.read_csv(args.scds, sep="\t", header=0, index_col=0)
doubletdetection = pd.read_csv(args.doubletdetection, sep="\t", header=0, index_col=0)

# if demuxlet.shape[0] != souporcell.shape[0] != scrublet.shape[0] != scds.shape[0] != doubletdetection.shape[0]:
if scrublet.shape[0] != scds.shape[0] != doubletdetection.shape[0]:
    print("Error, not identical input sizes.")
    exit()

# df = demuxlet.merge(souporcell, left_index=True, right_index=True).merge(scrublet, left_index=True, right_index=True).merge(scds, left_index=True, right_index=True).merge(doubletdetection, left_index=True, right_index=True)
df = scrublet.merge(scds, left_index=True, right_index=True).merge(doubletdetection, left_index=True, right_index=True)

# if df.shape[0] != demuxlet.shape[0]:
if df.shape[0] != scrublet.shape[0]:
    print("Error, not identical input indices.")

df.to_csv(args.output, sep = "\t", index = True, header=True)
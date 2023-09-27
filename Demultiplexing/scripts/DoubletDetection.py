#!/usr/bin/env python
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/DoubletDetection.py
import argparse
import os


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-c", "--counts", required = True, help = "Path to the 10x filtered h5 file.")
parser.add_argument("-b", "--barcodes", required = True, help = "Path to the 10x droplet barcodes.")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, type = int, help = "Number of iterations to use; default is 50.")
parser.add_argument("-p", "--phenograph", required = False, default = False, help = "Whether to use phenograph (True) or not (False); default is False.")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True.")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, type = float, help = "P-value threshold for doublet calling; default is 1e-16.")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.5, type = float, help = "Voter threshold for doublet calling; default is 0.5.")
parser.add_argument("-o", "--out", required = True, help = "The output directory where results will be saved.")
args = parser.parse_args()

import numpy as np
import doubletdetection
import matplotlib
matplotlib.use('PDF')
import pandas as pd
import scanpy

if args.phenograph == 'True':
    pheno = True
elif args.phenograph == 'False':
    pheno = False
else:
    pheno = args.phenograph

if args.standard_scaling == 'True':
    standard_scaling = True
elif args.standard_scaling == 'False':
    standard_scaling = False
else:
    standard_scaling = args.standard_scaling

if not os.path.isdir(args.out):
    os.mkdir(args.out)

# Read in data
raw_counts = scanpy.read_10x_h5(args.counts)
barcodes_df = pd.read_csv(args.barcodes, sep="\t", header=None, names=["Barcode"])

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, use_phenograph=pheno, standard_scaling=standard_scaling, verbose=True)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes_df, results], axis=1)
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

print("Writing results to {}.".format(os.path.join(args.out, 'DoubletDetection_doublets_singlets.tsv')))
dataframe.to_csv(os.path.join(args.out, 'DoubletDetection_doublets_singlets.tsv'), sep = "\t", index = False)

# Figures
doubletdetection.plot.convergence(clf, save=os.path.join(args.out, 'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.out, 'threshold_test.pdf'), show=False, p_step=6)

# Make summary of singlets and doublets and write to file
summary = pd.DataFrame(dataframe.DoubletDetection_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'DoubletDetection_DropletType': 'Droplet N'}, axis=1)

print("Writing summary to {}.".format(os.path.join(args.out, 'DoubletDetection_summary.tsv')))
summary.to_csv(os.path.join(args.out, 'DoubletDetection_summary.tsv'), sep="\t", index=False)


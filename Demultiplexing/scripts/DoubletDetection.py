#!/usr/bin/env python
import numpy as np
import doubletdetection
import tarfile
import matplotlib.pyplot as plt
import os
import argparse
import sys
import pandas as pd


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix.mtx")
parser.add_argument("-b", "--barcodes", required = True, help = "cell ranger barcodes.tsv or barcodes.tsv.gz")
parser.add_argument("-d", "--mods_dir", required = True, help = "The mods directory")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("-p", "--phenograph", required = False, default = False, help = "Whether to use phenograph (True) or not (False); default is False")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.5, type = float, help = "Voter threshold for doublet calling; default is 0.5")
args = parser.parse_args()

print(args.p_thresh)
print(args.voter_thresh)
print(args.phenograph)
if args.phenograph == 'True':
    pheno = True
elif args.phenograph == 'False':
    pheno = False
else:
    pheno = args.phenograph
print(pheno)

if args.standard_scaling == 'True':
    standard_scaling = True
elif args.standard_scaling == 'False':
    standard_scaling = False
else:
    standard_scaling = args.standard_scaling
print(standard_scaling)

sys.path.append(args.mods_dir) 
import read10x

### Read in data ###
raw_counts = read10x.import_cellranger_mtx(args.counts_matrix)
barcodes_df = read10x.read_barcodes(args.barcodes)

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, use_phenograph=pheno, standard_scaling=standard_scaling, verbose = True)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes_df, results], axis=1)
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

dataframe.to_csv(os.path.join(args.outdir,'DoubletDetection_results.txt'), sep = "\t", index = False)


### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=True, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

f2, umap_coords = doubletdetection.plot.umap_plot(raw_counts, doublets, random_state=1, 
                                                       save=os.path.join(args.outdir,'umap_test.pdf'), show=True)
f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=True, p_step=6)
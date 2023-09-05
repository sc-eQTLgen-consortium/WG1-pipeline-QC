#!/usr/bin/env python
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/Scrublet.py
import argparse
import sys
import os
import scanpy

parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-c", "--counts", required = True, help = "Path to the 10x filtered h5 file.")
parser.add_argument("-b", "--barcodes", required=False, help="barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("-f", "--filtered_barcodes", required=False, default=None, help="File containing a filtered list of droplet barcodes. This may be used if you want to use a filtered list of barcodes for doublet detection (ie need to remove droplets that are empty or high in ambient RNA).")
parser.add_argument("-r", "--sim_doublet_ratio", required=False, default=2, type=int, help="Number of doublets to simulate relative to the number of observed transcriptomes.")
parser.add_argument("-c", "--min_counts", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("-e", "--min_cells", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than are excluded.")
parser.add_argument("-v", "--min_gene_variability_pctl", required=False, default=85, type=int, help="Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].")
parser.add_argument("-p", "--n_prin_comps", required=False, default=30, type=int, help="Number of principal components used to embed the transcriptomes priorto k-nearest-neighbor graph construction.")
parser.add_argument("-t", "--scrublet_doublet_threshold", required=False, default=None, type=float, help="Manually set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time")
parser.add_argument("-o", "--out", required = True, help = "The output directory where results will be saved.")
args = parser.parse_args()

import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
import numba
import numba.typed

# Get path of mods directory from current script directory
mods_path = "/opt/WG1-pipeline-QC/Demultiplexing/mods"
sys.path.append(mods_path)
import read10x

if not os.path.isdir(args.out):
    os.mkdir(args.out)


plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

### Read in data ###
counts_matrix = scanpy.read_10x_h5(args.counts)
barcodes_df = read10x.read_barcodes(args.barcodes)


dbl_rate = counts_matrix.shape[0]/1000 * 0.008
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio=args.sim_doublet_ratio)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts,
                                                          min_cells=args.min_cells,
                                                          min_gene_variability_pctl=args.min_gene_variability_pctl,
                                                          n_prin_comps=args.n_prin_comps)

### If running for the first time, the threshold will not manaually set so should continue but if manually set for a following run, need to assign it here


if args.scrublet_doublet_threshold is not None:
  print(args.scrublet_doublet_threshold)
  scrub.call_doublets(threshold=args.scrublet_doublet_threshold)

### Plotting and saving
scrub.plot_histogram()
plt.savefig(os.path.join(args.out, 'doublet_score_histogram.png'))
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig(os.path.join(args.out, 'UMAP.png'))

results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([barcodes_df, results, scores], axis=1)
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

print("Writing results to {}.".format(os.path.join(args.out, 'scrublet_doublets_singlets.tsv')))
dataframe.to_csv(os.path.join(args.out, 'scrublet_doublets_singlets.tsv'), sep="\t", index=False)

### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

print("Writing summary to {}.".format(os.path.join(args.out, 'scrublet_summary.tsv')))
summary.to_csv(os.path.join(args.out, 'scrublet_summary.tsv'), sep="\t", index=False)
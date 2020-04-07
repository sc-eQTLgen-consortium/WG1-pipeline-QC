#!/usr/bin/env python

import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix.mtx")
parser.add_argument("-g", "--genes", required = True, help = "genes.tsv from cellranger")
parser.add_argument("-d", "--expected_doublet_rate", required = False, default = 0.1, type = float, help = "expected doublet rate")
parser.add_argument("-r", "--sim_doublet_ratio", required = False, default = 2, type = int, help = "Number of doublets to simulate relative to the number of observed transcriptomes.")
parser.add_argument("-c", "--min_counts", required = False, default = 3, type = int, help = "Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("-e", "--min_cells", required = False, default = 3, type = int, help = "Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than are excluded.")
parser.add_argument("-v", "--min_gene_variability_pctl", required = False, default = 85, type = int, help = "Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].")
parser.add_argument("-p", "--n_prin_comps", required = False, default = 30, type = int, help = "Number of principal components used to embed the transcriptomes priorto k-nearest-neighbor graph construction.")
parser.add_argument("-t", "--scrublet_doublet_threshold", required = False, default = None, type = float, help = "Manually Set the scrublet doublet threshold location. For running a second time if scrublet incorreclty places the threshold the first time")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory")
args = parser.parse_args()


import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
print('scrublet' in sys.modules)

plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
counts_matrix = scipy.io.mmread(args.counts_matrix).T.tocsc()
genes = np.array(scr.load_genes(args.genes, delimiter='\t', column=1)) 
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate, sim_doublet_ratio = args.sim_doublet_ratio)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_cells, 
                                                          min_cells=args.min_cells, 
                                                          min_gene_variability_pctl=args.min_gene_variability_pctl, 
                                                          n_prin_comps=args.n_prin_comps)

### If running for the first time, the threshold will not manaually set so should continue but if manually set for a following run, need to assign it here

if args.scrublet_doublet_threshold is None:
  ### Plotting and saving
  scrub.plot_histogram();
  plt.savefig(os.path.join(args.outdir,'doublet_score_histogram.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(args.outdir,'UMAP.png'))

  np.savetxt(os.path.join(args.outdir,'predicted_doublet_mask.txt'), scrub.predicted_doublets_, fmt='%s')
  np.savetxt(os.path.join(args.outdir,'doublet_scores.txt'), scrub.doublet_scores_obs_, fmt='%.4f')   

else:
  print(args.scrublet_doublet_threshold)
  scrub.call_doublets(threshold=args.scrublet_doublet_threshold)
  ### Plotting and saving
  scrub.plot_histogram(); 
  plt.savefig(os.path.join(args.outdir,'doublet_score_histogram_manual_threshold.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(args.outdir,'UMAP_manual_threshold.png'))

  np.savetxt(os.path.join(args.outdir,'predicted_doublet_mask_manual_threshold.txt'), scrub.predicted_doublets_, fmt='%s')
  np.savetxt(os.path.join(args.outdir,'doublet_scores_manual_threshold.txt'), scrub.doublet_scores_obs_, fmt='%.4f')   
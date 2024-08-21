#!/usr/bin/env python
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/DoubletDetection.py
import argparse
import json
import os


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--counts", required=True, help="Path to the 10x filtered h5 file.")
parser.add_argument("--barcodes", required=True, help="Path to the 10x droplet barcodes.")
parser.add_argument("--boost_rate", required=False, default=0.25, type=float, help="Proportion of cell population size to produce as synthetic doublets.")
parser.add_argument("--n_components", required=False, default=30, type=int, help="Number of principal components used for clustering.")
parser.add_argument("--n_top_var_genes", required=False, default=10000, type=int, help="Number of highest variance genes to use; other genes discarded. Will use all genes when zero.")
parser.add_argument("--replace", action="store_true", required=False, default=False, help="If False, a cell will be selected as a synthetic doublet's parent no more than once. Defaults to False.")
parser.add_argument("--clustering_algorithm", choices=["louvain", "leiden", "phenograph"], default="phenograph", type=str, help="")
parser.add_argument("--n_iters", required=False, default=10, type=int, help="Number of fit operations from which to collect p-values.")
parser.add_argument("--pseudocount", required=False, default=0.1, type=float, help="Pseudocount used in normalize_counts. If `1` is used, and `standard_scaling=False`, the classifier is much more memory efficient; however, this may result in fewer doublets detected.")
parser.add_argument("--standard_scaling", action="store_true", required=False, default=False, help="Set to True to enable standard scaling of normalized count matrix prior to PCA. Recommended when not using Phenograph. Defaults to False.")
parser.add_argument("--n_jobs", required=False, default=1, type=int, help="Number of jobs to use. Speeds up neighbor computation.")
parser.add_argument("--p_thresh", required=False, default=1e-7, type=float, help="Hypergeometric test p-value threshold that determines per iteration doublet calls")
parser.add_argument("--voter_thresh", required=False, default=0.9, type=float, help="Fraction of iterations a cell must be called a doublet")
parser.add_argument("--out", required=True, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

with open(os.path.join(args.out, 'DoubletDetection_settings.json'), 'w') as f:
    json.dump(arguments, f, indent=4)
f.close()

import numpy as np
import doubletdetection
import matplotlib
matplotlib.use('PDF')
import pandas as pd
import scanpy

# Read in data
raw_counts = scanpy.read_10x_h5(args.counts).X
barcodes_df = pd.read_csv(args.barcodes, sep="\t", header=None, names=["Barcode"])

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Defaults: https://github.com/JonathanShor/DoubletDetection/blob/master/doubletdetection/doubletdetection.py
# boost_rate=0.25
# n_components=30
# n_top_var_genes=10000
# replace=False
# clustering_algorithm="phenograph"
# clustering_kwargs=None
# n_iters=10
# normalizer=None
# pseudocount=0.1
# random_state=0
# verbose=False
# standard_scaling=False
# n_jobs=1
clf = doubletdetection.BoostClassifier(
    boost_rate=args.boost_rate,
    n_components=args.n_components,
    n_top_var_genes=args.n_top_var_genes,
    replace=args.replace,
    clustering_algorithm=args.clustering_algorithm,
    clustering_kwargs=None,
    n_iters=args.n_iters,
    normalizer=None,
    pseudocount=args.pseudocount,
    random_state=0,
    verbose=True,
    standard_scaling=args.standard_scaling,
    n_jobs=args.n_jobs
)
# Defaults: https://github.com/JonathanShor/DoubletDetection/blob/master/doubletdetection/doubletdetection.py
# p_thresh=1e-7
# voter_thresh=0.9
labels = clf.fit(raw_counts).predict(
    p_thresh=args.p_thresh,
    voter_thresh=args.voter_thresh
)

# Figures
doubletdetection.plot.convergence(clf, save=os.path.join(args.out, 'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
doubletdetection.plot.threshold(clf, save=os.path.join(args.out, 'threshold_test.pdf'), show=False, p_step=6)

# Attributes:
# self.all_scores_ = np.zeros((self.n_iters, self._num_cells))
# self.all_log_p_values_ = np.zeros((self.n_iters, self._num_cells)) > required for plotting
# self.communities_ = np.zeros((self.n_iters, self._num_cells))
# labels_ = np.zeros(self._num_cells)
# suggested_score_cutoff_ = float, might not exist
# synth_communities_
# top_var_genes_
# voting_average_ = np.zeros(self._num_cells)


def save_df(df, fpath, name, sep="\t", index=False, header=True):
    compression = 'infer'
    if fpath.endswith('.gz'):
        compression = 'gzip'
    print("Writing {} to {}.".format(name, fpath))
    df.to_csv(fpath, sep=sep, header=header, index=index, compression=compression)

dataframe = pd.DataFrame({
    "Barcode": barcodes_df.iloc[:, 0].to_numpy(),
    "DoubletDetection_Labels": labels
})
dataframe["DoubletDetection_DropletType"] = dataframe["DoubletDetection_Labels"].map({1.0: "doublet", 0.0: "singlet"})
save_df(
    df=dataframe,
    fpath=os.path.join(args.out, 'DoubletDetection_doublets_singlets.tsv.gz'),
    name="results"
)

# Make summary of singlets and doublets and write to file
summary = pd.DataFrame(dataframe.DoubletDetection_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'DoubletDetection_DropletType': 'Droplet N'}, axis=1)
save_df(
    df=summary,
    fpath=os.path.join(args.out, 'DoubletDetection_summary.tsv.gz'),
    name="summary"
)

# Save the simulated info for plotting later.
save_df(
    df=pd.DataFrame(clf.all_log_p_values_),
    fpath=os.path.join(args.out, 'DoubletDetection_log_p_values.tsv.gz'),
    name="log p-values",
    header=False
)

#!/usr/bin/env python
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/Scrublet.py
import argparse
import os
import json

parser = argparse.ArgumentParser(
    description="wrapper for scrublet for doublet detection of transcriptomic data.")
parser.add_argument("--counts", required=True, help="Path to the 10x filtered h5 file.")
parser.add_argument("--barcodes", required=True, help="barcodes.tsv or barcodes.tsv.gz from cellranger")
parser.add_argument("--sim_doublet_ratio", required=False, default=2.0, type=float, help="Number of doublets to simulate relative to the number of observed transcriptomes.")
parser.add_argument("--n_neighbors", required=False, default=None, type=int, help="Number of neighbors used to construct the KNN graph of observed transcriptomes and simulated doublets. If `None`, this is  set to round(0.5 * sqrt(n_cells))")
parser.add_argument("--expected_doublet_scaling_factor", required=False, default=8e-06, type=float, help="The fraction of droublets expected based on the number of nuclei recovered.")
parser.add_argument("--stdev_doublet_rate", required=False, default=0.02, type=float, help="Uncertainty in the expected doublet rate.")

parser.add_argument("--synthetic_doublet_umi_subsampling", required=False, default=1.0, type=float, help="Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMIs from two randomly  sampled observed transcriptomes. For values less than 1, the  UMI counts are added and then randomly sampled at the specified rate.")
parser.add_argument("--get_doublet_neighbor_parents", action="store_true", required=False, default=False, help="If True, return the parent transcriptomes that generated the  doublet neighbors of each observed transcriptome. This information can  be used to infer the cell states that generated a given  doublet state.")
parser.add_argument("--min_counts", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than min_cells are excluded.")
parser.add_argument("--min_cells", required=False, default=3, type=int, help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_counts in fewer than are excluded.")
parser.add_argument("--min_gene_variability_pctl", required=False, default=85, type=int, help="Used for gene filtering prior to PCA. Keep the most highly variable genes in the top min_gene_variability_pctl percentile), as measured by the v-statistic [Klein et al., Cell 2015].")
parser.add_argument("--log_transform", action="store_true", default=False, help="If True, log-transform the counts matrix (log10(1+TPM)). `sklearn.decomposition.TruncatedSVD` will be used for dimensionality reduction, unless `mean_center` is True.")
parser.add_argument("--no_mean_center", dest="mean_center", action="store_false", default=True, help="If True, center the data such that each gene has a mean of 0. `sklearn.decomposition.PCA` will be used for dimensionality reduction.")
parser.add_argument("--no_normalize_variance", dest="normalize_variance", action="store_false", default=True, help="If True, normalize the data such that each gene has a variance of 1. `sklearn.decomposition.TruncatedSVD` will be used for dimensionality reduction, unless `mean_center` is True.")
parser.add_argument("--n_prin_comps", required=False, default=30, type=int, help="Number of principal components used to embed the transcriptomes priorto k-nearest-neighbor graph construction.")
parser.add_argument("--doublet_threshold", required=False, default=None, type=float, help="Manually set the scrublet doublet threshold location. For running a second time if scrublet incorrectly places the threshold the first time.")
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

with open(os.path.join(args.out, 'Scrublet_settings.json'), 'w') as f:
    json.dump(arguments, f, indent=4)
f.close()

import scanpy
import scrublet as scr
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pandas as pd

### Read in data ###
counts_matrix = scanpy.read_10x_h5(args.counts).X
barcodes_df = pd.read_csv(args.barcodes, sep="\t", header=None, names=["Barcode"])


print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
# Defaults: https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py
# total_counts=None
# sim_doublet_ratio=2.0
# n_neighbors=None
# expected_doublet_rate=0.1
# stdev_doublet_rate=0.02
# random_state=0
scrub = scr.Scrublet(counts_matrix,
                     total_counts=None,
                     sim_doublet_ratio=args.sim_doublet_ratio,
                     n_neighbors=args.n_neighbors,
                     expected_doublet_rate=counts_matrix.shape[0] * args.expected_doublet_scaling_factor,
                     stdev_doublet_rate=args.stdev_doublet_rate,
                     random_state=0)
# Defaults: https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py
# synthetic_doublet_umi_subsampling=1.0
# use_approx_neighbors=True
# distance_metric='euclidean'
# get_doublet_neighbor_parents=False
# min_counts=3
# min_cells=3
# min_gene_variability_pctl=85
# log_transform=False
# mean_center=True
# normalize_variance=True
# n_prin_comps=30
# svd_solver='arpack'
# verbose=True
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    synthetic_doublet_umi_subsampling=args.synthetic_doublet_umi_subsampling,
    use_approx_neighbors=True,
    distance_metric='euclidean',
    get_doublet_neighbor_parents=args.get_doublet_neighbor_parents,
    min_counts=args.min_counts,
    min_cells=args.min_cells,
    min_gene_variability_pctl=args.min_gene_variability_pctl,
    log_transform=args.log_transform,
    mean_center=args.mean_center,
    normalize_variance=args.normalize_variance,
    n_prin_comps=args.n_prin_comps,
    svd_solver='arpack',
    verbose=True
)

### If running for the first time, the threshold will not manaually set so should continue but if manually set for a following run, need to assign it here


if args.doublet_threshold is not None:
  print(args.doublet_threshold)
  scrub.call_doublets(threshold=args.doublet_threshold)

### Plotting and saving
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
scrub.plot_histogram()
plt.savefig(os.path.join(args.out, 'doublet_score_histogram.png'))
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig(os.path.join(args.out, 'UMAP.png'))

def save_df(df, fpath, name, sep="\t", index=False, header=True):
    compression = 'infer'
    if fpath.endswith('.gz'):
        compression = 'gzip'
    print("Writing {} to {}.".format(name, fpath))
    df.to_csv(fpath, sep=sep, header=header, index=index, compression=compression)

obs_dataframe = pd.DataFrame({
    "Barcode": barcodes_df.iloc[:, 0].to_numpy(),
    "scrublet_PredictedDoublets": scrub.predicted_doublets_,
    "scrublet_DoubletScores": scrub.doublet_scores_obs_,
    "scrublet_ZScores": scrub.z_scores_
})
obs_dataframe["scrublet_DropletType"] = obs_dataframe["scrublet_PredictedDoublets"].map({True: "doublet", False: "singlet"})
save_df(
    df=obs_dataframe,
    fpath=os.path.join(args.out, 'Scrublet_doublets_singlets.tsv.gz'),
    name="observed results"
)

### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(obs_dataframe.scrublet_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)
save_df(
    df=summary,
    fpath=os.path.join(args.out, 'Scrublet_summary.tsv.gz'),
    name="summary"
)

# Save the simulated info for plotting later.
sim_dataframe = pd.DataFrame({
    "scrublet_DoubletScores": scrub.doublet_scores_sim_,
    "scrublet_DoubletErrors": scrub.doublet_errors_sim_,
})
save_df(
    df=sim_dataframe,
    fpath=os.path.join(args.out, 'Scrublet_doublets_singlets_sim.tsv.gz'),
    name="simulated results"
)

# Save the observed UMAP info for plotting later.
save_df(
    df=pd.DataFrame(scrub.manifold_obs_),
    fpath=os.path.join(args.out, 'Scrublet_manifold.tsv.gz'),
    name="observed UMAP",
    header=False
)
save_df(
    df=pd.DataFrame(scrub.manifold_sim_),
    fpath=os.path.join(args.out, 'Scrublet_manifold_sim.tsv.gz'),
    name="simulated UMAP",
    header=False
)

print("Writing stats to {}.".format(os.path.join(args.out, 'Scrublet_stats.json')))
with open(os.path.join(args.out, 'Scrublet_stats.json'), 'w') as f:
    json.dump({
        "threshold": scrub.threshold_,
        "detected_doublet_rate": scrub.detected_doublet_rate_,
        "detectable_doublet_fraction": scrub.detectable_doublet_fraction_,
        "overall_doublet_rate": scrub.overall_doublet_rate_
    }, f, indent=4)
f.close()


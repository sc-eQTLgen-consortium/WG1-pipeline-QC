#!/usr/bin/env python
# Adapted from https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py
# and https://github.com/AllonKleinLab/SPRING_dev/blob/master/data_prep/spring_helper.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--settings", required=True, nargs="+", type=str, help="")
parser.add_argument("--results", required=True, nargs="+", type=str, help="")
parser.add_argument("--stats", required=True, nargs="+", type=str, help="")
parser.add_argument("--manifolds", required=True, nargs="+", type=str, help="")
parser.add_argument("--n_jobs", required=False, default=-1, type=int, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import umap

"""
singularity exec --bind /groups/umcg-biogen/tmp02/,/groups/umcg-biogen/tmp02/users/umcg-mvochteloo/simulated_home:/home/umcg-mvochteloo /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/20231103-2-WG1-pipeline-QC.sif python /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-27-DemultiplexingAndDoubletRemoval/scripts/plot_Scrublet.py             --settings /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun1/Scrublet_settings.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun2/Scrublet_settings.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun3/Scrublet_settings.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun4/Scrublet_settings.json             --results /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun1/Scrublet_doublets_singlets.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun2/Scrublet_doublets_singlets.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun3/Scrublet_doublets_singlets.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun4/Scrublet_doublets_singlets.tsv.gz             --manifolds /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun1/Scrublet_manifold.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun2/Scrublet_manifold.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun3/Scrublet_manifold.tsv.gz /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun4/Scrublet_manifold.tsv.gz             --stats /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun1/Scrublet_stats.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun2/Scrublet_stats.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun3/Scrublet_stats.json /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/test_dataset/ScrubletRun4/Scrublet_stats.json             --out /groups/umcg-biogen/tmp02/output/2022-09-01-scMetaBrainConsortium/2023-10-16-scMetaBrain-WorkGroup1QC/2023-10-31-ImputationTestDataset/Step2-DemultiplexingAndDoubletRemoval/figures/test_dataset/ --pool test_dataset
"""

setting_abbr = {
    "sim_doublet_ratio": "sim_dbl_ratio",
    "n_neighbors": "n_neighbors",
    "expected_doublet_scaling_factor": "dbl_sf",
    "stdev_doublet_rate": "stdev_dbl_rate",
    "synthetic_doublet_umi_subsampling": "dbl_umi_sub",
    "get_doublet_neighbor_parents": "dbl_neighbor",
    "min_counts": "min_counts",
    "min_cells": "min_cells",
    "min_gene_variability_pctl": "min_gene_var",
    "log_transform": "log_trans",
    "mean_center": "mean_center",
    "normalize_variance": "norm_var",
    "n_prin_comps": "n_prin_comps",
    "doublet_threshold": "dblt_thresh"
}


def plot_settings(ax, settings):
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame({setting_abbr[key]: value for key, value in settings.items() if key in setting_abbr}, index=[0]).T
    max_key_length = max([len(str(key)) for key in df.index])
    max_value_length = max([len(str(value)) for value in df[0]])
    total_length = max_key_length + max_value_length
    table = ax.table(cellText=df.values, colWidths=[total_length / max_key_length, total_length / max_value_length], rowLabels=df.index, loc='center', edges='open')
    table.auto_set_column_width(col=list(range(len(df.columns))))
    table.scale(1, 1.2)


def plot_histogram(ax, y, threshold, scale_hist='log', title=''):
    ax.hist(y, np.linspace(0, 1, 50), color='gray', linewidth=0, density=True)
    ax.set_yscale(scale_hist)
    yl = ax.get_ylim()
    ax.set_ylim(yl)
    ax.plot(threshold * np.ones(2), yl, c='black', linewidth=1)
    ax.set_title(title)
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Prob. density')


def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii, 0] = curcol[0] * .9
        cdat[ii, 1] = curcol[1] * .9
        cdat[ii, 2] = curcol[2] * .9
        cdat[ii, 3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap


def custom_cmap(rgb_list):
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0],rgb_list)
    return cmap


def plot_embedding(axs, embedding, embedding_name, color_dat, called_doubs):
    x = embedding[:, 0]
    y = embedding[:, 1]
    xl = (x.min() - x.ptp() * .05, x.max() + x.ptp() * 0.05)
    yl = (y.min() - y.ptp() * .05, y.max() + y.ptp() * 0.05)

    ax = axs[1]
    vmin = color_dat.min()
    vmax = color_dat.max()
    cmap_use = darken_cmap(plt.cm.Reds, 0.9)
    o = np.argsort(color_dat)
    pp = ax.scatter(x[o], y[o], s=5, edgecolors=None, c=color_dat[o],
                    cmap=cmap_use, vmin=vmin, vmax=vmax)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Doublet score')
    ax.set_xlabel(embedding_name + ' 1')
    ax.set_ylabel(embedding_name + ' 2')
    fig.colorbar(pp, ax=ax)

    ax = axs[0]
    ax.scatter(x[o], y[o], s=5, edgecolors=None, c=called_doubs[o], cmap=custom_cmap([[.7, .7, .7], [0, 0, 0]]))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Predicted doublets')
    ax.set_xlabel(embedding_name + ' 1')
    ax.set_ylabel(embedding_name + ' 2')

################################################################################

nrows = len(args.settings)

print("Plotting")
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(nrows, 5, figsize=(20, 4 * nrows), gridspec_kw={"width_ratios": [0.1, 0.225, 0.225, 0.225, 0.225]})
if nrows == 1:
    axs = axs[np.newaxis, ...]

for row_index, (settings_path, results_path, stats_path, manifold_path) in enumerate(list(zip(args.settings, args.results, args.stats, args.manifolds))):
    print("\tRow {}:".format(row_index))
    print("\t  --settings {}".format(settings_path))
    print("\t  --result {}".format(results_path))
    print("\t  --stats {}".format(stats_path))
    print("\t  --manifold {}".format(manifold_path))
    print("")

    fh = open(settings_path)
    settings = json.load(fh)
    fh.close()
    plot_settings(
        ax=axs[row_index, 0],
        settings=settings
    )

    fh = open(stats_path)
    stats = json.load(fh)
    fh.close()

    results_df = pd.read_csv(results_path, sep="\t", header=0, index_col=None)
    plot_histogram(
        ax=axs[row_index, 1],
        y=results_df["scrublet_DoubletScores"].to_numpy(),
        threshold=stats["threshold"],
        scale_hist='log',
        title='Observed transcriptomes'
    )

    results_sim_df = pd.read_csv(results_path, sep="\t", header=0, index_col=None)
    plot_histogram(
        ax=axs[row_index, 2],
        y=results_sim_df["scrublet_DoubletScores"].to_numpy(),
        threshold=stats["threshold"],
        scale_hist='linear',
        title='Simulated doublets'
    )

    manifold_df = pd.read_csv(manifold_path, sep="\t", header=None, index_col=None)
    umap_embedding = umap.UMAP(n_neighbors=10, min_dist=0.3, metric='euclidean', random_state=0, n_jobs=args.n_jobs).fit_transform(manifold_df.to_numpy())

    results_df = pd.read_csv(results_path, sep="\t", header=0, index_col=None)

    plot_embedding(
        axs=axs[row_index, 3:],
        embedding=umap_embedding,
        embedding_name='UMAP',
        color_dat=results_df["scrublet_DoubletScores"],
        called_doubs=results_df["scrublet_PredictedDoublets"],
    )

fig.suptitle(args.pool)
fig.tight_layout()
plt.savefig(os.path.join(args.out, 'Scrublet_histograms_and_UMAPs.png'), bbox_inches="tight")
print("")

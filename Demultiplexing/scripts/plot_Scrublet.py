#!/usr/bin/env python
# Adapted from https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py
# and https://github.com/AllonKleinLab/SPRING_dev/blob/master/data_prep/spring_helper.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--infolers", required=True, nargs="+", type=str, help="")
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


settings_info = {
    "sim_doublet_ratio": ("sim_dbl_ratio", 2.0),
    "n_neighbors": ("n_neighbors", None),
    "expected_doublet_scaling_factor": ("dbl_sf", 8e-06),
    "stdev_doublet_rate": ("stdev_dbl_rate", 0.02),
    "synthetic_doublet_umi_subsampling": ("dbl_umi_sub", 1.0),
    "get_doublet_neighbor_parents": ("dbl_neighbor", False),
    "min_counts": ("min_counts", 3),
    "min_cells": ("min_cells", 3),
    "min_gene_variability_pctl": ("min_gene_var", 85),
    "log_transform": ("log_trans", False),
    "mean_center": ("mean_center", True),
    "normalize_variance": ("norm_var", True),
    "n_prin_comps": ("n_prin_comps", 30),
    "doublet_threshold": ("dblt_thresh", None)
}


def plot_table(ax, values, colors=None, title=""):
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(values, index=[0]).T
    max_key_length = max([len(str(key)) for key in df.index])
    max_value_length = max([len(str(value)) for value in df[0]])
    total_length = max_key_length + max_value_length
    table = ax.table(cellText=df.values, colWidths=[total_length / max_key_length, total_length / max_value_length], rowLabels=df.index, loc='center', edges='open')
    table.auto_set_column_width(col=list(range(len(df.columns))))
    table.scale(1, 1.2)
    if colors is not None:
        for row_index in range(len(values)):
            key = table[(row_index, -1)].get_text()
            key.set_color(colors[key.get_text()])
            value = table[(row_index, 0)].get_text()
            value.set_color(colors[key.get_text()])
    ax.set_title(title)


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
        cdat[ii, 0] = curcol[0] * scale_factor
        cdat[ii, 1] = curcol[1] * scale_factor
        cdat[ii, 2] = curcol[2] * scale_factor
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

nrows = len(args.infolders)

print("Plotting")
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(nrows, 5, figsize=(20, 4 * nrows), gridspec_kw={"width_ratios": [0.1, 0.225, 0.225, 0.225, 0.225]})
if nrows == 1:
    axs = axs[np.newaxis, ...]

for row_index, (infolder) in enumerate(args.infolders):
    print("\tRow {}:".format(row_index))
    print("\t  --infolder {}".format(infolder))
    print("")

    settings_path = os.path.join(infolder, "Scrublet_settings.json")
    results_path = os.path.join(infolder, "Scrublet_doublets_singlets.tsv.gz")
    sim_results_path = os.path.join(infolder, "Scrublet_doublets_singlets_sim.tsv.gz")
    stats_path = os.path.join(infolder, "Scrublet_stats.json")
    manifold_path = os.path.join(infolder, "Scrublet_manifold.tsv.gz")

    fh = open(settings_path)
    settings = {}
    colors = {}
    for key, value in json.load(fh).items():
        if key in settings_info:
            new_key, default_value = settings_info[key]
            settings[new_key] = value
            if str(value) == str(default_value):
                colors[new_key] = "black"
            else:
                colors[new_key] = "red"
    fh.close()
    plot_table(
        ax=axs[row_index, 0],
        values=settings,
        colors=colors
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

    results_sim_df = pd.read_csv(sim_results_path, sep="\t", header=0, index_col=None)
    plot_histogram(
        ax=axs[row_index, 2],
        y=results_sim_df["scrublet_DoubletScores"].to_numpy(),
        threshold=stats["threshold"],
        scale_hist='linear',
        title='Simulated doublets'
    )

    manifold_df = pd.read_csv(manifold_path, sep="\t", header=None, index_col=None)
    umap_embedding = umap.UMAP(n_neighbors=10, min_dist=0.3, metric='euclidean', random_state=0, n_jobs=args.n_jobs).fit_transform(manifold_df.to_numpy())

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

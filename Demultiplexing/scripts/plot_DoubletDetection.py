#!/usr/bin/env python
# Adapted from https://github.com/JonathanShor/DoubletDetection/blob/master/doubletdetection/plot.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--infolders", required=True, nargs="+", type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import warnings

settings_info = {
    "boost_rate": ("boost_rate", 0.25),
    "n_components": ("n_comp", 30),
    "n_top_var_genes": ("n_var_genes", 10000),
    "replace": ("replace", False),
    "clustering_algorithm": ("clust_algo", "louvain"),
    "n_iters": ("n_iters", 50),
    "pseudocount": ("pseudocount", 0.1),
    "standard_scaling": ("stand_scaling", True),
    "p_thresh": ("p_thresh", 1e-16),
    "voter_thresh": ("voter_thresh", 0.5)
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
    table.scale(1, 1.1)
    if colors is not None:
        for row_index in range(len(values)):
            key = table[(row_index, -1)].get_text()
            key.set_color(colors[key.get_text()])
            value = table[(row_index, 0)].get_text()
            value.set_color(colors[key.get_text()])
    ax.set_title(title)


def convergence(ax, n_iters, all_log_p_values, p_thresh, voter_thresh):
    log_p_thresh = np.log(p_thresh)
    doubs_per_run = []
    # Ignore numpy complaining about np.nan comparisons
    with np.errstate(invalid="ignore"):
        for i in range(n_iters):
            cum_log_p_values = all_log_p_values[: i + 1]
            cum_vote_average = np.mean(
                np.ma.masked_invalid(cum_log_p_values) <= log_p_thresh, axis=0
            )
            cum_doublets = np.ma.filled((cum_vote_average >= voter_thresh).astype(float), np.nan)
            doubs_per_run.append(np.nansum(cum_doublets))

    ax.plot(np.arange(len(doubs_per_run)), doubs_per_run)
    ax.set_xlabel("Number of Iterations")
    ax.set_ylabel("Number of Predicted Doublets")
    ax.set_title("Predicted Doublets per Iteration")


def threshold(f, ax, all_log_p_values, log10=True, log_p_grid=None, voter_grid=None, v_step=2, p_step=5):
    # Ignore numpy complaining about np.nan comparisons
    with np.errstate(invalid="ignore"):
        all_log_p_values_ = np.copy(all_log_p_values)
        if log10:
            all_log_p_values_ /= np.log(10)
        if log_p_grid is None:
            log_p_grid = np.arange(-100, -1)
        if voter_grid is None:
            voter_grid = np.arange(0.3, 1.0, 0.05)
        doubs_per_t = np.zeros((len(voter_grid), len(log_p_grid)))
        for i in range(len(voter_grid)):
            for j in range(len(log_p_grid)):
                voting_average = np.mean(
                    np.ma.masked_invalid(all_log_p_values_) <= log_p_grid[j], axis=0
                )
                labels = np.ma.filled((voting_average >= voter_grid[i]).astype(float), np.nan)
                doubs_per_t[i, j] = np.nansum(labels)

    cax = ax.imshow(doubs_per_t, cmap="hot", aspect="auto")
    ax.set_xticks(np.arange(len(log_p_grid))[::p_step])
    ax.set_xticklabels(np.around(log_p_grid, 1)[::p_step], rotation="vertical")
    ax.set_yticks(np.arange(len(voter_grid))[::v_step])
    ax.set_yticklabels(np.around(voter_grid, 2)[::v_step])
    cbar = f.colorbar(cax)
    cbar.set_label("Predicted Doublets")
    if log10 is True:
        ax.set_xlabel("Log10 p-value")
    else:
        ax.set_xlabel("Log p-value")
    ax.set_ylabel("Voting Threshold")
    ax.set_title("Threshold Diagnostics")

################################################################################

nrows = len(args.infolders)

print("Plotting")
# Ignore warning for convergence plot
with warnings.catch_warnings():
    warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")
    f, axs = plt.subplots(nrows, 3, figsize=(12, 4 * nrows), dpi=150, gridspec_kw={"width_ratios": [0.1, 0.4, 0.4]})
    if nrows == 1:
        axs = axs[np.newaxis, ...]

    for row_index, (infolder) in enumerate(args.infolders):
        print("\tRow {}:".format(row_index))
        print("\t  --infolder {}".format(infolder))
        print("")

        settings_path = os.path.join(infolder, "DoubletDetection_settings.json")
        log_p_values_path = os.path.join(infolder, "DoubletDetection_log_p_values.tsv.gz")

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

        log_p_values_df = pd.read_csv(log_p_values_path, sep="\t", header=None, index_col=None)
        convergence(
            ax=axs[row_index, 1],
            n_iters=settings["n_iters"],
            all_log_p_values=log_p_values_df.to_numpy(),
            p_thresh=settings["p_thresh"],
            voter_thresh=settings["voter_thresh"])
        threshold(
            f=f,
            ax=axs[row_index, 2],
            all_log_p_values=log_p_values_df.to_numpy(),
            p_step=6
        )

    f.suptitle(args.pool)
    f.tight_layout()
    f.savefig(os.path.join(args.out, 'DoubletDetection_convergence_and_threshold_test.png'), format="png", bbox_inches="tight")

#!/usr/bin/env python
# Adapted from https://github.com/JonathanShor/DoubletDetection/blob/master/doubletdetection/plot.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--settings", required=True, nargs="+", type=str, help="")
parser.add_argument("--log_p_values", required=True, nargs="+", type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import warnings

setting_abbr = {
    "boost_rate": "boost_rate",
    "n_components": "n_comp",
    "n_top_var_genes": "n_var_genes",
    "replace": "replace",
    "clustering_algorithm": "clust_algo",
    "n_iters": "n_iters",
    "pseudocount": "pseudocount",
    "standard_scaling": "stand_scaling",
    "p_thresh": "p_thresh",
    "voter_thresh": "voter_thresh"
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
    table.scale(1, 1.1)


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

nrows = len(args.settings)

print("Plotting")
# Ignore warning for convergence plot
with warnings.catch_warnings():
    warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")
    f, axs = plt.subplots(nrows, 3, figsize=(12, 4 * nrows), dpi=150, gridspec_kw={"width_ratios": [0.1, 0.4, 0.4]})
    if nrows == 1:
        axs = axs[np.newaxis, ...]

    for row_index, (settings_path, log_p_values_path) in enumerate(list(zip(args.settings, args.log_p_values))):
        print("\tRow {}:".format(row_index))
        print("\t  --settings {}".format(settings_path))
        print("\t  --p_values {}".format(log_p_values_path))
        print("")

        fh = open(settings_path)
        settings = json.load(fh)
        fh.close()
        plot_settings(
            ax=axs[row_index, 0],
            settings=settings
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
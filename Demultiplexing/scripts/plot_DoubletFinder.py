#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--settings", required=True, nargs="+", type=str, help="")
parser.add_argument("--bcmvns", required=True, nargs="+", type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

setting_abbr = {
    "dims": "dims",
    "resolution": "resolution",
    "expected_doublet_scaling_factor": "dbl_sf",
    "pn": "pn"
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
    table.scale(1, 2)


def plot_scatterplot(ax, x, y):
    ax.scatter(x, y, color="black")
    ax.set_xlabel('pK')
    ax.set_ylabel('BCmetric')

################################################################################

nrows = len(args.settings)

print("Plotting")
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(nrows, 2, figsize=(8, 4 * nrows), gridspec_kw={"width_ratios": [0.1, 0.9]})
if nrows == 1:
    axs = axs[np.newaxis, ...]

for row_index, (settings_path, bcmvn_path) in enumerate(list(zip(args.settings, args.bcmvns))):
    print("\tRow {}:".format(row_index))
    print("\t  --settings {}".format(settings_path))
    print("\t  --bcmvn {}".format(bcmvn_path))
    print("")

    fh = open(settings_path)
    settings = json.load(fh)
    fh.close()
    plot_settings(
        ax=axs[row_index, 0],
        settings=settings
    )

    bcmvn_df = pd.read_csv(bcmvn_path, sep="\t", header=0, index_col=None)
    plot_scatterplot(
        ax=axs[row_index, 1],
        x=bcmvn_df["pK"].to_numpy(),
        y=bcmvn_df["BCmetric"].to_numpy()
    )

fig.suptitle(args.pool)
fig.tight_layout()
plt.savefig(os.path.join(args.out, 'DoubletFinder_pKvBCmetrics.png'), bbox_inches="tight")
print("")
#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--infolders", required=True, nargs="+", type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

settings_info = {
    "dims": ("dims", 10),
    "resolution": ("resolution", 0.1),
    "expected_doublet_scaling_factor": ("dbl_sf", 8e-06),
    "pn": ("pn", 0.25),
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
    table.scale(1, 2)
    if colors is not None:
        for row_index in range(len(values)):
            key = table[(row_index, -1)].get_text()
            key.set_color(colors[key.get_text()])
            value = table[(row_index, 0)].get_text()
            value.set_color(colors[key.get_text()])
    ax.set_title(title)


def plot_scatterplot(ax, x, y):
    ax.scatter(x, y, color="black")
    ax.set_xlabel('pK')
    ax.set_ylabel('BCmetric')

################################################################################

nrows = len(args.infolders)

print("Plotting")
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(nrows, 2, figsize=(8, 4 * nrows), gridspec_kw={"width_ratios": [0.1, 0.9]})
if nrows == 1:
    axs = axs[np.newaxis, ...]

for row_index, (infolder) in enumerate(args.infolders):
    print("\tRow {}:".format(row_index))
    print("\t  --infolder {}".format(infolder))
    print("")

    settings_path = os.path.join(infolder, "DoubletFinder_settings.json")
    bcmvn_path = os.path.join(infolder, "DoubletFinder_bcmvn.tsv.gz")

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
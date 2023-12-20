#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--assignments", required=True, type=str, help="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_barplot(ax, labels, values, color="black", ylabel="", title=""):
    ax.bar(labels, values, color=color)
    ax.set(ylabel=ylabel, title=title, ylim=(0, max(values) * 1.1))
    ax.set_xticks(range(len(labels)), labels, rotation=90)

################################################################################

print("Loading data")

assignments = pd.read_csv(args.assignments, sep="\t", dtype=str)
poolsheet = pd.read_csv(args.poolsheet, sep="\t", dtype=str)

data = {}
for column in assignments.columns:
    if column.endswith("_DropletType"):
        for _, row in poolsheet.iterrows():
            method = column.split("_")[0]
            results = assignments.loc[assignments["Pool"] == row["Pool"], column].value_counts().to_frame()
            results.reset_index(inplace=True)
            results.columns = ["index", row["Pool"]]
            results.set_index("index", inplace=True)
            results.index.name = None
            if method in data:
                data[method].append(results)
            else:
                data[method] = [results]

ncols = len(data)
npools = poolsheet.shape[0]

print("Plotting")
plt.rcParams['pdf.fonttype'] = 42
fig, axs = plt.subplots(3, ncols, figsize=(0.2 * npools * ncols, 12), sharey='row')
if ncols == 1:
    axs = axs[np.newaxis, ...]

for col_index, (method, dfs) in enumerate(data.items()):
    print("\tMethod: {}".format(method))
    df = pd.concat(dfs, axis=1)
    df.loc["%doublet", :] = (df.loc["doublet", :] / df.sum(axis=0)) * 100

    plot_barplot(
        ax=axs[0, col_index],
        labels=df.columns,
        values=df.loc["singlet", :].values,
        color="#E65A25",
        ylabel="number of singlets" if col_index == 0 else "",
        title=method
    )

    plot_barplot(
        ax=axs[1, col_index],
        labels=df.columns,
        values=df.loc["doublet", :].values,
        color="#00A174",
        ylabel="number of doublets" if col_index == 0 else ""
    )

    plot_barplot(
        ax=axs[2, col_index],
        labels=df.columns,
        values=df.loc["%doublet", :].values,
        color="#00A174",
        ylabel="% of doublets" if col_index == 0 else ""
    )

fig.suptitle("Number of doublet and singlets per pool")
fig.tight_layout()
plt.savefig(os.path.join(args.out, 'doublets_singlets.png'), bbox_inches="tight")
print("")

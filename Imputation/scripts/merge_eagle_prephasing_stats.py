#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import re
import os

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--logs", nargs="+", type=str, required=True, help="The input log files.")
parser.add_argument("--out", type=str, required=True, help="The output file.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

info = {}
samples = set()
chromosomes = set()
for log_file in args.logs:
    chr = os.path.basename(log_file).split("_")[-1].split(".")[0]
    chromosomes.add(chr)
    print("Loading phase_confidence from chr{}".format(chr))
    flag = False
    with open(log_file, "r") as f:
        for line in f:
            if line == "ID\tPHASE_CONFIDENCE\n":
                flag = True
                continue

            if flag:
                values = line.strip("\n").split("\t")

                if len(values) != 2:
                    flag = False
                    break

                if values[0] not in info:
                    info[values[0]] = {}
                    samples.add(values[0])
                info[values[0]][chr] = values[1]
    f.close()


def natural_keys(text):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]


samples = list(samples)
samples.sort(key=natural_keys)

chromosomes = list(chromosomes)
chromosomes.sort(key=natural_keys)

all_values = []

print("Calculating averages per sample:")
for sample in samples:
    sample_values = []
    for chr in chromosomes:
        sample_values.append(float(info[sample][chr]))
        all_values.append(float(info[sample][chr]))

    average = round(sum(sample_values) / len(sample_values), 6)
    print("\t{}: {}".format(sample, average))
    info[sample]["avg"] = str(average)
print("")

print("Calculating averages per chromosome:")
info["avg"] = {}
for chr in chromosomes:
    chr_values = []
    for sample in samples:
        chr_values.append(float(info[sample][chr]))
        all_values.append(float(info[sample][chr]))

    average = round(sum(chr_values) / len(chr_values), 6)
    print("\t{}: {}".format(chr, average))
    info["avg"][chr] = str(average)
print("")

info["avg"]["avg"] = str(round(sum(all_values) / len(all_values), 6))

chromosomes += ["avg"]
samples += ["avg"]

print("Saving stats to {}".format(args.out))
with open(args.out, "w") as f:
    f.write("\t".join(chromosomes) + "\n")
    for sample in samples:
        row = []
        for chr in chromosomes:
            value = ""
            if chr in info[sample]:
                value = info[sample][chr]
            else:
                print(sample, chr)
            row.append(value)
        f.write(sample + "\t" + "\t".join(row) + "\n")
f.close()

print("\tsaved matrix of shape {} rows by {} columns".format(len(samples), len(chromosomes)))
print("End")

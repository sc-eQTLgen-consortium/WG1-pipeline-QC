#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import numpy as np
import pandas as pd
import gzip
import os

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-i", "--input", nargs="*", type=str, dest='input_paths', required=True, help="The input VCF file.")
parser.add_argument("-o", "--output", type=str, dest='output_path', required=True, help="The output VCF file.")
parser.add_argument("-maf", "--minor_allele_frequency", type=float, dest='thresh_maf', default=0.01, help="The minor allele frequency threshold. Default: 0.01.")
parser.add_argument("-cr", "--call_rate", type=float, dest='thresh_cr', default=0.99, help="The call rate threshold. Default: 0.99.")
parser.add_argument("-hwe", "--hardy_weinberg_equilibrium", type=float, dest='thresh_hwe', default=1E-6, help="The hardy weinberg equilibrium threshold. Default: 1e-6.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

if not os.path.isdir(os.path.dirname(args.output_path)):
    os.mkdir(os.path.dirname(args.output_path))


def read_filter_logfile(filepath, thresh_maf, thresh_cr, thresh_hwe):
    variants = -1
    multi_allelic = 0
    indel_below_vqsr = 0
    indel_non_pass = 0
    snv_below_vqsr = 0
    snv_non_pass = 0
    incorrect_inbreeding_coeff = 0
    below_inbreeding_coeff = 0
    no_gt_col = 0
    failed_pre_filter_var_stats = 0
    fail_qc_post_filter = 0
    n_failed_cr = 0
    n_failed_maf = 0
    n_failed_hwe = 0
    monomorphic_post_filter = 0
    pass_qc = 0
    try:
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                variants += 1
                (_, reason, pre_filter_stats, post_filter_stats) = line.split("\t")
                if reason == "MultiAllelic":
                    multi_allelic += 1
                elif reason == "IndelBelowVQSR":
                    indel_below_vqsr += 1
                elif reason == "IndelNonPass":
                    indel_non_pass += 1
                elif reason == "SNVBelowVQSR":
                    snv_below_vqsr += 1
                elif reason == "SNVNonPass":
                    snv_non_pass += 1
                elif reason.startswith("IncorrectInbreedingCoeff"):
                    incorrect_inbreeding_coeff += 1
                elif reason.startswith("BelowInbreedingCoeff"):
                    below_inbreeding_coeff += 1
                elif reason == "NoGTCol":
                    no_gt_col += 1
                elif reason == "FailedPrefilterVarStats":
                    failed_pre_filter_var_stats += 1

                    splitted_stats = pre_filter_stats.split(";")
                    cr = float(splitted_stats[0].replace("CR=", ""))
                    maf = float(splitted_stats[1].replace("MAF=", ""))
                    hwe = float(splitted_stats[2].replace("HWE=", ""))

                    if cr <= thresh_cr:
                        n_failed_cr += 1
                    if maf <= thresh_maf:
                        n_failed_maf += 1
                    if hwe != -1 and hwe <= thresh_hwe:
                        n_failed_hwe += 1
                elif reason == "FailQCPostFilter":
                    fail_qc_post_filter += 1

                    splitted_stats = post_filter_stats.split(";")
                    cr = float(splitted_stats[0].replace("CR=", ""))
                    maf = float(splitted_stats[1].replace("MAF=", ""))
                    hwe = float(splitted_stats[2].replace("HWE=", ""))

                    if cr <= thresh_cr:
                        n_failed_cr += 1
                    if maf <= thresh_maf:
                        n_failed_maf += 1
                    if hwe != -1 and hwe <= thresh_hwe:
                        n_failed_hwe += 1
                elif reason == "MonomorphicPostFilter":
                    monomorphic_post_filter += 1
                elif reason == "PASSQC":
                    pass_qc += 1
                else:
                    pass
        f.close()
    except EOFError:
        pass

    return [variants, multi_allelic, indel_below_vqsr, indel_non_pass,
            snv_below_vqsr, snv_non_pass, incorrect_inbreeding_coeff,
            below_inbreeding_coeff, no_gt_col, failed_pre_filter_var_stats, fail_qc_post_filter,
            n_failed_cr, n_failed_maf, n_failed_hwe, monomorphic_post_filter, pass_qc,
            np.round(pass_qc / variants, 2)]


df = pd.DataFrame(np.nan, index=[i for i in range(len(args.input_paths) + 1)], columns=["Chromosome", "Variants", "MultiAllelic", "IndelBelowVQSR", "IndelNonPass", "SNVBelowVQSR", "SNVNonPass", "IncorrectInbreedingCoeff", "BelowInbreedingCoeff", "NoGTCol", "FailedPrefilterVarStats", "FailQCPostFilter", "FailedCR", "FailedMAF", "FailedHWE", "MonomorphicPostFilter", "PASSQC", "PctKept"])
df["Chromosome"] = df["Chromosome"].astype(str)
prev_filter_thresh = None
prev_failed_var_thresh = None
for i, input_path in enumerate(args.input_paths):
    df.iloc[i, 0] = os.path.basename(input_path).split("_")[1]

    df.iloc[i, 1:] = read_filter_logfile(
        filepath=input_path,
        thresh_maf=args.thresh_maf,
        thresh_cr=args.thresh_cr,
        thresh_hwe=args.thresh_hwe
    )
df.iloc[df.shape[0] - 1, 0] = "Total"
df.iloc[df.shape[0] - 1, 1:17] = df.iloc[:, 1:17].sum(axis=0)
df.iloc[df.shape[0] - 1, 17] = np.round(df.iloc[df.shape[0] - 1, 16] / df.iloc[df.shape[0] - 1, 1], 2)
df.iloc[:, 1:17] = df.iloc[:, 1:17].astype(int)
df.to_csv(args.output_path, sep="\t", header=True, index=False)
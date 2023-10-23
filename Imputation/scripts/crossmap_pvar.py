#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--pvar", type=str, required=True, help="The input pvar file.")
parser.add_argument("--bed", type=str, required=True, help="The crossmapped output bed file.")
parser.add_argument("--out_pvar", type=str, required=True, help="The crossmapped pvar file.")
parser.add_argument("--out_exclude", type=str, required=True, help="The variants to exclude.")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

crossmapped_info = {}
error = False
print("Loading variants loaded from {}".format(args.bed))
with open(args.bed, "r") as f:
    for line in f:
        values = line.strip("\n").split("\t")
        if values[3] in crossmapped_info:
            print("\tError, {} is duplicated".format(values[3]))
            error = True
        crossmapped_info[values[3]] = [values[0], values[1], values[4], values[5]]
f.close()
print("\t{} variants loaded".format(len(crossmapped_info)))

if error:
    exit()

print("Loading variants loaded from {}".format(args.pvar))
fh = open(args.pvar, "r")
fho_pvar = open(args.out_pvar, "w")
fho_excl = open(args.out_exclude, "w")

variants = set()
n_written = 0
n_unmapped = 0
n_duplicated = 0
n_exclude = 0
for line in fh:
    if line.startswith("#"):
        fho_pvar.write(line)
        continue

    chr, pos, variant_id, ref, alt = line.strip("\n").split("\t")

    crossmapped = False
    if variant_id in crossmapped_info:
        crossmapped = True
        chr, pos, ref, alt = crossmapped_info[variant_id]
        variant_id = chr + ":" + pos + ":" + ref + "_" + alt
    else:
        n_unmapped += 1

    duplicated = False
    if variant_id in variants:
        for i in range(1, 100):
            if variant_id + "_" + str(i) not in variants:
                variant_id = variant_id + "_" + str(i)
                duplicated = True
                n_duplicated += 1
                break

        if not duplicated:
            print("\tError, could not make unique ID for {}".format(variant_id))
            exit()

    if not crossmapped or duplicated:
        fho_excl.write(variant_id + "\n")
        n_exclude += 1

    fho_pvar.write("\t".join([chr, pos, variant_id, ref, alt]) + "\n")
    n_written += 1

    variants.add(variant_id)

fh.close()
fho_pvar.close()
fho_excl.close()

print("\t{} variants written to output".format(n_written))
print("\t{} variants written to exclude ({} unmapped, {} duplicated)".format(n_exclude, n_unmapped, n_duplicated))
print("End")

if error:
    exit()

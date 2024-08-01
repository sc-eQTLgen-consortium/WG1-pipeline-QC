#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import gzip

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


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


crossmapped_info = {}
error = False
print("Loading variants from {}".format(args.bed))
with gzopen(args.bed, mode="r") as f:
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

print("Loading variants from {}".format(args.pvar))
fh = gzopen(args.pvar, mode="r")
fho_pvar = gzopen(args.out_pvar, mode="w")
fho_excl = gzopen(args.out_exclude, mode="w")

variants = set()
n_written = 0
n_unmapped = 0
n_duplicated = 0
n_exclude = 0
for line in fh:
    if line.startswith("##"):
        fho_pvar.write(line)
        continue
    elif line.startswith("#CHROM"):
        # Validate header.
        values = line.strip("\n").split("\t")
        if (values[0] != "#CHROM") or (values[1] != "POS") or (values[2] != "ID") or (values[3] != "REF") or (values[4] != "ALT"):
            print("\tError, unexpected header.")
            exit()
        fho_pvar.write(line)
        continue

    values = line.strip("\n").split("\t")

    crossmapped = False
    if values[2] in crossmapped_info:
        crossmapped = True
        chr, pos, ref, alt = crossmapped_info[values[2]]
        values[0] = chr
        values[1] = pos
        values[2] = chr + ":" + pos + ":" + ref + "_" + alt
        values[3] = ref
        values[4] = alt
    else:
        n_unmapped += 1

    duplicated = False
    if values[2] in variants:
        for i in range(1, 100):
            if values[2] + "_" + str(i) not in variants:
                values[2] = values[2] + "_" + str(i)
                duplicated = True
                n_duplicated += 1
                break

        if not duplicated:
            print("\tError, could not make unique ID for {}".format(values[2]))
            exit()

    if not crossmapped or duplicated:
        fho_excl.write(values[2] + "\n")
        n_exclude += 1

    fho_pvar.write("\t".join(values) + "\n")
    n_written += 1

    variants.add(values[2])

fh.close()
fho_pvar.close()
fho_excl.close()

print("\t{} variants written to output".format(n_written))
print("\t{} variants written to exclude ({} unmapped, {} duplicated)".format(n_exclude, n_unmapped, n_duplicated))
print("End")

if error:
    exit()

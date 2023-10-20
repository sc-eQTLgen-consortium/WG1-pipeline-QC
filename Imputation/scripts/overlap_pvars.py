#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("--pvar1", type=str, required=True, help="The primary pvar file.")
parser.add_argument("--pvar2", type=str, required=True, help="The secondary pvar file.")
parser.add_argument("--variants1", type=str, required=True, help="The overlapping variants from the primary.")
parser.add_argument("--variants2", type=str, required=True, help="The overlapping variants from the secondary.")
args = parser.parse_args()

print("Options in effect:")
print("  --pvar1 {}".format(args.pvar1))
print("  --pvar2 {}".format(args.pvar2))
print("  --variants1 {}".format(args.variants1))
print("  --variants2 {}\n".format(args.variants2))

variants1 = {}
error = False
print("Loading variants loaded from {}".format(args.pvar1))
with open(args.pvar1, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        values = line.strip("\n").split("\t")
        variant1 = values[0] + ":" + values[1] + ":" + values[3] + "_" + values[4]
        if variant1 in variants1:
            print("Error, {} is duplicated".format(variant1))
            error = True
        variants1[variant1] = values[2]
f.close()
print("{} variants loaded".format(len(variants1)))

print("Loading variants loaded from {}".format(args.pvar2))
fh = open(args.pvar2, "r")
fho1 = open(args.variants1, "w")
fho2 = open(args.variants2, "w")

variants2 = set()
n_overlap = 0
for line in fh:
    if line.startswith("#"):
        continue

    values = line.strip("\n").split("\t")
    variant2 = values[0] + ":" + values[1] + ":" + values[3] + "_" + values[4]
    if variant2 in variants2:
        print("Error, {} is duplicated".format(variant1))
        error = True
    variants2.add(variant2)

    if variant2 in variants1:
        fho1.write(variants1[variant2] + "\n")
        fho2.write(values[2] + "\n")
        n_overlap += 1

fh.close()
fho1.close()
fho2.close()

print("{} variants loaded".format(len(variants2)))
print("{} variants written to output".format(n_overlap))
print("End")

if error:
    exit()

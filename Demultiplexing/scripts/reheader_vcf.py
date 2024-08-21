#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--vcf", required=True, type=str, help="")
parser.add_argument("--fai", required=True, type=str, help="")
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import gzip
import re

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


# Read the fasta index.
fai_contig_order = []
with gzopen(args.fai, mode="r") as f:
    for line in f:
        name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")
        fai_contig_order.append(name)
f.close()

# Load and save the original header.
fhi = gzopen(args.vcf, mode="r")
fho_ori = gzopen(args.out + ".old.hr", mode="w")
vcf_contigs = {}
valid = False
for line in fhi:
    if line.startswith("##contig"):
        contig_info = re.match("##contig=<(.+)>", line.rstrip("\n")).group(1).split(",")
        contig_id = None
        for contig_info_pair in contig_info:
            key, value = contig_info_pair.split("=")
            if key == "ID":
                contig_id = value
                break
        if contig_id is None:
            print("Error, contig '{}' has no ID.".format(line.rstrip("\n")))
            break
        if contig_id not in fai_contig_order:
            print("Error, contig '{}' not in fai.".format(line.rstrip("\n")))
            break
        vcf_contigs[contig_id] = line
        fho_ori.write(line)
    elif line.startswith("#"):
        fho_ori.write(line)
    else:
        valid = True
        break
fhi.close()
fho_ori.close()

if not valid:
    exit()

fho_ori = gzopen(args.out + ".old.hr", mode="r")
fho_new = gzopen(args.out + ".new.hr", mode="w")
ordered_contigs = [vcf_contigs[contig_id] for contig_id in fai_contig_order if contig_id in vcf_contigs]
contig_index = 0
for line in fho_ori:
    if line.startswith("##contig"):
        fho_new.write(ordered_contigs[contig_index])
        contig_index += 1
    elif line.startswith("##"):
        fho_new.write(line)
    else:
        # Add the command just before the header. SHould only happen once.
        fho_new.write("##reheader_vcf.py --vcf {} --fai {} --out {}\n".format(args.vcf, args.fai, args.out))
        fho_new.write(line)
fho_ori.close()
fho_new.close()

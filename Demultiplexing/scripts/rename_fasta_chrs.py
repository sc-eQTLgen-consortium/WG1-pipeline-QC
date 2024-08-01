#!/usr/bin/env python
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--fasta", required=True, type=str, help="")
parser.add_argument("--chr_name_conv", required=True, type=str, help="")
parser.add_argument("--ignore_missing_conv", action="store_true", default=False)
parser.add_argument("--out", required=True, type=str, help="")
args = parser.parse_args()


def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)


# Load the rename dictionary.
chr_trans = {}
fhi = gzopen(args.chr_name_conv, mode='r')
for line in fhi:
    old, new = line.rstrip("\n").split(" ")
    chr_trans[old] = new
fhi.close()

# Parse the fasta.
fhi = gzopen(args.fasta, mode='r')
fho = gzopen(args.out, mode='w')
valid = True
for line in fhi:
    if line.startswith(">"):
        line_data = line.lstrip(">").rstrip("\n").split(" ")
        chr = line_data[0]
        if chr not in chr_trans:
            print("{}, could not translate chr '{}'".format('Warning' if args.ignore_missing_conv else 'Error', chr))
            valid = False
        else:
            line_data[0] = chr_trans[chr]
        fho.write(">{}\n".format(" ".join(line_data)))
    else:
        fho.write(line)
fhi.close()
fho.close()

if not valid and not args.ignore_missing_conv:
    exit(1)
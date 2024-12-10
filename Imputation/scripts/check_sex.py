#!/usr/bin/env python3
# Author: M.Vochteloo
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", type=str, required=True, help="The primary pvar file.")
parser.add_argument("--outfile", type=str, required=True, help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

# Open input / output file handles.
fh_in = gzopen(args.input, mode="r")
fh_tab_out = gzopen(args.input + ".tsv", mode="w")
fh_manselect_out = gzopen(args.outfile, mode="w")

# Which columns we want in the manual select output file.
output_columns = ["FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F"]
# the output header name must be slightly different than the input column name
trans_column = {"FID": "#FID"}

print("Loading and reformatting input from {}".format(args.input))
pos = {}
nrows = 0
nproblem = 0
for i, line in enumerate(fh_in):
    values = line.rstrip("\n").split()

    # Write all data to the tab seperated output file.
    fh_tab_out.write("\t".join(values) + "\n")

    # Parse the header.
    if i == 0:
        pos = {column:index for index, column in enumerate(values)}
        if len(set(output_columns).difference(set(values))) != 0:
            print("Error, missing required columns.")
            exit()

        # Write the header of the manual select output file.
        fh_manselect_out.write("\t".join([trans_column[column] if column in trans_column else column for column in output_columns] + ["UPDATE/REMOVE/KEEP"]) + "\n")
        continue

    # Write the selected columns to the manual select output file, only if it is a problem.
    if values[pos["STATUS"]] == "PROBLEM":
        fh_manselect_out.write("\t".join([values[pos[column]] for column in output_columns] + [""]) + "\n")
        nproblem += 1

    nrows += 1

fh_in.close()
fh_tab_out.close()
fh_manselect_out.close()
print("\t{:,} rows loaded, {:,} problems found".format(nrows, nproblem))

print("End")

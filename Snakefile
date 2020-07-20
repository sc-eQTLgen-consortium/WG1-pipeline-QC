#!/usr/local/envs/py36/bin python3

import os
import sys
sys.path.append('mods') 
import pandas as pd
import prepareArguments
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["refs"]
popscle_dict = config["popscle"]
souporcell_dict = config["souporcell"]
DoubletDetection_dict = config["DoubletDetection"]
DoubletDetection_extra_dict = config["DoubletDetection"]
scrublet_dict = config["scrublet"]
scrublet_extra_dict = config["scrublet"]
scds_dict = config["scds"]
CombineResults_dict = config["CombineResults"]


# Use prepareArguments.py script to retrieve exact directories of single cell files
scrnaseq_libs = prepareArguments.get_scrnaseq_dirs(config)

# Get list of pools to process
samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")
samples.columns = ["Pool", "N"]



    

# Includes
include: 
    includes/Snakefile_popscle.smk
    includes/Snakefile_souporcell.smk
    includes/Snakefile_scrublet.smk
    includes/Snakefile_scds.smk
    includes/Snakefile_DoubletDetection.smk


rule all:
    input:
        # expand(output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/predicted_doublet_mask.txt", pool=samples.Pool, pctl=config["percentile"]),
        # expand(output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best",  pool=pools),
        # expand(output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv", pool=pools),
        if scrublet_dict["run_scrublet_manual"] == True:
            expand(output_dict["output_dir"] + , zip, pool = scrublet_dict["scrublet_manual_threshold_pools"], pctl = scrublet_dict["scrublet_manual_threshold_percentiles"])
        elif scrublet_dict["run_scrublet_manual"] == False:
            expand(output_dict["output_dir"] +  "/{pool}/CombinedResults/CombinedDropletAssignments.tsv", pool=pools),
            # expand(output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf", pool=pools),
            expand(output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz", pool=pools)


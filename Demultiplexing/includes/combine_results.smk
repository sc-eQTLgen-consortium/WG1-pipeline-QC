#!/usr/bin/env python
import pandas as pd
import os

def get_scrublet_input(wildcards):
    """
    This function selection the correct run of Scrublet we want to use for the combine_results rule.
    """
    if "Scrublet" in METHODS:
        man_select_path = config["outputs"]["output_dir"] + "manual_selections/Scrublet_percentile_manual_selection.tsv"
        if os.path.exists(man_select_path):
            # Find which thresholds were selected.
            logger.info("Read in the Scrublet manual selection file.")
            selection = pd.read_csv(man_select_path,sep="\t")
            selection["Pool"] = selection["Pool"].astype(str)
            selection["GeneVariabilityPctl"] = selection["GeneVariabilityPctl"].astype(str)
            select_dict = dict(zip(selection["Pool"], selection["GeneVariabilityPctl"]))

            return config["outputs"]["output_dir"] + "" + wildcards.pool + "/Scrublet_" + select_dict[wildcards.pool] + "/Scrublet_doublets_singlets.tsv"
    return []

rule combine_results:
    input:
        demuxlet = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best" if "popscle" in METHODS else [],
        souporcell = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv" if "souporcell" in METHODS else [],
        souporcell_assignments = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt" if "souporcell" in METHODS else [],
        doubletfinder = config["outputs"]["output_dir"] + "{pool}/DoubletFinder/DoubletFinder_doublets_singlets.tsv" if "DoubletFinder" in METHODS else [],
        scdblfinder = config["outputs"]["output_dir"] + "{pool}/scDblFinder/scDblFinder_doublets_singlets.tsv" if "scDblFinder" in METHODS else [],
        doubletdetection =  config["outputs"]["output_dir"] + "{pool}/DoubletDetection/DoubletDetection_doublets_singlets.tsv" if "DoubletDetection" in METHODS else [],
        scds = config["outputs"]["output_dir"] + "{pool}/scds/scds_doublets_singlets.tsv" if "scds" in METHODS else [],
        scrublet = get_scrublet_input,
    output:
        summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["combine_results"]["combine_results_time"]]
    threads: config["combine_results"]["combine_results_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Combine_Results.R",
        souporcell_correlation_limit = config["souporcell"]["souporcell_genotype_correlation_threshold"],
        out = config["outputs"]["output_dir"] + "{pool}/CombinedResults/"
    log: config["outputs"]["output_dir"] + "log/combine_results.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --demuxlet {input.demuxlet} \
            --souporcell {input.souporcell} \
            --souporcell_assignments {input.souporcell_assignments} \
            --souporcell_correlation_limit {params.souporcell_correlation_limit} \
            --doubletfinder {input.doubletfinder} \
            --scdblfinder {input.scdblfinder} \
            --doubletdetection {input.doubletdetection} \
            --scds {input.scds} \
            --scrublet {input.scrublet} \
            --out {params.out}
        """

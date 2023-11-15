#!/usr/bin/env python
import math

# def get_output_files():
#     output_files = []
#     if len([method for method in ["demuxlet", "freemuxlet", "scSplit", "souporcell", "vireo"] if method in METHODS]) > 0:
#         output_files.append(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results.tsv")
#         output_files.append(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv")
#
#     if len([method for method in ["Demuxlet", "Freemuxlet", "scSplit", "Souporcell", "Vireo", "DoubletDecon", "DoubletDetection", "DoubletFinder", "scDblFinder", "scds", "scrublet", "solo"] if method in METHODS]) > 0 and not math.isnan(config["combine_results_extra"]["method"]):
#         output_files.append(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_w_combined_assignments.tsv")
#         output_files.append(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Singlets_upset.pdf")
#
#     return output_files



rule combine_results:
    input:
        demuxlet = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best" if "popscle" in METHODS else [],
        souporcell = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz" if "souporcell" in METHODS else [],
        souporcell_assignments = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt.gz" if "souporcell" in METHODS else [],
        doubletfinder = config["outputs"]["output_dir"] + "{pool}/DoubletFinder/DoubletFinder_doublets_singlets.tsv.gz" if "DoubletFinder" in METHODS else [],
        scdblfinder = config["outputs"]["output_dir"] + "{pool}/scDblFinder/scDblFinder_doublets_singlets.tsv.gz" if "scDblFinder" in METHODS else [],
        doubletdetection = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=DOUBLETDETECTION_SELECTION[wildcards.pool]) if "DoubletDetection" in METHODS else [],
        scds = config["outputs"]["output_dir"] + "{pool}/scds/scds_doublets_singlets.tsv.gz" if "scds" in METHODS else [],
        scrublet = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCRUBLET_SELECTION[wildcards.pool]) if "Scrublet" in METHODS else []
    output:
        combined_results = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_summary.tsv",
        w_combined_assignments = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_w_combined_assignments.tsv",
        singlets_upset = config["outputs"]["output_dir"] + "{pool}/CombinedResults/Singlets_upset.pdf"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["combine_results_time"]]
    threads: config["combine_results"]["combine_results_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/scripts/Combine_Results.R",
        souporcell_correlation_limit = config["combine_results_extra"]["souporcell_correlation_limit"],
        pct_agreement = config["combine_results_extra"]["pct_agreement"],
        method = config["combine_results_extra"]["method"],
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
            --pct_agreement {params.pct_agreement} \
            --method {params.method} \
            --out {params.out}
        """

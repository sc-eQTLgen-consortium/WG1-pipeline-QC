#!/usr/bin/env python


rule combine_results:
    input:
        demuxlet = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best" if "popscle" in METHODS else [],
        souporcell = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz" if "souporcell" in METHODS else [],
        souporcell_assignments = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt.gz" if "souporcell" in METHODS else [],
        doubletfinder = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=DOUBLETFINDER_SELECTION[wildcards.pool]) if "DoubletFinder" in METHODS else [],
        scdblfinder = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDBLFINDER_SELECTION[wildcards.pool]) if "scDblFinder" in METHODS else [],
        doubletdetection = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=DOUBLETDETECTION_SELECTION[wildcards.pool]) if "DoubletDetection" in METHODS else [],
        scds = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDS_SELECTION[wildcards.pool]) if "scds" in METHODS else [],
        scrublet = lambda wildcards: config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCRUBLET_SELECTION[wildcards.pool]) if "Scrublet" in METHODS else []
    output:
        combined_results = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_summary.tsv",
        demultiplex_summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv" if "popscle" in METHODS or "souporcell" in METHODS else [],
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
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Combine_Results.R",
        demuxlet = "--demuxlet " + config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best" if "popscle" in METHODS else "",
        souporcell = "--souporcell " + config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz" if "souporcell" in METHODS else [],
        souporcell_assignments = "--souporcell_assignments " + config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt.gz" if "souporcell" in METHODS else[],
        souporcell_correlation_limit = config["combine_results_extra"]["souporcell_correlation_limit"],
        doubletfinder = lambda wildcards: "--DoubletFinder " + config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=DOUBLETFINDER_SELECTION[wildcards.pool]) if "DoubletFinder" in METHODS else[],
        scdblfinder = lambda wildcards: "--scDblFinder " + config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDBLFINDER_SELECTION[wildcards.pool]) if "scDblFinder" in METHODS else[],
        doubletdetection = lambda wildcards: "--DoubletDetection " + config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_doublets_singlets.tsv.gz".format(pool=wildcards.pool,run=DOUBLETDETECTION_SELECTION[wildcards.pool]) if "DoubletDetection" in METHODS else[],
        scds = lambda wildcards: "--scds " + config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDS_SELECTION[wildcards.pool]) if "scds" in METHODS else[],
        scrublet = lambda wildcards: "--scrublet " + config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz".format(pool=wildcards.pool,run=SCRUBLET_SELECTION[wildcards.pool]) if "Scrublet" in METHODS else[],
        pct_agreement = config["combine_results_extra"]["pct_agreement"],
        method = config["combine_results_extra"]["method"],
        out = config["outputs"]["output_dir"] + "{pool}/CombinedResults/"
    log: config["outputs"]["output_dir"] + "log/combine_results.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            {params.demuxlet} \
            {params.souporcell} \
            {params.souporcell_assignments} \
            --souporcell_correlation_limit {params.souporcell_correlation_limit} \
            {params.doubletfinder} \
            {params.scdblfinder} \
            {params.doubletdetection} \
            {params.scds} \
            {params.scrublet} \
            --pct_agreement {params.pct_agreement} \
            --method {params.method} \
            --out {params.out}
        """

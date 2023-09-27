#!/usr/bin/env python

#####################################
############ DoubletFinder ##########
#####################################
rule DoubletFinder:
    input:
        counts = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Counts"]
    output:
        doublets = config["outputs"]["output_dir"] + "{pool}/DoubletFinder/DoubletFinder_doublets_singlets.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/DoubletFinder/DoubletFinder_doublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["doubletfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["doubletfinder_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletfinder"]["doubletfinder_time"]]
    threads: config["doubletfinder"]["doubletfinder_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/DoubletFinder.R",
        out = config["outputs"]["output_dir"] + "{pool}/DoubletFinder/",
        find_clusters_resolution = config["doubletfinder_extra"]["find_clusters_resolution"],
        n_generated_artificial_doublets = config["doubletfinder_extra"]["n_generated_artificial_doublets"]
    log: config["outputs"]["output_dir"] + "log/DoubletFinder.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts} \
            --resolution {params.find_clusters_resolution} \
            --pN {params.n_generated_artificial_doublets}
        """

#####################################
############ scDblFinder ############
#####################################
rule scDblFinder:
    input:
        counts = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Counts"]
    output:
        doublets = config["outputs"]["output_dir"] + "{pool}/scDblFinder/scDblFinder_doublets_singlets.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/scDblFinder/scDblFinder_doublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scdblfinder"]["scfblfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scdblfinder"]["scfblfinder_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scdblfinder"]["scfblfinder_time"]]
    threads: config["scdblfinder"]["scdblfinder_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/scDblFinder.R",
        out = config["outputs"]["output_dir"] + "{pool}/scDblFinder/"
    log: config["outputs"]["output_dir"] + "log/scDblFinder.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts}
        """

###########################################
############ DOUBLET DETECTION ############
###########################################
# These steps uses different settings based on the user input in the yaml.
if config["doubletdetection_manual"]["run_doubletdetection_manual"]:
    logger.info("Running DoubletDetection with manual settings since 'run_doubletdetection_manual' is set to True.")
    doubletdetection_params_dict = {
        "log": "manual_rerun_variables.txt",
        "step": "manual",
        "n_iterations": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["doubletdetection_manual"]["doubletdetection_manual_pools"],config["doubletdetection_manual"]["doubletdetection_manual_n_iterations"])},
        "phenograph": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["doubletdetection_manual"]["doubletdetection_manual_pools"],config["doubletdetection_manual"]["doubletdetection_manual_phenograph"])},
        "standard_scaling": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["doubletdetection_manual"]["doubletdetection_manual_pools"],config["doubletdetection_manual"]["doubletdetection_manual_standard_scaling"])},
        "p_thresh": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["doubletdetection_manual"]["doubletdetection_manual_pools"],config["doubletdetection_manual"]["doubletdetection_manual_p_thresh"])},
        "voter_thresh": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["doubletdetection_manual"]["doubletdetection_manual_pools"],config["doubletdetection_manual"]["doubletdetection_manual_voter_thresh"])}
    }
else:
    logger.info("Running DoubletDetection with default settings since 'run_doubletdetection_manual' is set to False.")
    doubletdetection_params_dict = {
        "log": "default_run_variables.txt",
        "step": "default",
        "n_iterations": {pool: config["doubletdetection_extra"]["n_iterations"] for pool in SAMPLES},
        "phenograph": {pool: config["doubletdetection_extra"]["phenograph"] for pool in SAMPLES},
        "standard_scaling": {pool: config["doubletdetection_extra"]["standard_scaling"] for pool in SAMPLES},
        "p_thresh": {pool: config["doubletdetection_extra"]["p_thresh"] for pool in SAMPLES},
        "voter_thresh": {pool: config["doubletdetection_extra"]["voter_thresh"] for pool in SAMPLES}
    }


rule DoubletDetection:
    input:
        counts = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Counts"],
        barcodes = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"],
    output:
        doublets = config["outputs"]["output_dir"] + "{pool}/DoubletDetection/DoubletDetection_doublets_singlets.tsv",
        figure = report(config["outputs"]["output_dir"] + "{pool}/DoubletDetection/convergence_test.pdf", category = "DoubletDetection", subcategory = "{pool}", caption = "../report_captions/DoubletDetection.rst"),
        summary = config["outputs"]["output_dir"] + "{pool}/DoubletDetection/DoubletDetection_summary.tsv",
        log = config["outputs"]["output_dir"] + "{pool}/DoubletDetection/" + doubletdetection_params_dict["log"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["doubletdetection_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["doubletdetection_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletdetection"]["doubletdetection_time"]]
    threads: config["doubletdetection"]["doubletdetection_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/DoubletDetection.py",
        n_iterations = lambda wildcards: doubletdetection_params_dict["n_iterations"][wildcards.pool],
        phenograph = lambda wildcards: doubletdetection_params_dict["phenograph"][wildcards.pool],
        standard_scaling = lambda wildcards: doubletdetection_params_dict["standard_scaling"][wildcards.pool],
        p_thresh = lambda wildcards: doubletdetection_params_dict["p_thresh"][wildcards.pool],
        voter_thresh = lambda wildcards: doubletdetection_params_dict["voter_thresh"][wildcards.pool],
        step = doubletdetection_params_dict["step"],
        out = config["outputs"]["output_dir"] + "{pool}/DoubletDetection/"
    log: config["outputs"]["output_dir"] + "log/DoubletDetection.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --n_iterations {params.n_iterations} \
            --phenograph {params.phenograph} \
            --standard_scaling {params.standard_scaling} \
            --p_thresh {params.p_thresh} \
            --voter_thresh {params.voter_thresh} \
            --out {params.out} > {output.log}

		singularity exec --bind {params.bind} {params.sif} echo "The pool:" {wildcards.pool} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "This was a" {params.step} "run" >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The number of iterations used to determine doublets:" {params.n_iterations} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The phenograph was was used:" {params.phenograph} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The standard scaling was used:" {params.standard_scaling} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The p threshold was used:" {params.p_thresh} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The voter threshold is:" {params.voter_thresh} >> {output.log}
        """

##############################
############ SCDS ############
##############################
rule scds:
    input:
        counts = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Counts"]
    output:
        doublets = config["outputs"]["output_dir"] + "{pool}/scds/scds_doublets_singlets.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/scds/scds_doublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scds"]["scds_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scds"]["scds_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scds"]["scds_time"]]
    threads: config["scds"]["scds_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/scds.R",
        out = config["outputs"]["output_dir"] + "{pool}/scds/"
    log: config["outputs"]["output_dir"] + "log/scds.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts}
        """

##################################
############ SCRUBLET ############
##################################
# These steps uses different settings based on the user input in the yaml.
if config["scrublet_manual"]["run_scrublet_manual"]:
    logger.info("Running Scrublet with manual settings since 'run_scrublet_manual' is set to True.")
    pool_threshold_dict = dict(zip(config["scrublet_manual"]["scrublet_manual_threshold_pools"], config["scrublet_manual"]["scrublet_manual_threshold_thresholds"]))
    scrublet_params_dict = {
        "step": "manual",
        "scrublet_doublet_threshold": {str(pool): doublet_threshold for pool, doublet_threshold in zip(config["scrublet_manual"]["scrublet_manual_threshold_pools"], config["scrublet_manual"]["scrublet_manual_threshold_thresholds"])}
    }
else:
    logger.info("Running Scrublet with default settings since 'run_scrublet_manual' is set to False.")
    scrublet_params_dict = {
        "step": "default",
        "scrublet_doublet_threshold": {pool: None for pool in SAMPLES}
    }

rule Scrublet:
    input:
        counts = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Counts"],
        barcodes = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"]
    output:
        figure = report(config["outputs"]["output_dir"] + "{pool}/Scrublet_{pctl}/doublet_score_histogram.png", category = "Scrublet", caption = "../report_captions/Scrublet.rst", subcategory = "{pool}"),
        umap = report(config["outputs"]["output_dir"] + "{pool}/Scrublet_{pctl}/UMAP.png", category = "Scrublet", caption = "../report_captions/Scrublet.rst", subcategory = "{pool}"),
        results = config["outputs"]["output_dir"] + "{pool}/Scrublet_{pctl}/Scrublet_doublets_singlets.tsv",
        summary = config["outputs"]["output_dir"] + "{pool}/Scrublet_{pctl}/Scrublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["scrublet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["scrublet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scrublet"]["scrublet_time"]]
    threads: config["scrublet"]["scrublet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Scrublet_pipeline.py",
        sim_dbl = config["scrublet_extra"]["sim_dbl"],
        min_counts = config["scrublet_extra"]["min_counts"],
        min_cells = config["scrublet_extra"]["min_cells"],
        n_prin_comps = config["scrublet_extra"]["n_prin_comps"],
        out = config["outputs"]["output_dir"] + "{pool}/Scrublet_{pctl}/",
        step = lambda wildcards: scrublet_params_dict["step"],
        scrublet_doublet_threshold = lambda wildcards: scrublet_params_dict["scrublet_doublet_threshold"][wildcards.pool],
    log: config["outputs"]["output_dir"] + "log/Scrublet.{pool}.pctl{pctl}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --sim_doublet_ratio {params.sim_dbl} \
            --min_counts {params.min_counts} \
            --min_cells {params.min_cells} \
            --n_prin_comps {params.n_prin_comps} \
            --min_gene_variability_pctl {wildcards.pctl} \
            --scrublet_doublet_threshold {params.scrublet_doublet_threshold}
            --out {params.out}

            singularity exec --bind {params.bind} {params.sif} echo "The pool:" {wildcards.pool} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "This was a" {params.step} "run" >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of doublets simulated per droplet:" {params.sim_dbl} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The min number of counts used for filtering cells prior to PCA:" {params.min_counts} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of cells for a gene to be expressed in for filtering cells prior to PCA:" {params.min_cells} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of principle components used to embed the trnscriptomes prior to k-nearest-neighbor graph:" {params.n_prin_comps} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The manual doublet threshold set:" {params.scrublet_doublet_threshold} >> {log}
        [[ -s {output.results} ]]
        echo $?
        """
#!/usr/bin/env python
import math

#####################################
############ DoubletFinder ##########
#####################################
rule DoubletFinder:
    input:
        counts = lambda wildcards: POOL_DF.loc[wildcards.pool, "Counts"]
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_settings.json",
        doublets = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_doublets_singlets.tsv.gz",
        bcmvn = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_bcmvn.tsv.gz",
        figure = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/pKvBCmetric.png",
        unwated_figure = temp(config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/Rplots.pdf"),
        summary = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_doublet_summary.tsv",
        stats = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_stats.json"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["doubletfinder"]["doubletfinder_memory"] * config["doubletfinder"]["doubletfinder_threads"] - config["settings_extra"]["r_memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["doubletfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["doubletfinder_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletfinder"]["doubletfinder_time"]]
    threads: config["doubletfinder"]["doubletfinder_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/DoubletFinder.R",
        dims = lambda wildcards: "" if math.isnan(DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["dims"]) else "--dims " + str(DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["dims"]),
        resolution = lambda wildcards: DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["resolution"],
        expected_doublet_scaling_factor = lambda wildcards: "" if math.isnan(DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["expected_doublet_scaling_factor"]) else "--expected_doublet_scaling_factor " + str(DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["expected_doublet_scaling_factor"]),
        pn = lambda wildcards: DOUBLETFINDER_SETTINGS[wildcards.pool][wildcards.run]["pn"],
        out = config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/",
    log: config["outputs"]["output_dir"] + "log/DoubletFinder.{pool}.run{run}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --counts {input.counts} \
            {params.dims} \
            --resolution {params.resolution} \
            {params.expected_doublet_scaling_factor} \
            --num.cores {threads} \
            --pn {params.pn} \
            --mem {resources.mem} \
            --out {params.out}
        """

rule plot_DoubletFinder:
    input:
        settings = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_settings.json", run=DOUBLETFINDER_SETTINGS[wildcards.pool].keys(), allow_missing=True),
        bcmvns = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_bcmvn.tsv.gz", run=DOUBLETFINDER_SETTINGS[wildcards.pool].keys(), allow_missing=True)
    output:
        figure = report(config["outputs"]["output_dir"] + "QC_figures/{pool}/DoubletFinder_pKvBCmetrics.png", category="DoubletFinder", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/DoubletFinder.rst"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["plot_doubletfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletfinder"]["plot_doubletfinder_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletfinder"]["plot_doubletfinder_time"]]
    threads: config["doubletfinder"]["plot_doubletfinder_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/plot_DoubletFinder.py",
        out = config["outputs"]["output_dir"] + "QC_figures/{pool}/"
    log: config["outputs"]["output_dir"] + "log/plot_DoubletFinder.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --settings {input.settings} \
            --bcmvns {input.bcmvns} \
            --pool {wildcards.pool} \
            --out {params.out}
        """

#####################################
############ scDblFinder ############
#####################################
rule scDblFinder:
    input:
        counts = lambda wildcards: POOL_DF.loc[wildcards.pool, "Counts"]
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_settings.json",
        doublets = config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_doublets_singlets.tsv.gz",
        summary = config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_doublet_summary.tsv"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["scdblfinder"]["scfblfinder_memory"] * config["scdblfinder"]["scdblfinder_threads"] - config["settings_extra"]["r_memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scdblfinder"]["scfblfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scdblfinder"]["scfblfinder_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scdblfinder"]["scfblfinder_time"]]
    threads: config["scdblfinder"]["scdblfinder_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/scDblFinder.R",
        expected_doublet_scaling_factor = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["expected_doublet_scaling_factor"],
        stdev_doublet_rate = lambda wildcards: "" if math.isnan(SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["stdev_doublet_rate"]) else "--stdev_doublet_rate " + str(SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["stdev_doublet_rate"]),
        nfeatures = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["nfeatures"],
        keepUnidentifiable = lambda wildcards: "--keepUnidentifiable" if not SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["removeUnidentifiable"] else "",
        include_pcs = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["include_pcs"],
        prop_markers = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["prop_markers"],
        score = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["score"],
        processing = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["processing"],
        metric = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["metric"],
        nrounds = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["nrounds"],
        max_depth = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["max_depth"],
        iter = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["iter"],
        multi_sample_mode = lambda wildcards: SCDBLFINDER_SETTINGS[wildcards.pool][wildcards.run]["multi_sample_mode"],
        out = config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/"
    log: config["outputs"]["output_dir"] + "log/scDblFinder.{pool}.run{run}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --counts {input.counts} \
            --expected_doublet_scaling_factor {params.expected_doublet_scaling_factor} \
            {params.stdev_doublet_rate} \
            --nfeatures {params.nfeatures} \
            {params.keepUnidentifiable} \
            --include_pcs {params.include_pcs} \
            --prop_markers {params.prop_markers} \
            --score {params.score} \
            --processing {params.processing} \
            --metric {params.metric} \
            --nrounds {params.nrounds} \
            --max_depth {params.max_depth} \
            --iter {params.iter} \
            --multi_sample_mode {params.multi_sample_mode} \
            --mem {resources.mem} \
            --out {params.out}
        """

##########################################
############ DoubletDetection ############
##########################################
rule DoubletDetection:
    input:
        counts = lambda wildcards: POOL_DF.loc[wildcards.pool, "Counts"],
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"]
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_settings.json",
        doublets = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_doublets_singlets.tsv.gz",
        figure = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/convergence_test.pdf",
        summary = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_summary.tsv",
        pvalues = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_log_p_values.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["doubletdetection_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["doubletdetection_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletdetection"]["doubletdetection_time"]]
    threads: config["doubletdetection"]["doubletdetection_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/DoubletDetection.py",
        boost_rate = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["boost_rate"],
        n_components = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["n_components"],
        n_top_var_genes = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["n_top_var_genes"],
        replace = lambda wildcards: "--replace" if DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["replace"] else "",
        clustering_algorithm = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["clustering_algorithm"],
        n_iters = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["n_iters"],
        pseudocount = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["pseudocount"],
        standard_scaling = lambda wildcards: "--standard_scaling" if DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["standard_scaling"] else "",
        p_thresh = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["p_thresh"],
        voter_thresh = lambda wildcards: DOUBLETDETECTION_SETTINGS[wildcards.pool][wildcards.run]["voter_thresh"],
        out = config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/"
    log: config["outputs"]["output_dir"] + "log/DoubletDetection.{pool}.run{run}.log"
    shell:
        """        
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --boost_rate {params.boost_rate} \
            --n_components {params.n_components} \
            --n_top_var_genes {params.n_top_var_genes} \
            {params.replace} \
            --clustering_algorithm {params.clustering_algorithm} \
            --n_iters {params.n_iters} \
            --pseudocount {params.pseudocount} \
            {params.standard_scaling} \
            --n_jobs {threads} \
            --p_thresh {params.p_thresh} \
            --voter_thresh {params.voter_thresh} \
            --out {params.out}
        """


rule plot_DoubletDetection:
    input:
        settings = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_settings.json", run=DOUBLETDETECTION_SETTINGS[wildcards.pool].keys(), allow_missing=True),
        log_p_values = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_log_p_values.tsv.gz", run=DOUBLETDETECTION_SETTINGS[wildcards.pool].keys(), allow_missing=True)
    output:
        figure = report(config["outputs"]["output_dir"] + "QC_figures/{pool}/DoubletDetection_convergence_and_threshold_test.png", category="DoubletDetection", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/DoubletDetection.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["plot_doubletdetection_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["doubletdetection"]["plot_doubletdetection_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["doubletdetection"]["plot_doubletdetection_time"]]
    threads: config["doubletdetection"]["plot_doubletdetection_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/plot_DoubletDetection.py",
        out = config["outputs"]["output_dir"] + "QC_figures/{pool}/"
    log: config["outputs"]["output_dir"] + "log/plot_DoubletDetection.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --settings {input.settings} \
            --log_p_values {input.log_p_values} \
            --pool {wildcards.pool} \
            --out {params.out}
        """


##############################
############ SCDS ############
##############################
# Note: this method gives different results each time you run it.
rule scds:
    input:
        counts = lambda wildcards: POOL_DF.loc[wildcards.pool, "Counts"]
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_settings.json",
        doublets = config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_doublets_singlets.tsv.gz",
        summary = config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_doublet_summary.tsv"
    resources:
        mem = lambda wildcards, attempt: (attempt * config["scds"]["scds_memory"] * config["scds"]["scds_threads"] - config["settings_extra"]["r_memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scds"]["scds_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scds"]["scds_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scds"]["scds_time"]]
    threads: config["scds"]["scds_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/scds.R",
        bcds_ntop = lambda wildcards: SCDS_SETTINGS[wildcards.pool][wildcards.run]["bcds_ntop"],
        bcds_srat = lambda wildcards: SCDS_SETTINGS[wildcards.pool][wildcards.run]["bcds_srat"],
        bcds_nmax = lambda wildcards: SCDS_SETTINGS[wildcards.pool][wildcards.run]["bcds_nmax"],
        cxds_ntop = lambda wildcards: SCDS_SETTINGS[wildcards.pool][wildcards.run]["cxds_ntop"],
        cxds_binthresh = lambda wildcards: SCDS_SETTINGS[wildcards.pool][wildcards.run]["cxds_binthresh"],
        out = config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/"
    log: config["outputs"]["output_dir"] + "log/scds.{pool}.run{run}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --counts {input.counts} \
            --bcds_ntop {params.bcds_ntop} \
            --bcds_srat {params.bcds_srat} \
            --bcds_nmax {params.bcds_nmax} \
            --cxds_ntop {params.cxds_ntop} \
            --cxds_binthresh {params.cxds_binthresh} \
            --mem {resources.mem} \
            --out {params.out}
        """


##################################
############ SCRUBLET ############
##################################
# Gives slightly different results due to new software version.
rule Scrublet:
    input:
        counts = lambda wildcards: POOL_DF.loc[wildcards.pool, "Counts"],
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"]
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_settings.json",
        figure = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/doublet_score_histogram.png",
        umap = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/UMAP.png",
        results = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz",
        summary = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_summary.tsv",
        results_sim = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets_sim.tsv.gz",
        manifold = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_manifold.tsv.gz",
        manifold_sim = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_manifold_sim.tsv.gz",
        stats = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_stats.json"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["scrublet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["scrublet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scrublet"]["scrublet_time"]]
    threads: config["scrublet"]["scrublet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/Scrublet_pipeline.py",
        sim_doublet_ratio = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["sim_doublet_ratio"],
        n_neighbors = lambda wildcards: "" if math.isnan(SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["n_neighbors"]) else "--n_neighbors " + str(SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["n_neighbors"]),
        expected_doublet_scaling_factor = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["expected_doublet_scaling_factor"],
        stdev_doublet_rate = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["stdev_doublet_rate"],
        synthetic_doublet_umi_subsampling = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["synthetic_doublet_umi_subsampling"],
        get_doublet_neighbor_parents = lambda wildcards: "--get_doublet_neighbor_parents" if SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["get_doublet_neighbor_parents"] else "",
        min_counts = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["min_counts"],
        min_cells = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["min_cells"],
        min_gene_variability_pctl = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["min_gene_variability_pctl"],
        log_transform = lambda wildcards: "--log_transform" if SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["log_transform"] else "",
        no_mean_center = lambda wildcards: "--no_mean_center" if not SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["mean_center"] else "",
        no_normalize_variance = lambda wildcards: "--no_normalize_variance" if not SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["normalize_variance"] else "",
        n_prin_comps = lambda wildcards: SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["n_prin_comps"],
        doublet_threshold = lambda wildcards: "" if math.isnan(SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["doublet_threshold"]) else "--doublet_threshold " + str(SCRUBLET_SETTINGS[wildcards.pool][wildcards.run]["doublet_threshold"]),
        out = config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/",
    log: config["outputs"]["output_dir"] + "log/Scrublet.{pool}.run{run}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --sim_doublet_ratio {params.sim_doublet_ratio} \
            {params.n_neighbors} \
            --expected_doublet_scaling_factor {params.expected_doublet_scaling_factor} \
            --stdev_doublet_rate {params.stdev_doublet_rate} \
            --synthetic_doublet_umi_subsampling {params.synthetic_doublet_umi_subsampling} \
            {params.get_doublet_neighbor_parents} \
            --min_counts {params.min_counts} \
            --min_cells {params.min_cells} \
            --min_gene_variability_pctl {params.min_gene_variability_pctl} \
            {params.log_transform} \
            {params.no_mean_center} \
            {params.no_normalize_variance} \
            --n_prin_comps {params.n_prin_comps} \
            {params.doublet_threshold} \
            --out {params.out}
        """


rule plot_Scrublet:
    input:
        settings = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_settings.json", run=SCRUBLET_SETTINGS[wildcards.pool].keys(), allow_missing=True),
        results = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz", run=SCRUBLET_SETTINGS[wildcards.pool].keys(), allow_missing=True),
        stats = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_stats.json", run=SCRUBLET_SETTINGS[wildcards.pool].keys(), allow_missing=True),
        manifolds = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_manifold.tsv.gz", run=SCRUBLET_SETTINGS[wildcards.pool].keys(), allow_missing=True),
    output:
        figure = report(config["outputs"]["output_dir"] + "QC_figures/{pool}/Scrublet_histograms_and_UMAPs.png", category="Scrublet", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/Scrublet.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["plot_scrublet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["scrublet"]["plot_scrublet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["scrublet"]["plot_scrublet_time"]]
    threads: config["scrublet"]["plot_scrublet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/plot_Scrublet.py",
        out = config["outputs"]["output_dir"] + "QC_figures/{pool}/"
    log: config["outputs"]["output_dir"] + "log/plot_Scrublet.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --settings {input.settings} \
            --results {input.results} \
            --stats {input.stats} \
            --manifolds {input.manifolds} \
            --n_jobs {threads} \
            --pool {wildcards.pool} \
            --out {params.out}
        """
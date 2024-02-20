#!/usr/bin/env python


#################################
######## COMBINE RESULTS ########
#################################
# TODO: --assignment only works if there are no demultiplexing methods applied. Not sure how to implement this for multiplexed data yet.
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
        combined_results = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results.tsv.gz",
        summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_summary.tsv.gz",
        demultiplex_summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv.gz" if "popscle" in METHODS or "souporcell" in METHODS else [],
        w_combined_assignments = config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_w_combined_assignments.tsv.gz",
        droplet_assign_summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv.gz",
        droplet_assign_plot = report(config["outputs"]["output_dir"] + "{pool}/CombinedResults/DropletType_Assignment_BarPlot.png", category="Combine Results", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/final_assignments.rst"),
        singlets_upset = report(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Singlets_upset.pdf", category="Combine Results", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/combine_results.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["combine_results_time"]]
    threads: config["combine_results"]["combine_results_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/Combine_Results.R",
        demuxlet = "--demuxlet " + config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best" if "popscle" in METHODS else [],
        souporcell = "--souporcell " + config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz" if "souporcell" in METHODS else [],
        souporcell_assignments = "--souporcell_assignments " + config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt.gz" if "souporcell" in METHODS else [],
        souporcell_correlation_limit = config["combine_results_extra"]["souporcell_correlation_limit"],
        doubletfinder = lambda wildcards: "--DoubletFinder " + config["outputs"]["output_dir"] + "{pool}/DoubletFinderRun{run}/DoubletFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=DOUBLETFINDER_SELECTION[wildcards.pool]) if "DoubletFinder" in METHODS else [],
        scdblfinder = lambda wildcards: "--scDblFinder " + config["outputs"]["output_dir"] + "{pool}/scDblFinderRun{run}/scDblFinder_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDBLFINDER_SELECTION[wildcards.pool]) if "scDblFinder" in METHODS else [],
        doubletdetection = lambda wildcards: "--DoubletDetection " + config["outputs"]["output_dir"] + "{pool}/DoubletDetectionRun{run}/DoubletDetection_doublets_singlets.tsv.gz".format(pool=wildcards.pool,run=DOUBLETDETECTION_SELECTION[wildcards.pool]) if "DoubletDetection" in METHODS else [],
        scds = lambda wildcards: "--scds " + config["outputs"]["output_dir"] + "{pool}/scdsRun{run}/scds_doublets_singlets.tsv.gz".format(pool=wildcards.pool, run=SCDS_SELECTION[wildcards.pool]) if "scds" in METHODS else [],
        scrublet = lambda wildcards: "--scrublet " + config["outputs"]["output_dir"] + "{pool}/ScrubletRun{run}/Scrublet_doublets_singlets.tsv.gz".format(pool=wildcards.pool,run=SCRUBLET_SELECTION[wildcards.pool]) if "Scrublet" in METHODS else [],
        pct_agreement = config["combine_results_extra"]["pct_agreement"],
        method = config["combine_results_extra"]["method"],
        assignment = lambda wildcards: "--assignment " + ASSIGNMENT_COUPLING[wildcards.pool] if wildcards.pool in ASSIGNMENT_COUPLING else "",
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
            --pool {wildcards.pool} \
            {params.assignment} \
            --out {params.out}
        """


rule combine_pools:
    input:
        poolsheet = config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv",
        combined_results = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results.tsv.gz", pool=POOLS),
        summary = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_summary.tsv.gz", pool=POOLS),
        demultiplex_summary = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv.gz", pool=POOLS) if "popscle" in METHODS or "souporcell" in METHODS else [],
        w_combined_assignments = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_w_combined_assignments.tsv.gz", pool=POOLS),
        droplet_assign_summary = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv.gz", pool=POOLS)
    output:
        combined_results = config["outputs"]["output_dir"] + "CombinedResults/combined_results.tsv.gz",
        summary = config["outputs"]["output_dir"] + "CombinedResults/combined_results_summary.tsv.gz",
        demultiplex_summary = config["outputs"]["output_dir"] + "CombinedResults/combined_results_demultiplexing_summary.tsv.gz" if "popscle" in METHODS or "souporcell" in METHODS else [],
        w_combined_assignments = config["outputs"]["output_dir"] + "CombinedResults/combined_results_w_combined_assignments.tsv.gz",
        droplet_assign_summary = config["outputs"]["output_dir"] + "CombinedResults/Final_Assignments_demultiplexing_doublets.tsv.gz",
        done = config["outputs"]["output_dir"] + "CombinedResults/combine_pools.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_pools_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_pools_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["combine_pools_time"]]
    threads: config["combine_results"]["combine_pools_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/combine_pools.py",
        main_dir = config["outputs"]["output_dir"],
        out = config["outputs"]["output_dir"] + "CombinedResults/"
    log: config["outputs"]["output_dir"] + "log/combine_pools.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --poolsheet {input.poolsheet} \
            --main_dir {params.main_dir} \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """


####################################################
############ SCRIPTS TO PRODUCE QC PLOTS ############
####################################################
rule expected_observed_numbers:
    input:
        assignments = config["outputs"]["output_dir"] + "CombinedResults/Final_Assignments_demultiplexing_doublets.tsv.gz",
        poolsheet = config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv"
    output:
        figure = report(config["outputs"]["output_dir"] + "QC_figures/expected_observed_individuals_classifications.png", category="Number Individuals Summary", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/expected_observed_numbers.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["expected_observed_numbers_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["expected_observed_numbers_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["expected_observed_numbers_time"]]
    threads: config["combine_results"]["expected_observed_numbers_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/expected_observed_individuals_doublets.R",
        expected_doublet_scaling_factor = config["settings_extra"]["expected_doublet_scaling_factor"],
        out = config["outputs"]["output_dir"] + "QC_figures/"
    log: config["outputs"]["output_dir"] + "log/expected_observed_numbers.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --assignments {input.assignments} \
            --poolsheet {input.poolsheet} \
            --expected_doublet_scaling_factor {params.expected_doublet_scaling_factor} \
            --out {params.out}
        """


rule plot_singlet_doublet:
    input:
        assignments = config["outputs"]["output_dir"] + "CombinedResults/combined_results_w_combined_assignments.tsv.gz",
        poolsheet = config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv"
    output:
        figure = report(config["outputs"]["output_dir"] + "QC_figures/doublets_singlets.png", category="Number Individuals Summary", subcategory="doublets_singlets", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/singlet_doublet.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["plot_singlet_doublet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["plot_singlet_doublet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["plot_singlet_doublet_time"]]
    threads: config["combine_results"]["plot_singlet_doublet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/plot_singlet_doublet.py",
        basedir = config["outputs"]["output_dir"],
        out = config["outputs"]["output_dir"] + "QC_figures/"
    log: config["outputs"]["output_dir"] + "log/plot_singlet_doublet.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --assignments {input.assignments} \
            --poolsheet {input.poolsheet} \
            --out {params.out}
        """

#!/usr/bin/env python


#################################
######## COMBINE RESULTS ########
#################################
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
        droplet_assign_summary = config["outputs"]["output_dir"] + "{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv",
        droplet_assign_plot = report(config["outputs"]["output_dir"] + "{pool}/CombinedResults/DropletType_Assignment_BarPlot.png", category="Number Individuals Summary", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/final_assignments.rst"),
        singlets_upset = report(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Singlets_upset.pdf", category="Number Individuals Summary", subcategory="{pool}", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/combine_results.rst")
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
            --out {params.out}
        """


rule combine_pools:
    input:
        assignments = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/combined_results_w_combined_assignments.tsv", pool=POOLS),
        final_assignments = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv", pool=POOLS),
        poolsheet = config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv",
        rb_genes = config["refs"]["ref_dir"] + config["refs_extra"]["ribosomal_genes"],
        mt_genes = config["refs"]["ref_dir"] + config["refs_extra"]["mitochondrial_genes"]
    output:
        seurat_all = config["outputs"]["output_dir"] + "CombinedResults/seurat_object_all_pools_all_barcodes_all_metadata.rds",
        seurat_final = config["outputs"]["output_dir"] + "CombinedResults/seurat_object_all_pools_all_barcodes_final_assignments.rds",
        seurat_singlet_qc = config["outputs"]["output_dir"] + "CombinedResults/seurat_object_all_pools_singlet_barcodes_final_assignments.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_pools_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_pools_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["combine_pools_time"]]
    threads: config["combine_results"]["combine_pools_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/combine_pools.R",
        main_dir = config["outputs"]["output_dir"],
        out = config["outputs"]["output_dir"] + "CombinedResults/"
    log: config["outputs"]["output_dir"] + "log/combine_pools.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --main_dir {params.main_dir} \
            --rb_genes {input.rb_genes} \
            --mt_genes {input.mt_genes} \
            --out {params.out}
        """


####################################################
############ SCRIPTS TO PRODUCE QC PLOTS ############
####################################################
rule expected_observed_numbers:
    input:
        assignments = expand(config["outputs"]["output_dir"] + "{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv", pool=POOLS),
        poolsheet = config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv"
    output:
        report(config["outputs"]["output_dir"] + "QC_figures/expected_observed_individuals_classifications.png", category="Number Individuals Summary", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/expected_observed_numbers.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["expected_observed_numbers_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["expected_observed_numbers_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["expected_observed_numbers_time"]]
    threads: config["combine_results"]["expected_observed_numbers_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/expected_observed_individuals_doublets.R",
        basedir = config["outputs"]["output_dir"],
        out = config["outputs"]["output_dir"] + "QC_figures/"
    log: config["outputs"]["output_dir"] + "log/expected_observed_numbers.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --basedir {params.basedir} \
            --out {params.out}
        """

rule qc_plots:
    input:
        poolsheet=config["outputs"]["output_dir"] + "manual_selection/poolsheet.tsv",
        seurat = config["outputs"]["output_dir"] + "CombinedResults/seurat_object_all_pools_singlet_barcodes_final_assignments.rds",
    output:
        fig1 = report(config["outputs"]["output_dir"] + "QC_figures/nCount_RNA_violin_MAD_All.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_nCount.rst"),
        fig2 = report(config["outputs"]["output_dir"] + "QC_figures/nCount_RNA_violin_MADper_Pool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_nCount_RNA_MADall.rst"),
        fig3 = report(config["outputs"]["output_dir"] + "QC_figures/nCount_RNA_violin_noMADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_nCount.rst"),
        fig4 = report(config["outputs"]["output_dir"] + "QC_figures/nFeature_RNA_violin_MAD_All.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_feature_MADall.rst"),
        fig5 = report(config["outputs"]["output_dir"] + "QC_figures/nFeature_RNA_violin_MADper_Pool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_feature_MADperPool.rst"),
        fig6 = report(config["outputs"]["output_dir"] + "QC_figures/nFeature_RNA_violin_noMADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_feature.rst"),
        fig7 = report(config["outputs"]["output_dir"] + "QC_figures/nFeatures_vs_percentMT_QC_scatter_colorPool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_features_mt_pool.rst"),
        fig8 = report(config["outputs"]["output_dir"] + "QC_figures/nFeatures_vs_percentMT_QC_scatter_w_MADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_features_mt_MAD.rst"),
        fig9 = report(config["outputs"]["output_dir"] + "QC_figures/percent.mt_violin_MAD_All.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_mt_MADall.rst"),
        fig10 = report(config["outputs"]["output_dir"] + "QC_figures/percent.mt_violin_MADper_Pool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_mt_MADpool.rst"),
        fig11 = report(config["outputs"]["output_dir"] + "QC_figures/percent.mt_violin_noMADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_mt.rst"),
        fig12 = report(config["outputs"]["output_dir"] + "QC_figures/percent.rb_violin_MAD_All.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_rb_MADall.rst"),
        fig13 = report(config["outputs"]["output_dir"] + "QC_figures/percent.rb_violin_MADper_Pool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_rb_MADperPool.rst"),
        fig14 = report(config["outputs"]["output_dir"] + "QC_figures/percent.rb_violin_noMADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_rb.rst"),
        fig15 = report(config["outputs"]["output_dir"] + "QC_figures/UMI_vs_Genes_QC_scatter.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_UMI_features.rst"),
        fig16 = report(config["outputs"]["output_dir"] + "QC_figures/UMI_vs_Genes_QC_scatter_w_MADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_UMI_features_MADall.rst"),
        fig17 = report(config["outputs"]["output_dir"] + "QC_figures/UMI_vs_percentMT_QC_scatter_colorPool.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_UMI_mt_pool.rst"),
        fig18 = report(config["outputs"]["output_dir"] + "QC_figures/UMI_vs_percentMT_QC_scatter_w_MADlines.png", category="QC", caption=config["inputs"]["repo_dir"] + "Demultiplexing/report_captions/QC_plots_UMI_mt_MDA_all.rst"),
        done = config["outputs"]["output_dir"] + "QC_figures/qc_plots.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["qc_plots_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["qc_plots_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["qc_plots_time"]]
    threads: config["combine_results"]["qc_plots_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Demultiplexing/scripts/Singlet_QC_Figures.R",
        out = config["outputs"]["output_dir"] + "QC_figures/"
    log: config["outputs"]["output_dir"] + "log/QC_plots.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --poolsheet {input.poolsheet} \
            --seurat {input.seurat} \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """

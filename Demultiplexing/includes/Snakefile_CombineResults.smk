#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

#################################
######## COMBINE RESULTS ########
#################################
if os.path.exists(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"):
    scrublet_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv", sep = "\t")
    scrublet_selection["Pool"] = scrublet_selection["Pool"].astype(str)

rule join_results:
    input:
        demuxlet = output_dict["output_dir"] + "/{pool}/CombinedResults/demuxlet_results.txt",
        souporcell = output_dict["output_dir"] + "/{pool}/CombinedResults/souporcell_results.txt",
        scrublet = lambda wildcards: expand(output_dict["output_dir"] + "/{pool}/CombinedResults/{pctl}_scrublet_results.txt", zip, pool = wildcards.pool, pctl = scrublet_selection.scrublet_Percentile[scrublet_selection.Pool == wildcards.pool]),
        scds = output_dict["output_dir"] + "/{pool}/CombinedResults/scds_results.txt",
        DoubletDetection = output_dict["output_dir"] + "/{pool}/CombinedResults/DoubletDetection_results.txt"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path
    log: output_dict["output_dir"] + "/logs/join_results.{pool}.log"
    shell:
        """
        join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5" {input.demuxlet} {input.souporcell} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2,2.3" - {input.scrublet} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3" - {input.scds} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2" - {input.DoubletDetection} > {output}
        """



#####################################################################
############ SCRIPT FOR CALLING FINAL BARCODE ASSIGNMENT ############
#####################################################################
rule final_assignments:
    input:
        assignments = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"
    output:
        figure = report(output_dict["output_dir"] + "/{pool}/CombinedResults/DropletType_Assignment_BarPlot.png", category = "Number Individuals Summary", caption =  "../report_captions/final_assignments.rst"),
        table = output_dict["output_dir"] + "/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt",
        variables = output_dict["output_dir"] + "/{pool}/CombinedResults/variables.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalAssignments_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalAssignments_memory"]
    threads: CombineResults_dict["FinalAssignments_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path,
        out = output_dict["output_dir"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/FinalBarcodeAssignments.R"
    log: output_dict["output_dir"] + "/logs/final_assignments.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.out} > {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {wildcards.pool} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {output.variables}
        [[ -s {output.figure} ]]
        echo $?
        """

####################################################
############ SCRIPT TO PRODUCE QC PLOTS ############
####################################################
rule echo_final_assignments:
    input:
        expand(output_dict["output_dir"] + "/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt", pool = samples.Pool)
    output:
        output_dict["output_dir"] + "/QC_figures/final_assignments.txt"
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path
    log: output_dict["output_dir"] + "/logs/echo_final_assignments.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {input} | singularity exec --bind {params.bind} {params.sif} tr ' ' '\n' >> {output}
        """

rule final_assignments_check:
    input:
        assignment_list = output_dict["output_dir"] + "/QC_figures/final_assignments.txt",
        meta=input_dict["samplesheet_filepath"]
    output:
        assignment_list = output_dict["output_dir"] + "/QC_figures/final_assignments_comparison.txt",
        meta = output_dict["output_dir"] + "/QC_figures/meta_comparison.txt"
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path
    log: output_dict["output_dir"] + "/logs/final_assignments_check.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cat {input.assignment_list} > {output.assignment_list}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1}}' {input.meta} | singularity exec --bind {params.bind} {params.sif} tail -n+2 > {output.meta}
        if [ "$(wc -l < {output.meta})" -eq "$(wc -l < {output.assignment_list})" ]
        then 
            echo 0
        else 
            echo "The number of pools in the final_assignments.txt file don't match the number of pools in your sample sheet" > {log}
            rm {input.assignment_list}
            rm {output.assignment_list}
            rm {output.meta}
            echo 1
        fi
        """


rule expected_observed_numbers:
    input:
        final_assignment = output_dict["output_dir"] + "/QC_figures/final_assignments.txt",
        sample_sheet = input_dict["samplesheet_filepath"]
    output:
        report(output_dict["output_dir"] + "/QC_figures/expected_observed_individuals_classifications.png", category = "Number Individuals Summary", caption = "../report_captions/expected_observed_numbers.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["expected_observed_numbers_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["expected_observed_numbers_memory"]
    threads: CombineResults_dict["expected_observed_numbers_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path,
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/expected_observed_individuals_doublets.R",
        out = output_dict["output_dir"] + "/QC_figures/",
        basedir = output_dict["output_dir"]
    log: output_dict["output_dir"] + "/logs/expected_observed_numbers.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {input.sample_sheet} {params.out} {params.basedir}
        """


rule QC_plots:
    input:
        assignment_list = output_dict["output_dir"] + "/QC_figures/final_assignments_comparison.txt",
        pools = output_dict["output_dir"] + "/QC_figures/meta_comparison.txt"
    output:
        variables = temp(output_dict["output_dir"] + "/QC_figures/R_variables.txt"),
        fig1 = report(output_dict["output_dir"] + "/QC_figures/nCount_RNA_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_nCount.rst"),
        fig2 = report(output_dict["output_dir"] + "/QC_figures/nCount_RNA_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_nCount_RNA_MADall.rst"),
        fig3 = report(output_dict["output_dir"] + "/QC_figures/nCount_RNA_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_nCount.rst"),
        fig4 = report(output_dict["output_dir"] + "/QC_figures/nFeature_RNA_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_feature_MADall.rst"),
        fig5 = report(output_dict["output_dir"] + "/QC_figures/nFeature_RNA_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_feature_MADperPool.rst"),
        fig6 = report(output_dict["output_dir"] + "/QC_figures/nFeature_RNA_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_feature.rst"),
        fig7 = report(output_dict["output_dir"] + "/QC_figures/nFeatures_vs_percentMT_QC_scatter_colorPool.png", category = "QC", caption = "../report_captions/QC_plots_features_mt_pool.rst"),
        fig8 = report(output_dict["output_dir"] + "/QC_figures/nFeatures_vs_percentMT_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_features_mt_MAD.rst"),
        fig9 = report(output_dict["output_dir"] + "/QC_figures/percent.mt_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_mt_MADall.rst"),
        fig10 = report(output_dict["output_dir"] + "/QC_figures/percent.mt_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_mt_MADpool.rst"),
        fig11 = report(output_dict["output_dir"] + "/QC_figures/percent.mt_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_mt.rst"),
        fig12 = report(output_dict["output_dir"] + "/QC_figures/percent.rb_violin_MAD_All.png", category = "QC", caption = "../report_captions/QC_plots_rb_MADall.rst"),
        fig13 = report(output_dict["output_dir"] + "/QC_figures/percent.rb_violin_MADper_Pool.png", category = "QC", caption = "../report_captions/QC_plots_rb_MADperPool.rst"),
        fig14 = report(output_dict["output_dir"] + "/QC_figures/percent.rb_violin_noMADlines.png", category = "QC", caption = "../report_captions/QC_plots_rb.rst"),
        fig15 = report(output_dict["output_dir"] + "/QC_figures/UMI_vs_Genes_QC_scatter.png", category = "QC", caption = "../report_captions/QC_plots_UMI_features.rst"),
        fig16 = report(output_dict["output_dir"] + "/QC_figures/UMI_vs_Genes_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_UMI_features_MADall.rst"),
        fig17 = report(output_dict["output_dir"] + "/QC_figures/UMI_vs_percentMT_QC_scatter_colorPool.png", category = "QC", caption = "../report_captions/QC_plots_UMI_mt_pool.rst"),
        fig18 = report(output_dict["output_dir"] + "/QC_figures/UMI_vs_percentMT_QC_scatter_w_MADlines.png", category = "QC", caption = "../report_captions/QC_plots_UMI_mt_MDA_all.rst"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalQC_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * CombineResults_dict["FinalQC_memory"]
    threads: CombineResults_dict["FinalQC_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path,
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Singlet_QC_Figures.R",
        main_dir = output_dict["output_dir"],
        dirs10x = output_dict["output_dir"] + '/file_directories.txt',
        out = output_dict["output_dir"] + "/QC_figures/",
        rb_genes = "/opt/WG1-pipeline-QC/Demultiplexing/Ribosomal_genes.txt",
        mt_genes = "/opt/WG1-pipeline-QC/Demultiplexing/Mitochondrial_genes.txt"
    log: output_dict["output_dir"] + "/logs/QC_plots.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.main_dir} > {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.pools} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.dirs10x} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.out} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.rb_genes} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.mt_genes} >> {output.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {output.variables}
        [[ -s {output.fig18} ]]
        echo $?
        """

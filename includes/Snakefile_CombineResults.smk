#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["ref_dir"]
CombineResults_dict = config["CombineResults"]

# Use prepareArguments.py script to retrieve exact directories of single cell files
scrnaseq_filedict = prepareArguments.get_scrnaseq_dir(config)


#################################
######## COMBINE RESULTS ########
#################################
rule join_results:
    input:
        demuxlet = output_dict["output_dir"] + "/{pool}/CombinedResults/demuxlet_temp.txt",
        souporcell = output_dict["output_dir"] + "/{pool}/CombinedResults/souporcell_temp.txt",
        scrublet = output_dict["output_dir"] + "/{pool}/CombinedResults/scrublet_temp.txt",
        scds = output_dict["output_dir"] + "/{pool}/CombinedResults/scds_temp.txt",
        DoubletDetection = output_dict["output_dir"] + "/{pool}/CombinedResults/DoubletDetection_temp.txt"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
         join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5" {input.demuxlet} {input.souporcell} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2" - {input.scrublet} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3" - {input.scds} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.2" - {input.DoubletDetection} > {output}
        """



#####################################################################
############ SCRIPT FOR CALLING FINAL BARCODE ASSIGNMENT ############
#####################################################################
rule final_assignments:
    input:
        assignments = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"
    output:
        figure = output_dict["output_dir"] + "/{pool}/CombinedResults/DropletType_Assignment_BarPlot.png"
        table = output_dict["output_dir"] + "/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt"
    resources:
        mem_per_thread_gb = CombineResults_dict["FinalAssignments_memory"],
        disk_per_thread_gb = CombineResults_dict["FinalAssignments_memory"]
    threads: CombineResults_dict["FinalAssignments_threads"]
    params:
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"],
        script = input_dict["pipeline_dir"] + "/scripts/FinalBarcodeAssignments.R"
    shell:
        """
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {wildcards.pool} >> {output.variables}
        singularity exec {params.sif} Rscript {input.assignments} {output.variables}
        [[ -s {output.figure} ]]
        echo $?
        """

####################################################
############ SCRIPT TO PRODUCE QC PLOTS ############
####################################################
rule echo_final_assignments:
    input:
        output_dict["output_dir"] + "/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv"
    output:
        temp(output_dict["output_dir"] + "/QC_figures/final_assignments.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} echo {input} >> {output}
        """

rule final_assignments_check:
    input:
        assignment_list = output_dict["output_dir"] + "/QC_figures/final_assignments.txt"
        meta=input_dict["samplesheet_filepath"]
    output:
        assignment_list = temp(output_dict["output_dir"] + "/QC_figures/final_assignments_comparison.txt")
        meta = temp(output_dict["output_dir"] + "/QC_figures/meta_comparison.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} cat {input.assignment_list} > {output.assignment_list}
        singularity exec {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{print $1}' {input.meta} | singularity exec {params.sif} tail -n+2 > {output.meta}
        if [ "$(wc -l < {output.meta})" -eq "$(wc -l < {output.assignment_list})" ]
        then 
            echo 0
        else 
            rm {input.assignment_list}
            rm {output.assignment_list}
            rm {output.meta}
            echo 1
        fi
        """

rule QC_plots:
    input:
        assignment_list = output_dict["output_dir"] + "/QC_figures/final_assignments_comparison.txt"
        pools = output_dict["output_dir"] + "/QC_figures/meta_comparison.txt"
    output:
        variables = temp(output_dict["output_dir"] + "/QC_figures/R_variables.txt")
        fig = output_dict["output_dir"] + "/QC_figures/UMI_vs_Genes_QC_scatter.png"
    resources:
        mem_per_thread_gb = CombineResults_dict["FinalQC_memory"],
        disk_per_thread_gb = CombineResults_dict["FinalQC_memory"]
    threads: CombineResults_dict["FinalQC_threads"]
    params:
        sif = input_dict["singularity_image"],
        script = input_dict["pipeline_dir"] + "/scripts/Singlet_QC_Figures.R",
        main_dir = output_dict["output_dir"]
        dirs10x = scrnaseq_libs
        out = output_dict["output_dir"] + "/QC_figures/"
        rb_genes = input_dict["pipeline_dir"] + "/Ribosomal_genes.txt"
        mt_gene = input_dict["pipeline_dir"] + "/Mitochondrial_genes.txt
    shell:
        """
        singularity exec {params.sif} echo {params.main_dir} > {output.variables}
        singularity exec {params.sif} echo {input.pools} >> {output.variables}
        singularity exec {params.sif} echo {params.dirs10x} >> {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.rb_genes} >> {output.variables}
        singularity exec {params.sif} echo {params.mt_genes} >> {output.variables}
        singularity exec {params.sif} Rscript {input.script} {output.variables}
        [[ -s {output.fig} ]]
        echo $?
        """

#!/usr/bin/env python
import os
import pandas as pd
from glob import glob


##############################
############ SCDS ############
##############################
rule scds:
    input:
        script = input_dict["pipeline_dir"] + "/scripts/scds.R",
    output: 
        doublets= output_dict["output_dir"] + "/{pool}/scds/scds_doublets.txt",
        variables = output_dict["output_dir"] + "/{pool}/scds/scds_variables.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * scds_dict["scds_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * scds_dict["scds_memory"]
    threads: scds_dict["scds_threads"]
    params:
        matrix_dir = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
        out = output_dict["output_dir"] + "/{pool}/scds/",
        sif = input_dict["singularity_image"]
    log: output_dict["output_dir"] + "/logs/scds.{pool}.log"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} Rscript {input.script} {output.variables}
        """


rule scds_results_temp:
    input:
        output_dict["output_dir"] + "/{pool}/scds/scds_doublets.txt"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/scds_results.txt"
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input}  > {output}
        """

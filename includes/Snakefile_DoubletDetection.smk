#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["ref_dir"]
DoubletDetection_dict = config["scrublet"]
DoubletDetection_extra_dict = config["scrublet"]

# Use prepareArguments.py script to retrieve exact directories of single cell files
scrnaseq_filedict = prepareArguments.get_scrnaseq_dir(config)


###########################################
############ DOUBLET DETECTION ############
###########################################
if DoubletDetection_dict["run_DoubletDetection_manual"] == False:
    rule DoubletDetection:
        input:
            script = input_dict["pipeline_dir"] + "/scripts/DoubletDetection.py",
            matrix = output_dict["output_dir"] + "/{pool}/matrix_out/matrix.mtx.gz"
        output:
            doublets = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt",
            variables = temp(output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_variables.txt")
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"]
        threads: DoubletDetection_dict["DoubletDetection_threads"]
        params:
            matrix_dir = ,
            out = output_dict["output_dir"] + "/{pool}/DoubletDetection/",
            sif = input_dict["singularity_image"],
            n_iterations = DoubletDetection_extra_dict["n_iterations"],
            phenograph = DoubletDetection_extra_dict["phenograph"],
            standard_scaling = DoubletDetection_extra_dict["standard_scaling"],
            p_thresh = DoubletDetection_extra_dict["p_thresh"],
            voter_thresh = DoubletDetection_extra_dict["voter_thresh"],
        shell:
            """
            singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
            singularity exec {params.sif} echo {params.out} >> {output.variables}
            singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
            singularity exec {params.sif} python {input.script} \
                --counts_matrix {input.matrix} \
                --n_iterations {params.n_iterations} \
                --phenograph {params.phenograph}
                --standard_scaling {params.standard_scaling}
                --p_thresh {params.p_thresh}
                --voter_thresh {params.voter_thresh}
                -o {params.out}

            [[ -s {output.doublets} ]]
            echo $?
            """

elif DoubletDetection_dict["run_DoubletDetection_manual"] == True:
    rule DoubletDetection_manual_threshold:
        input:
            script = input_dict["pipeline_dir"] + "/scripts/DoubletDetection.py",
            matrix = output_dict["output_dir"] + "/{pool}/matrix_out/matrix.mtx.gz"
        output:
            doublets = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt",
            variables = temp(output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_variables.txt")
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * DoubletDetection_dict["DoubletDetection_memory"]
        threads: DoubletDetection_dict["DoubletDetection_threads"]
        params:
            matrix_dir = ,
            out = output_dict["output_dir"] + "/{pool}/DoubletDetection/",
            sif = input_dict["singularity_image"]
        shell:
            """
            singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
            singularity exec {params.sif} echo {params.out} >> {output.variables}
            singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
            singularity exec {params.sif} python {input.script} --counts_matrix {input.matrix} -o {params.out}
            [[ -s {output.doublets} ]]
            echo $?
            """
else:
    sys.exit()

rule DoubletDetection_create_selection_df:
    input:
        output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt"
    output:
        temporary = temp(output_dict["output_dir"] + "/manual_selections/DoubletDetection_PassFail_temp.txt"),
        temporary2 = temp(output_dict["output_dir"] + "/manual_selections/scrublet_gene_pctl_temp2.txt"),
        final = output_dict["output_dir"] + "/manual_selections/scrublet_gene_pctl.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} >> {output.temporary}
        singularity exec {params.sif} sort -u {output.temporary} > {output.temporary2}
        singularity exec {params.sif} sed "1 i Pool\tscrublet_Percentile" {output.temporary2} > {output.final}
        """

rule DoubletDetection_check_user_input:
    input:
        output_dict["output_dir"] + "/manual_selections/scrublet_gene_pctl.txt"
    output:
        user_input = temp(output_dict["output_dir"] + "/manual_selections/scrublet_user_input_temp.txt")
        pools = temp(output_dict["output_dir"] + "/manual_selections/pools")
        check_file = output_dict["output_dir"] + "/manual_selections/scrublet_check.done"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
    """
    singularity exec {params.sif} tail -n+2 {input} | \
        singularity exec {params.sif} awk '$2 > 0 && $2 <= 100' > {output.user_input}
    singularity exec {params.sif} tail -n+2 {input} > {output.pools}

    if [ "$(wc -l < {output.meta})" -eq "$(wc -l < {output.user_input})" ]
    then 
        singularity exec {params.sif} echo "Looks like you put percentile selections into the scrublet_user_input_temp.txt file."
        singularity exec {params.sif} echo 0
    else 
        singularity exec {params.sif} echo "You haven't put the scrublet gene percentile selection for each of the pools into the scrublet_gene_pctl.txt"
        singularity exec {params.sif} echo "Please check the scrublet outputs and choose the best variable genes percentile - rerun any of the pools where the thresholding failed (see the docs) to choose a manual threshold"
        singularity exec {params.sif} echo "Once you are happy with the thresholding, input the correct gene percentiles (as numbers between 0 and 100) into the second column of the scrublet_gene_pctl.txt file and restart the snakemake pipeline"
        singularity exec {params.sif} echo 1
    fi
    """

############################################################
############ REFORMAT DOUBLET DETECTION RESULTS ############
############################################################
rule DoubletDetection_results_temp:
    input:
        barcode = output_dict["output_dir"] + "/{pool}/matrix_out/barcodes.tsv",
        results = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt"
    output:
        DoubletDetection_temp=temp(output_dict["output_dir"] + "/{pool}/CombinedResults/DoubletDetection_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    params:
        sif = input_dict["singularity_image"]
    threads: 1
    shell:
        """
        sed "s/0.0/singlet/g" {input.results} | sed "s/1.0/doublet/g" | paste {input.barcode} - | sed "1 i Barcode\tDoubletDetection_DropletType" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.DoubletDetection_temp}
        """
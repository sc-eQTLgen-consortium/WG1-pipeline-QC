
#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["ref_dir"]
scrublet_dict = config["scrublet"]
scrublet_manual_dict = config["scrublet_manual"]
scrublet_extra_dict = config["scrublet_extra"]

# Use prepareArguments.py script to retrieve exact directories of single cell files
scrnaseq_filedict = prepareArguments.get_scrnaseq_dir(config)

# Get list of pools to process
samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")
pools = samples.iloc[:, 0]

def read_barcodes(pool, input_dict = None):
    # Use glob to grab the barcode filename within the filtered file directory, regardless of reference
    try:
        barcode_filename = glob(input_dict[pool] + "/outs/filtered_gene_bc_matrices/**/*/barcodes.tsv*", recursive = True)[0]
    except IndexError:
        print("Unable to find barcode files.")

    # Newer versions of Cell Ranger zip up the barcodes
    if (barcode_filename.endswith(".gz")):
        barcodes_df = pd.read_csv(barcode_filename, compression = "gz", header = None)
    else:
        barcodes_df = pd.read_csv(barcode_filename, header = None)

    # Turn this into a string - a bit convoluted but we had to deal with the gunzipped files just in case
    barcodes = barcodes_df.iloc[:,0]
    return(barcodes.str.cat(sep = "\n"))
    
    
##################################
############ SCRUBLET ############
##################################
if scrublet_dict["run_scrublet_manual"] == False:
    rule scrublet:
        input:
            matrix = ,
            genes = "
        output:
            output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/predicted_doublet_mask.txt"
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
        threads: scrublet_dict["scrublet_threads"]
        params:
            pctl = "{pctl}",
            sif = input_dict["singularity_image"],
            out = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/",
            script = input_dict["pipeline_dir"] + "/scripts/scrublet_pipeline.py",
            sim_dbl = scrublet_extra_dict["sim_dbl"],
            min_counts = scrublet_extra_dict["min_counts"],
            min_cells = scrublet_extra_dict["min_cells"],
            n_prin_comps = scrublet_extra_dict["n_prin_comps"]
        shell:
            """
            singularity exec {params.sif} python {params.script} \
                --counts_matrix {input.matrix} \
                --genes {input.genes} \
                --sim_doublet_ratio {params.sim_dbl} \
                --min_counts {params.min_counts} \
                --min_cells {params.min_cells} \
                --n_prin_comps {params.n_prin_comps} \
                --min_gene_variability_pctl {params.pctl} \
                -o {params.out}
            [[ -s {output} ]]
            echo $?
            """

elif scrublet_dict["run_scrublet_manual"] == True:
    rule scrublet_manual_threshold:
        input:
            matrix = output_dict["output_dir"] + "/{pool}/matrix_out/matrix.mtx",
            genes = output_dict["output_dir"] + "/{pool}/matrix_out/features.tsv"
        output:
            results = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/predicted_doublet_mask.txt"
            log = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/manual_rerun_variables.txt"
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
        threads: scrublet_dict["scrublet_threads"]
        params:
            pctl = "{pctl}",
            sif = input_dict["singularity_image"],
            out = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/",
            script = input_dict["pipeline_dir"] + "/scripts/scrublet_pipeline.py"
            sim_dbl = scrublet_extra_dict["sim_dbl"],
            min_counts = scrublet_extra_dict["min_counts"],
            min_cells = scrublet_extra_dict["min_cells"],
            n_prin_comps = scrublet_extra_dict["n_prin_comps"],
            scrublet_doublet_threshold = scrublet_manual_dict["scrublet_manual_threshold_thresholds"]
        shell:
            """
            singularity exec {params.sif} python {params.script} \
                --counts_matrix {input.matrix} \
                --genes {input.genes} \
                --sim_doublet_ratio {params.sim_dbl} \
                --min_counts {params.min_counts} \
                --min_cells {params.min_cells} \
                --n_prin_comps {params.n_prin_comps} \
                --min_gene_variability_pctl {params.pctl} \
                --scrublet_doublet_threshold {params.scrublet_doublet_threshold} \
                -o {params.out}
            [[ -s {output.log} ]]
            echo $?
            """
else:
    sys.exit()

rule scrublet_create_selection_df:
    input:
        output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/predicted_doublet_mask.txt"
    output:
        temporary = temp(output_dict["output_dir"] + "/manual_selections/scrublet_gene_pctl_temp.txt"),
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

rule scrublet_check_user_input:
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

rule scrublet_results_temp:
    input:
        barcode = output_dict["output_dir"] + "/{pool}/matrix_out/barcodes.tsv",
        scrublet_check = output_dict["output_dir"] + "/scrublet_check.done"
    output:
        scrublet_temp = temp(output_dict["output_dir"] + "/{pool}/CombinedResults/scrublet_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    params:
        path = output_dict["output_dir"] + "/{pool}/scrublet_",
        pctl = ,
        sif = input_dict["singularity_image"]
    threads: 1
    shell:
        """
        singularity exec {params.sif} paste {input.barcode} {params.path}{params.pctl}/predicted_doublet_mask.txt | \
            singularity exec {params.sif} sed "s/False/singlet/g" | \
            singularity exec {params.sif} sed "s/True/doublet/g" | \
            singularity exec {params.sif} sed "1 i Barcode\tscrublet_DropletType" | \
            singularity exec {params.sif} awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.scrublet_temp}
        """

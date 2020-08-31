
#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

    
    
##################################
############ SCRUBLET ############
##################################
rule make_scrublet_selection_df:
    input:
        input_dict["samplesheet_filepath"]
    output:
        output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"
    resources:
        mem_per_thread_gb = 1,
        disk_per_thread_gb = 1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    log: output_dict["output_dir"] + "/logs/make_scrublet_selection_df.log"
    shell:
        """
        singularity exec {params.sif} awk 'BEGIN{{OFS=FS="\t"}}{{print $1 "\t"}}' {input} | sed "1s/.*/Pool\tscrublet_Percentile/" > {output} 2> {log}
        """



if os.path.exists(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"):
    scrublet_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv", sep = "\t")
    if scrublet_selection["scrublet_Percentile"].count() != len(scrublet_selection):
        ready = False
    elif scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection):
        ready = True
        step = "ready"
    else:
        sys.exit()

    if (scrublet_manual_dict["run_scrublet_manual"] == False) or (scrublet_manual_dict["run_scrublet_manual"] == True and scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection)):
        log = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/default_run_variables.txt"
        sim_dbl = scrublet_extra_dict["sim_dbl"]
        min_counts = scrublet_extra_dict["min_counts"]
        min_cells = scrublet_extra_dict["min_cells"]
        n_prin_comps = scrublet_extra_dict["n_prin_comps"]
        step = "default"
        scrublet_doublet_threshold = None
    elif (scrublet_manual_dict["run_scrublet_manual"] == True): ### This deals with the possibility that the user still indicated that defaults need to be run but have completed the dataframe 
        log = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/manual_rerun_variables_" + datetime_now + ".txt"
        step = "manual"
        sim_dbl = scrublet_extra_dict["sim_dbl"]
        min_counts = scrublet_extra_dict["min_counts"]
        min_cells = scrublet_extra_dict["min_cells"]
        n_prin_comps = scrublet_extra_dict["n_prin_comps"]
        scrublet_doublet_threshold = lambda wildcards: scrublet_manual_dict["scrublet_manual_threshold_thresholds"][scrublet_manual_dict["scrublet_manual_threshold_pools"].index(wildcards.pool)]
    else:
        sys.exit()

    rule scrublet:
        input:
            matrix = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
            barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
            df = ancient(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv")
        output:
            results = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/scrublet_results.txt",
            log = log
        resources:
            mem_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
            disk_per_thread_gb = lambda wildcards, attempt: attempt * scrublet_dict["scrublet_memory"],
        threads: scrublet_dict["scrublet_threads"]
        params:
            sif = input_dict["singularity_image"],
            out = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/",
            script = input_dict["pipeline_dir"] + "/scripts/scrublet_pipeline.py",
            sim_dbl = sim_dbl,
            min_counts = min_counts,
            min_cells = min_cells,
            n_prin_comps = n_prin_comps,
            scrublet_doublet_threshold = scrublet_doublet_threshold,
            step = step,
            pipeline_dir = input_dict["pipeline_dir"],
            ready = ready
        log: output_dict["output_dir"] + "/logs/scrublet." + step + ".{pool}_{pctl}.log"
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                "No need to rerun scrublet since the parameters have alreeady been chosen. Will move on to sorting results and merging with results from all other softwares" 2> {log}
            elif [ {params.ready} == "False" ]
            then 
                if [ {params.step} == "default" ]
                then
                    singularity exec {params.sif} python {params.script} \
                        --counts_matrix {input.matrix} \
                        --barcodes {input.barcodes} \
                        --sim_doublet_ratio {params.sim_dbl} \
                        --min_counts {params.min_counts} \
                        --min_cells {params.min_cells} \
                        --n_prin_comps {params.n_prin_comps} \
                        --min_gene_variability_pctl {wildcards.pctl} \
                        -o {params.out} \
                        -d {params.pipeline_dir} 2> {log}
                elif [ {params.step} == "manual" ]
                then
                    singularity exec {params.sif} python {params.script} \
                        --counts_matrix {input.matrix} \
                        --barcodes {input.barcodes} \
                        --sim_doublet_ratio {params.sim_dbl} \
                        --min_counts {params.min_counts} \
                        --min_cells {params.min_cells} \
                        --n_prin_comps {params.n_prin_comps} \
                        --min_gene_variability_pctl {wildcards.pctl} \
                        -o {params.out} \
                        -d {params.pipeline_dir} \
                        --scrublet_doublet_threshold {params.scrublet_doublet_threshold} 2> {log}
                fi
                singularity exec {params.sif} echo "The pool:" {wildcards.pool} >> {log}
                singularity exec {params.sif} echo "This was a" {params.step} "run" >> {log}
                singularity exec {params.sif} echo "The number of doublets simulated per droplet:" {params.sim_dbl} >> {log}
                singularity exec {params.sif} echo "The min number of counts used for filtering cells prior to PCA:" {params.min_counts} >> {log}
                singularity exec {params.sif} echo "The number of cells for a gene to be expressed in for filtering cells prior to PCA:" {params.min_cells} >> {log}
                singularity exec {params.sif} echo "The number of principle components used to embed the trnscriptomes prior to k-nearest-neighbor graph:" {params.n_prin_comps} >> {log}
                singularity exec {params.sif} echo "The manual doublet threshold set:" {params.scrublet_doublet_threshold} >> {log}
            fi
            [[ -s {output.results} ]]
            echo $?
            """


    rule scrublet_check_user_input:
        input:
            df = ancient(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"),
            results = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/scrublet_results.txt"
        output:
            output_dict["output_dir"] + "/{pool}/CombinedResults/{pctl}_scrublet_results.txt"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
            disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
        threads: 1
        params:
            sif = input_dict["singularity_image"],
            ready = ready
        log: output_dict["output_dir"] + "/logs/scrublet_check_user_input.{pool}_{pctl}.log"
        shell:
            """
            if [ {params.ready} == "True" ]
            then 
                singularity exec {params.sif} echo "Looks like you put percentile selections into the scrublet_gene_pctl.txt file." 2> {log}
                singularity exec {params.sif} echo "The scrublet check is done and the next step of the pipeline will proceed" 2>> {log}
                singularity exec {params.sif} awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' {input.results} > {output}
            elif [ {params.ready} == "False" ]
            then
                singularity exec {params.sif} echo "You haven't put the scrublet gene percentile selection for each of the pools into the scrublet_gene_pctl.txt" 2> {log}
                singularity exec {params.sif} echo "Please check the scrublet outputs and choose the best variable genes percentile - rerun any of the pools where the thresholding failed (see the docs) to choose a manual threshold"
                singularity exec {params.sif} echo "Once you are happy with the thresholding, input the correct gene percentiles (as numbers between 0 and 100) into the second column of the scrublet_gene_pctl.txt file and restart the snakemake pipeline"
                singularity exec {params.sif} echo 1
            fi
            """

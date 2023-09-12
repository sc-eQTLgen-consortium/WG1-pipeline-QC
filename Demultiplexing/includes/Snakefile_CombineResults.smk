#################################
######## COMBINE RESULTS ########
#################################
rule combine_results:
    input:
        demuxlet = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["demuxlet"],
        souporcell = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["souporcell"],
        souporcell_assignments = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["souporcell_assignments"],
        doubletfinder = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["DoubletFinder"],
        scdblfinder = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["scDblFinder"],
        doubletdetection = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["DoubletDetection"],
        scds = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["scds"],
        scrublet = lambda wildcards: combine_results_extra_dict["combine_results_files"][wildcards.pool]["Scrublet"]
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * combine_results_dict["combine_results_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * combine_results_dict["combine_results_memory"]
    threads: combine_results_dict["combine_results_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Combine_Results.R",
        souporcell_correlation_limit = combine_results_extra_dict["souporcell_genotype_correlation_threshold"],
        out = output_dict["output_dir"] + "/{pool}/CombinedResults/"
    log: output_dict["output_dir"] + "/logs/combine_results.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --demuxlet {input.demuxlet} \
            --souporcell {input.souporcell} \
            --souporcell_assignments {input.souporcell_assignments} \
            --souporcell_correlation_limit {params.souporcell_correlation_limit} \
            --doubletfinder {input.doubletfinder} \
            --scdblfinder {input.scdblfinder} \
            --doubletdetection {input.doubletdetection} \
            --scds {input.scds} \
            --scrublet {input.scrublet} \
            --out {params.out}
        """
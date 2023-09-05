#################################
######## COMBINE RESULTS ########
#################################
rule combine_results:
    input:
        demuxlet = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best",
        souporcell = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv",
        souporcell_assignments = output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt",
        doubletdetection = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_doublets_singlets.tsv",
        doubletfinder = output_dict["output_dir"] + "/{pool}/DoubletFinder/DoubletFinder_doublets_singlets.tsv",
        scdblfinder = output_dict["output_dir"] + "/{pool}/scDblFinder/scDblFinder_doublets_singlets.tsv",
        scds = output_dict["output_dir"] + "/{pool}/scds/scds_doublets_singlets.tsv",
        scrublet = lambda wildcards: output_dict["output_dir"] + "/{pool}/scrublet_" + scrublet_select_dict[wildcards.pool] + "/scrublet_doublets_singlets.tsv"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/combined_results_demultiplexing_summary.tsv"
    resources:
        mem_per_thread_gb=lambda attempt: attempt * combine_results_dict["combine_results_memory"],
        disk_per_thread_gb=lambda attempt: attempt * combine_results_dict["combine_results_memory"]
    threads: combine_results_dict["combine_results_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Combine_Results.R",
        out = output_dict["output_dir"] + "/{pool}/CombinedResults/"
    log: output_dict["output_dir"] + "/logs/combine_results.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --demuxlet {input.demuxlet} \
            --souporcell {input.souporcell} \
            --souporcell_assignments {input.souporcell_assignments} \
            --doubletfinder {input.doubletfinder} \
            --scdblfinder {input.scdblfinder} \
            --scds {input.scds} \
            --scrublet {input.scrublet} \
            --out {params.out}
        """
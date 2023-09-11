#####################################
############ DoubletFinder ##########
#####################################
rule DoubletFinder:
    input:
        counts = lambda wildcards: scrnaseq_libs_df["CountH5File"][wildcards.pool]
    output: 
        doublets = output_dict["output_dir"] + "/{pool}/DoubletFinder/DoubletFinder_doublets_singlets.tsv",
        summary = output_dict["output_dir"] + "/{pool}/DoubletFinder/DoubletFinder_doublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * doubletfinder_dict["doubletfinder_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * doubletfinder_dict["doubletfinder_memory"]
    threads: doubletfinder_dict["doubletfinder_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/DoubletFinder.R",
        out = output_dict["output_dir"] + "/{pool}/DoubletFinder/",
        find_clusters_resolution = doubletfinder_extra_dict["find_clusters_resolution"],
        n_generated_artificial_doublets = doubletfinder_extra_dict["n_generated_artificial_doublets"]
    log: output_dict["output_dir"] + "/logs/DoubletFinder.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts} \
            --resolution {params.find_clusters_resolution} \
            --pN {params.n_generated_artificial_doublets}
        """

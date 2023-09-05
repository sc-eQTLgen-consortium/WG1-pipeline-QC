#####################################
############ scDblFinder ############
#####################################
rule scDblFinder:
    input:
        counts = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
    output: 
        doublets = output_dict["output_dir"] + "/{pool}/scDblFinder/scDblFinder_doublets.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * scdblfinder_dict["scfblfinder_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * scdblfinder_dict["scfblfinder_memory"]
    threads: scdblfinder_dict["scdblfinder_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/scDblFinder.R",
        out = output_dict["output_dir"] + "/{pool}/scDblFinder/"
    log: output_dict["output_dir"] + "/logs/scDblFinder.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts}
        """

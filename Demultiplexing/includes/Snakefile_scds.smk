##############################
############ SCDS ############
##############################
rule scds:
    input:
        counts = lambda wildcards: scrnaseq_libs_df["CountH5File"][wildcards.pool],
    output: 
        doublets = output_dict["output_dir"] + "/{pool}/scds/scds_doublets_singlets.tsv",
        summary = output_dict["output_dir"] + "/{pool}/scds/scds_doublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * scds_dict["scds_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * scds_dict["scds_memory"]
    threads: scds_dict["scds_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/scds.R",
        out = output_dict["output_dir"] + "/{pool}/scds/"
    log: output_dict["output_dir"] + "/logs/scds.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --out {params.out} \
            --counts {input.counts}
        """
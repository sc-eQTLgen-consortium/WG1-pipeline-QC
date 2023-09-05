##################################
############ SCRUBLET ############
##################################
rule scrublet:
    input:
        counts = lambda wildcards: scrnaseq_libs_df["Matrix_Directories"][wildcards.pool],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        df = ancient(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv")
    output:
        figure = report(output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/doublet_score_histogram.png", category = "Scrublet", caption = "../report_captions/scrublet.rst", subcategory = "{pool}")
        umap = report(output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/UMAP.png", category = "Scrublet", caption = "../report_captions/scrublet.rst", subcategory = "{pool}")
        results = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/scrublet_doublets_singlets.tsv",
        summary = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/scrublet_summary.tsv"
    resources:
        mem_per_thread_gb = lambda attempt: attempt * scrublet_dict["scrublet_memory"],
        disk_per_thread_gb = lambda attempt: attempt * scrublet_dict["scrublet_memory"],
    threads: scrublet_dict["scrublet_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/scrublet_pipeline.py",
        sim_dbl = scrublet_extra_dict["sim_dbl"],
        min_counts = scrublet_extra_dict["min_counts"],
        min_cells = scrublet_extra_dict["min_cells"],
        n_prin_comps = scrublet_extra_dict["n_prin_comps"],
        out = output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/",
        step = scrublet_params_dict["step"],
        scrublet_doublet_threshold = lambda wildcards: scrublet_params_dict["scrublet_doublet_threshold"][wildcards.pool]
    log: output_dict["output_dir"] + "/logs/scrublet.{pool}.pctl{pctl}.log"
    shell:
    """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --sim_doublet_ratio {params.sim_dbl} \
            --min_counts {params.min_counts} \
            --min_cells {params.min_cells} \
            --n_prin_comps {params.n_prin_comps} \
            --min_gene_variability_pctl {wildcards.pctl} \
            --out {params.out} \
            [[ ${params.step} == "manual" ]] && { '--scrublet_doublet_threshold {params.scrublet_doublet_threshold}'}

            singularity exec --bind {params.bind} {params.sif} echo "The pool:" {wildcards.pool} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "This was a" {params.step} "run" >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of doublets simulated per droplet:" {params.sim_dbl} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The min number of counts used for filtering cells prior to PCA:" {params.min_counts} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of cells for a gene to be expressed in for filtering cells prior to PCA:" {params.min_cells} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The number of principle components used to embed the trnscriptomes prior to k-nearest-neighbor graph:" {params.n_prin_comps} >> {log}
            singularity exec --bind {params.bind} {params.sif} echo "The manual doublet threshold set:" {params.scrublet_doublet_threshold} >> {log}
        [[ -s {output.results} ]]
        echo $?
        """




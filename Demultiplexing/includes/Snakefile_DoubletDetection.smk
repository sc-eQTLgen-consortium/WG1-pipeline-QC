###########################################
############ DOUBLET DETECTION ############
###########################################
rule DoubletDetection:
    input:
        counts = lambda wildcards: scrnaseq_libs_df["CountH5File"][wildcards.pool],
        barcodes = lambda wildcards: scrnaseq_libs_df["BarcodeFile"][wildcards.pool],
        df = ancient(output_dict["output_dir"] + "/manual_selections/DoubletDetection_manual_selection.tsv")
    output:
        doublets = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_doublets_singlets.tsv",
        figure = report(output_dict["output_dir"] + "/{pool}/DoubletDetection/convergence_test.pdf", category = "DoubletDetection", subcategory = "{pool}", caption = "../report_captions/DoubletDetection.rst"),
        summary = output_dict["output_dir"] + "/{pool}/DoubletDetection/DoubletDetection_summary.tsv",
        log = output_dict["output_dir"] + "/{pool}/DoubletDetection/" + doubletdetection_params_dict["log"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * doubletdetection_dict["doubletdetection_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * doubletdetection_dict["doubletdetection_memory"]
    threads: doubletdetection_dict["doubletdetection_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/DoubletDetection.py",
        n_iterations = lambda wildcards: doubletdetection_params_dict["n_iterations"][wildcards.pool],
        phenograph = lambda wildcards: doubletdetection_params_dict["phenograph"][wildcards.pool],
        standard_scaling = lambda wildcards: doubletdetection_params_dict["standard_scaling"][wildcards.pool],
        p_thresh = lambda wildcards: doubletdetection_params_dict["p_thresh"][wildcards.pool],
        voter_thresh = lambda wildcards: doubletdetection_params_dict["voter_thresh"][wildcards.pool],
        step = doubletdetection_params_dict["step"],
        out = output_dict["output_dir"] + "/{pool}/DoubletDetection/"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --n_iterations {params.n_iterations} \
            --phenograph {params.phenograph} \
            --standard_scaling {params.standard_scaling} \
            --p_thresh {params.p_thresh} \
            --voter_thresh {params.voter_thresh} \
            --out {params.out} > {output.log}

		singularity exec --bind {params.bind} {params.sif} echo "The pool:" {wildcards.pool} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "This was a" {params.step} "run" >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The number of iterations used to determine doublets:" {params.n_iterations} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The phenograph was was used:" {params.phenograph} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The standard scaling was used:" {params.standard_scaling} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The p threshold was used:" {params.p_thresh} >> {output.log}
		singularity exec --bind {params.bind} {params.sif} echo "The voter threshold is:" {params.voter_thresh} >> {output.log}
        """
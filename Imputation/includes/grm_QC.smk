#!/usr/local/envs/py36/bin python3


rule grm:
    input:
        bim = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pvar",
        bed = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pgen",
        fam = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.psam"
    output:
        bim = output_dict["output_dir"] + "/grm/{ancestry}_grm.pvar",
        bed = output_dict["output_dir"] + "/grm/{ancestry}_grm.pgen",
        fam = output_dict["output_dir"] + "/grm/{ancestry}_grm.psam",
        ypar = output_dict["output_dir"] + "/grm/{ancestry}_y_par.snps",
        Nbin = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.N.bin",
        id = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.id",
        bin = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.bin"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_memory"]
    threads: grm_QC_dict["grm_threads"]
    params:
        out = output_dict["output_dir"] + "/grm/{ancestry}_grm",
        het_filter = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk '$1 == "2324"' {input.bim} | singularity exec --bind {params.bind} {params.sif} cut -f2 > {output.ypar} ### Y SNPs are coded as XY => after changing to numbers = 2324
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.het_filter} --exclude {output.ypar} --make-pgen --out {params.out}
        singularity exec --bind {params.bind} {params.sif} gcta64 --pfile {params.out} --make-grm --out {params.out}
        """

rule report_relatedness:
    input:
        Nbin = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.N.bin",
        id = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.id",
        bin = output_dict["output_dir"] + "/grm/{ancestry}_grm.grm.bin"
    output:
        output_dict["output_dir"] + "/grm/{ancestry}_paired_relatedness.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["report_relatedness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["report_relatedness_memory"]
    threads: grm_QC_dict["report_relatedness_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/grm/{ancestry}_grm",
        out = output_dict["output_dir"] + "/grm/{ancestry}_",
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-pipeline-QC/Imputation/scripts/viewGRM.R"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.infile} {params.out}
        """

#!/usr/local/envs/py36/bin python3


rule grm:
    input:
        bim = output_dict["output_dir"] + "/het_filter/het_filter.bim",
        bed = output_dict["output_dir"] + "/het_filter/het_filter.bed",
        fam = output_dict["output_dir"] + "/het_filter/het_filter.fam"
    output:
        bim = output_dict["output_dir"] + "/grm/grm.bim",
        bed = output_dict["output_dir"] + "/grm/grm.bed",
        fam = output_dict["output_dir"] + "/grm/grm.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_memory"]
    threads: grm_QC_dict["grm_threads"]
    params:
        out = output_dict["output_dir"] + "/grm",
        het_filter = output_dict["output_dir"] + "/het_filter/het_filter",
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk '$1 == "2324"' {input.bim} | singularity exec --bind {params.bind} {params.sif} cut -f2 > {params.out}/y_par.snps ### Y SNPs are coded as XY => after changing to numbers = 2324
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.het_filter} --exclude {params.out}/y_par.snps --make-bed --out {params.out}/grm --noweb
        singularity exec --bind {params.bind} {params.sif} gcta64 --bfile {params.out}/grm --make-grm --out {params.out}/grm
        """

rule grm_filter:
    input:
        bim = output_dict["output_dir"] + "/grm/grm.bim",
        bed = output_dict["output_dir"] + "/grm/grm.bed",
        fam = output_dict["output_dir"] + "/grm/grm.fam"
    output:
        grm_id = output_dict["output_dir"] + "/grm_filter/grm_filter.grm.id",
        # bim = output_dict["output_dir"] + "/grm_filter/grm_filter.bim",
        # bed = output_dict["output_dir"] + "/grm_filter/grm_filter.bed",
        # fam = output_dict["output_dir"] + "/grm_filter/grm_filter.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_filter_memory"]
    threads: grm_QC_dict["grm_filter_threads"]
    params:
        grm = output_dict["output_dir"] + "/grm/grm",
        grm_filter = output_dict["output_dir"] + "/grm_filter/grm_filter",
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        relatedness = grm_QC_dict["grm_filter_relatedness"],
        pcs = grm_QC_dict["grm_filter_pcs"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gcta64 --grm {params.grm} --grm-cutoff 0.9 --make-grm --out {params.grm_filter}
        singularity exec --bind {params.bind} {params.sif} gcta64 --grm {params.grm_filter} --pca 6 --out {params.grm_filter}
        """

# rule grm_singleton:
#     input:
#         bim = output_dict["output_dir"] + "/grm/grm.bim",
#         bed = output_dict["output_dir"] + "/grm/grm.bed",
#         fam = output_dict["output_dir"] + "/grm/grm.fam"
#     output:
#         family = output_dict["output_dir"] + "/grm_singleton/grm_singleton.family.txt",
#         singleton = output_dict["output_dir"] + "/grm_singleton/grm_singleton.singleton.txt"
#     resources:
#         mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
#         disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
#     threads: 1
#     params:
#         sif = input_dict["singularity_image"],
#         bind = input_dict["bind_paths"],
#         grm_singleton = output_dict["output_dir"] + "/grm_singleton/grm_singleton",
#         grm = output_dict["output_dir"] + "/grm/grm"
#     shell:
#         """
#         singularity exec --bind {params.bind} {params.sif} gcta64 --grm {params.grm} --grm-singleton 0.05 --out {params.grm_singleton}
#         """


rule grm_subset:
    input:
        grm_filter_id = output_dict["output_dir"] + "/grm_filter/grm_filter.grm.id"
    output:
        bim = output_dict["output_dir"] + "/grm_subset/grm_subset.bim",
        bed = output_dict["output_dir"] + "/grm_subset/grm_subset.bed",
        fam = output_dict["output_dir"] + "/grm_subset/grm_subset.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_subset_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * grm_QC_dict["grm_subset_memory"]
    threads: grm_QC_dict["grm_subset_threads"]
    params:
        bind = input_dict["bind_paths"],
        het_filter = output_dict["output_dir"] + "/het_filter/het_filter",
        grm_subset = output_dict["output_dir"] + "/grm_subset/grm_subset",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.het_filter} --keep {input.grm_filter_id} --make-bed --out  {params.grm_subset}
        """

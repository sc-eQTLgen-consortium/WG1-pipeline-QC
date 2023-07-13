#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

###################################
############# POPSCLE #############
###################################

###### popscle bam filter ######
rule popscle_bam_filter:
    input:
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        vcf = input_dict["snp_genotypes_filepath"],
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
    output:
        output_dict["output_dir"] + "/{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam"
    resources:
        bam_filter_threads=popscle_dict["bam_filter_threads"],
        mem_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["bam_filter_memory"],
        disk_per_thread_gb=0
    threads: popscle_dict["bam_filter_threads"]
    params:
        tag_group = popscle_dict["tag_group"],
        sif = input_dict["singularity_image"],
        bind = bind_path,
        out_dir = output_dict["output_dir"] + "/{pool}/popscle/bam_filter/",
        out = output_dict["output_dir"] + "/{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam"
    log: output_dict["output_dir"] + "/logs/popscle_bam_filter.{pool}.log"
    shell:
        """
        mkdir -p {params.out_dir} && \
        singularity exec --bind {params.bind} {params.sif} bedtools merge -i {input.vcf} \
            | singularity exec --bind {params.bind} {params.sif} samtools view \
                -L - \
                -D {params.tag_group}:<(zcat {input.barcodes}) \
                -o {params.out} \
                --write-index \
                -@ {resources.bam_filter_threads} \
                {input.bam}
        """

###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        vcf = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        bam = lambda wildcards: output_dict["output_dir"] + "/{wildcards.pool}/popscle/bam_filter/{wildcards.pool}_snpfiltered_alignment.bam",
#        bam = expand(output_dict["output_dir"] + "/popscle/bam_filter/{pool}_snpfiltered_alignment.bam", pool=samples.Pool),
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"]
    threads: popscle_dict["pileup_threads"]
    params:
        out = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup",
        sif = input_dict["singularity_image"],
        bind = bind_path, 
        tag_group = popscle_dict["tag_group"],
        tag_UMI = popscle_dict["tag_UMI"],
        cap_BQ = popscle_extra_dict["cap_BQ"],
        min_BQ = popscle_extra_dict["min_BQ"],
        min_MQ = popscle_extra_dict["min_MQ"],
        min_TD = popscle_extra_dict["min_TD"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"]
    log: output_dict["output_dir"] + "/logs/popscle_pileup.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle dsc-pileup \
            --sam {input.bam} \
            --tag-group {params.tag_group} \
            --tag-UMI {params.tag_UMI} \
            --sm-list {input.individuals} \
            --vcf {input.vcf} \
            --cap-BQ {params.cap_BQ} \
            --min-BQ {params.min_BQ} \
            --min-MQ {params.min_MQ} \
            --min-TD {params.min_TD} \
            --excl-flag {params.excl_flag} \
            --group-list {input.barcodes} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}
        """

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz",
        snps = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_threads"]
    threads: popscle_dict["demuxlet_threads"]
    params:
        pileup = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup",
        out = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/",
        sif = input_dict["singularity_image"],
        bind = bind_path,
        field = popscle_dict["genotype_field"],
        geno_error_offset = popscle_extra_dict["geno_error_offset"],
        geno_error_coeff = popscle_extra_dict["geno_error_coeff"],
        r2_info = popscle_extra_dict["r2_info"],
        min_mac = popscle_extra_dict["min_mac"],
        min_callrate = popscle_extra_dict["min_callrate"],
        doublet_prior = popscle_extra_dict["doublet_prior"],
        cap_BQ = popscle_extra_dict["cap_BQ"],
        min_BQ = popscle_extra_dict["min_BQ"],
        min_MQ = popscle_extra_dict["min_MQ"],
        min_TD = popscle_extra_dict["min_TD"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"]
    log: output_dict["output_dir"] + "/logs/popscle_demuxlet.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle demuxlet \
            --plp {params.pileup} \
            --vcf {input.snps} \
            --field {params.field} \
            --geno-error-offset {params.geno_error_offset} \
            --geno-error-coeff {params.geno_error_coeff} \
            --r2-info {params.r2_info} \
            --min-mac {params.min_mac} \
            --min-callrate {params.min_callrate} \
            --group-list {input.barcodes} \
            --sm-list {input.individuals} \
            --doublet-prior {params.doublet_prior} \
            --cap-BQ {params.cap_BQ} \
            --min-BQ {params.min_BQ} \
            --min-MQ {params.min_MQ} \
            --min-TD {params.min_TD} \
            --excl-flag {params.excl_flag} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}demuxletOUT 
        [[ -s {output} ]]
        echo $?
        """

###################################################
############ REFORMAT DEMUXLET RESULTS ############
###################################################
rule demuxlet_results_temp:
    input:
        demuxlet = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/demuxlet_results.txt"
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"],
        bind = bind_path
    log: output_dict["output_dir"] + "/logs/demuxlet_results_temp.{pool}.log"
    shell:
        """
            singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$13,$14,$19,$20}}' {input.demuxlet} | \
            singularity exec --bind {params.bind} {params.sif} sed "s/SNG/singlet/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/DBL/doublet/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/AMB/unassigned/g" | \
            singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}} $3=="doublet" {{$4="doublet"}}1' | \
            singularity exec --bind {params.bind} {params.sif} sed "s/NUM.SNPS/nSNP/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/DROPLET.TYPE/DropletType/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/singlet.BEST.GUESS/Assignment/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/singlet.BEST.LLK/SingletLLK/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/doublet.BEST.LLK/DoulbetLLK/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "1s/\t/\tdemuxlet_/g" | \
            singularity exec --bind {params.bind} {params.sif} sed "s/BARCODE/Barcode/g" | \
            singularity exec --bind {params.bind} {params.sif} awk 'NR<2{{print $0;next}}{{print $0 | "sort -k1"}}'  > {output}
        """

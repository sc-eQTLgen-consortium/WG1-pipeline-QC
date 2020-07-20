#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["ref_dir"]
popscle_dict = config["popscle"]
popscle_extra_dict = config["popscle_extra"]


###################################
############# POPSCLE #############
###################################
###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        vcf = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: read_barcodes(wildcards.pool, input_dict = scrnaseq_libs),
        bam = 
        individuals = 
    output:
        directory(output_dict["output_dir"] + "/{pool}/popscle/pileup/")
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"]
    threads: popscle_dict["pileup_threads"]
    params:
        sif=input_dict["singularity_image"],
        tag_group = popscle_dict["tag_group"],
        tag_UMI = popscle_dict["tag_UMI"],
        cap_BQ = popscle_extra_dict["cap_BQ"],
        min_BQ = popscle_extra_dict["min_BQ"],
        min_MQ = popscle_extra_dict["min_MQ"],
        min_TD = popscle_extra_dict["min_TD"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_uniq = popscle_extra_dict["min_uniq"],
        min_snp = popscle_extra_dict["min_snp"]
    shell:
        """
        singularity exec {params.sif} popscle dsc-pileup \
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
            --min-uniq {params.min_uniq} \
            --min-snp {params.min_snp} \
            --out {output}pileup
        [[ -s {output}pileup.var.gz ]]
        echo $?
        """

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup = output_dict["output_dir"] + "/{pool}/popscle/pileup/",
        snps = output_dict["snp_genotypes_filepath"],
        barcodes = ,
        individuals =
    output:
        output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_threads"]
    threads: popscle_dict["demuxlet_threads"]
    params:
        out = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/",
        sif = input_dict["singularity_image"],
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
        min_uniq = popscle_extra_dict["min_uniq"],
        min_snp = popscle_extra_dict["min_snp"]
    shell:
        """
        singularity exec {params.sif} popscle demuxlet \
            --plp {input.pileup}pileup \
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
            --min-uniq {params.min_uniq} \
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
        demuxlet = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT_impute_vars.best"
    output:
        demuxlet_temp = temp(output_dict["output_dir"] + "/{pool}/CombinedResults/demuxlet_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$6,$14,$19,$20}}' {input.demuxlet} | \
            singularity exec {params.sif} sed "s/SNG/singlet/g" | \
            singularity exec {params.sif} sed "s/DBL/doublet/g" | \
            singularity exec {params.sif} awk 'BEGIN{{FS=OFS="\t"}} $3=="doublet" {{$4="doublet"}}1' | \
            singularity exec {params.sif} sed -E "s/,[0-9]+_[0-9]+,[0-9].[0-9]+\t/\t/g" | sed "s/NUM.SNPS/nSNP/g" | \
            singularity exec {params.sif} sed "s/DROPLET.TYPE/DropletType/g" | \
            singularity exec {params.sif} sed "s/BEST.GUESS/Assignment/g" | \
            singularity exec {params.sif} sed "s/singlet.BEST.LLK/SingletLLK/g" | \
            singularity exec {params.sif} sed "s/doublet.BEST.LLK/DoulbetLLK/g" | \
            singularity exec {params.sif} sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | \
            singularity exec {params.sif} sed "1s/\t/\tdemuxlet_/g" | \
            singularity exec {params.sif} sed "s/BARCODE/Barcode/g" | \
            singularity exec {params.sif} awk 'NR<2{{print $0;next}}{{print $0 | "sort -k1"}}'  > {output.demuxlet_temp}
        """


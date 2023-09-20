#!/usr/bin/env python

###################################
############# POPSCLE #############
###################################
rule popscle_bam_filter:
    input:
        barcodes = SAMPLES.loc["{pool}", "Barcodes"],
        vcf = config["inputs"]["vcf"],
        bam = SAMPLES.loc["{pool}", "Bam"]
    output:
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"]
    threads: config["popscle"]["popscle_bam_filter_threads"]
    params:
        out_dir = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        tag_group = config["popscle"]["popscle_tag_group"]
    log: config["outputs"]["output_dir"] + "logs/popscle_bam_filter.{pool}.log"
    shell:
        """
        mkdir -p {params.out_dir} && \
            singularity exec --bind {params.bind} {params.sif} bedtools merge \
                -i {input.vcf} | singularity exec \
                    --bind {params.bind} {params.sif} samtools view \
                        --target-file \
                        - \
                        --tag-file {params.tag_group}:<(zcat {input.barcodes}) \
                        --output {output.bam} \
                        --write-index \
                        --threads {threads} \
                        {input.bam}
        """


rule popscle_pileup:
    input:
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam",
        sm_list = config["inputs"]["individual_list_dir"] + "{pool}.txt",
        vcf = config["inputs"]["vcf"],
        barcodes = SAMPLES.loc["{pool}", "Barcodes"],
    output:
        pileup = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup.var.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"]
    threads: config["popscle"]["popscle_pileup_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        tag_group = config["popscle"]["popscle_tag_group"],
        tag_UMI = config["popscle"]["popscle_tag_UMI"],
        cap_bq = config["popscle_extra"]["cap_bq"],
        min_bq = config["popscle_extra"]["min_bq"],
        min_mq = config["popscle_extra"]["min_mq"],
        min_td = config["popscle_extra"]["min_td"],
        excl_flag = config["popscle_extra"]["excl_flag"],
        min_total = config["popscle_extra"]["min_total"],
        min_snp = config["popscle_extra"]["min_snp"],
        out = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup",
    log: config["outputs"]["output_dir"] + "logs/popscle_pileup.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle dsc-pileup \
            --sam {input.bam} \
            --tag-group {params.tag_group} \
            --tag-UMI {params.tag_UMI} \
            --sm-list {input.sm_list} \
            --vcf {input.vcf} \
            --cap-BQ {params.cap_bq} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --min-TD {params.min_td} \
            --excl-flag {params.excl_flag} \
            --group-list {input.barcodes} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}
        """


rule popscle_demuxlet:
    input:
        pileup = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup.var.gz",
        vcf = config["inputs"]["vcf"],
        group_list = SAMPLES.loc["{pool}", "Barcodes"],
        sm_list = config["inputs"]["individual_list_dir"] + "{pool}.txt"
    output:
        out = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"]
    threads: config["popscle"]["popscle_demuxlet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        plp = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup",
        field = config["popscle"]["popscle_genotype_field"],
        geno_error_offset = config["popscle_extra"]["geno_error_offset"],
        geno_error_coeff = config["popscle_extra"]["geno_error_coeff"],
        r2_info = config["popscle_extra"]["r2_info"],
        min_mac = config["popscle_extra"]["min_mac"],
        min_callrate = config["popscle_extra"]["min_callrate"],
        doublet_prior = config["popscle_extra"]["doublet_prior"],
        cap_bq = config["popscle_extra"]["cap_bq"],
        min_bq = config["popscle_extra"]["min_bq"],
        min_mq = config["popscle_extra"]["min_mq"],
        min_td = config["popscle_extra"]["min_td"],
        excl_flag = config["popscle_extra"]["excl_flag"],
        min_total = config["popscle_extra"]["min_total"],
        min_snp = config["popscle_extra"]["min_snp"],
        out = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT"
    log: config["outputs"]["output_dir"] + "logs/popscle_demuxlet.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle demuxlet \
            --plp {params.plp} \
            --vcf {input.vcf} \
            --field {params.field} \
            --geno-error-offset {params.geno_error_offset} \
            --geno-error-coeff {params.geno_error_coeff} \
            --r2-info {params.r2_info} \
            --min-mac {params.min_mac} \
            --min-callrate {params.min_callrate} \
            --group-list {input.group_list} \
            --sm-list {input.sm_list} \
            --doublet-prior {params.doublet_prior} \
            --cap-BQ {params.cap_bq} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --min-TD {params.min_td} \
            --excl-flag {params.excl_flag} \
            --min-total {params.min_total} \
            --min-snp {params.min_snp} \
            --out {params.out}
        [[ -s {output.out} ]]
        echo $?
        """

####################################
############ SOUPORCELL ############
####################################
rule souporcell:
    input:
        bam = SAMPLES.loc["{pool}", "Bam"],
        barcodes = SAMPLES.loc["{pool}", "Barcodes"],
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        vcf = config["inputs"]["vcf"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_memory"]
    threads: config["souporcell"]["souporcell_threads"]
    output:
        ambient_rna = config["outputs"]["output_dir"] + "{pool}/souporcell/ambient_rna.txt",
        clustering = config["outputs"]["output_dir"] + "{pool}/souporcell/clustering.done",
        genotypes = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf",
        clusters = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv",
        concensus = config["outputs"]["output_dir"] + "{pool}/souporcell/consensus.done",
        troublet = config["outputs"]["output_dir"] + "{pool}/souporcell/troublet.done"
    params:
        out = config["outputs"]["output_dir"] + "{pool}/souporcell/",
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        clusters = lambda wildcards: SAMPLES_DICT[wildcards.pool],
        min_alt = config["souporcell_extra"]["min_alt"],
        min_ref = config["souporcell_extra"]["min_ref"],
        max_loci = config["souporcell_extra"]["max_loci"]
    log: config["outputs"]["output_dir"] + "logs/souporcell.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} souporcell_pipeline.py \
            --bam {input.bam} \
            --barcodes {input.barcodes} \
            --fasta {input.fasta} \
            --threads {threads} \
            --out_dir {params.out} \
            --clusters {params.clusters} \
            --common_variants {input.vcf} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} 2> {log}
        [[ -s {output.genotypes} ]]
        echo $?
        """


rule souporcell_summary:
    input:
        clusters = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv"
    output:
        summary = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_summary_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_summary_memory"]
    threads: config["souporcell"]["souporcell_summary_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "logs/souporcell_summary.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} \
            awk 'BEGIN{FS=OFS="\t"} $2=="unassigned" {$3="unassigned"}1' $1 \
            | awk 'BEGIN{FS=OFS="\t"}{print $3}' \
            | sed -E 's|[0-9]+/[0-9]+|doublet|g' \
            | tail -n+2 \
            | sort \
            | uniq -c \
            | sed -E 's/^ +//g' \
            | sed 's/ /\t/g' \
            | sed '1 i\Assignment N\tClassification' \
            | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}' {input.clusters} > {output.summary}
        """


rule souporcell_pool_vcf:
    input:
        vcf = config["inputs"]["vcf"],
        cluster_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf"
    output:
        filtered_refs_temp = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf",
        filtered_refs = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"]
    threads: config["souporcell"]["souporcell_pool_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        individuals = config["inputs"]["individual_list_dir"] + "/{pool}.txt"
    log: config["outputs"]["output_dir"] + "logs/souporcell_pool_vcf.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools intersect -a {input.vcf} -b {input.cluster_vcf} -f 1.0 -r -wa -header > {output.filtered_refs_temp} 2> {log}
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {params.individuals} -Oz -o {output.filtered_refs} {output.filtered_refs_temp} 2>> {log}
        """


rule souporcell_correlate_genotypes:
    input:
        reference_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        cluster_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf",
    output:
        correlation_file = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlations.tsv",
        correlation_img = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlation.png",
        assignments = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt"
    resources:
        mem_per_thread_gb = config["souporcell"]["souporcell_correlations_memory"],
        disk_per_thread_gb = config["souporcell"]["souporcell_correlations_memory"]
    threads: config["souporcell"]["souporcell_correlations_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        out = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations"
    log: config["outputs"]["output_dir"] + "logs/souporcell_correlate_genotypes.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script}
            --reference_vcf {input.reference_vcf} \
            --cluster_vcf {input.cluster_vcf} \
            --out {params.out} \
            2> {log}
        [[ -s {output.assignments} ]]
        echo $?
        """
#!/usr/bin/env python

#########################################
############# PREPROCESSING #############
#########################################


# In case of multiple inputs
rule combine_vcfs_all:
    input:
        vcfs = config["inputs"]["vcf"],
        indices = [vcf + ".csi" for vcf in config["inputs"]["vcf"]]
    output:
        vcf = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["combine_vcfs_all_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["combine_vcfs_all_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["combine_vcfs_all_time"]]
    threads: config["demultiplex_preprocessing"]["combine_vcfs_all_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
    log: config["outputs"]["output_dir"] + "log/combine_vcfs_all.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools merge -Oz {input.vcfs} > {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

# Add all the info fields
# Filter the Imputed SNP Genotype by Minor Allele Frequency (MAF) and INFO scores
rule filter4demultiplexing:
    input:
        vcf = config["inputs"]["vcf"][0] if len(config["inputs"]["vcf"]) == 1 else config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38.vcf.gz",
        index = config["inputs"]["vcf"][0] + ".csi" if len(config["inputs"]["vcf"]) == 1 else config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38.vcf.gz.csi"
    output:
        info_filled = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_info_filled.vcf.gz",
        qc_filtered = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered.vcf.gz",
        location_filtered = temp(config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered_exons.recode.vcf"),
        complete_cases = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered_exons_complete_cases.recode.vcf"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["filter4demultiplexing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["filter4demultiplexing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["filter4demultiplexing_time"]]
    threads: config["demultiplex_preprocessing"]["filter4demultiplexing_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        maf = config["demultiplex_preprocessing_extra"]["filter4demultiplexing_maf"],
        r2 = config["demultiplex_preprocessing_extra"]["filter4demultiplexing_r2"],
        bed = "/opt/hg38exonsUCSC.bed", # TODO: this file is removed from the image, refactor as input file
        out = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered_exons",
        complete_out = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered_exons_complete_cases",
    log: config["outputs"]["output_dir"] + "log/filter4demultiplexing.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fill-tags -Oz --output {output.info_filled} {input}
        singularity exec --bind {params.bind} {params.sif} bcftools filter --include 'MAF>={params.maf} & R2>={params.r2}' -Oz --output {output.qc_filtered} {output.info_filled}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --gzvcf {output.qc_filtered} \
            --max-alleles 2 \
            --remove-indels \
            --bed {params.bed} \
            --recode \
            --recode-INFO-all \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --recode \
            --recode-INFO-all \
            --vcf {output.location_filtered} \
            --max-missing 1 \
            --out {params.complete_out}
        """


rule sort4demultiplexing:
    input:
        complete_cases = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered_exons_complete_cases.recode.vcf"
    output:
        complete_cases_sorted = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["sort4demultiplexing_time"]]
    threads: config["demultiplex_preprocessing"]["sort4demultiplexing_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        jar = "/opt/picard-3.1.0/build/libs/picard.jar"
    log: config["outputs"]["output_dir"] + "log/sort4demultiplexing.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -Xms{resources.java_mem}g -jar {params.jar} SortVcf \
            I={input.complete_cases} \
            O={output.complete_cases_sorted}
        """


rule count_snps:
    input:
        info_filled = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_info_filled.vcf.gz",
        qc_filtered = config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38_qc_filtered.vcf.gz",
        complete_cases_sorted = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf"
    output:
        number_snps = report(config["outputs"]["output_dir"] + "count_snps/Number_SNPs.png", category = "SNP Numbers", caption = "../report_captions/counts_snps.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["count_snps_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["count_snps_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["count_snps_time"]]
    threads: config["demultiplex_preprocessing"]["count_snps_threads"]
    params:
        bind=config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script="/opt/WG1-pipeline-QC/Demultiplexing/scripts/SNP_numbers.R",
        basedir = config["outputs"]["output_dir"],
        outdir = config["outputs"]["output_dir"] + "metrics/",
    log: config["outputs"]["output_dir"] + "log/count_snps.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --indir {params.basedir} \
            --out {params.outdir}
        """



###################################
############# POPSCLE #############
###################################
rule popscle_bam_filter:
    input:
        barcodes = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"],
        vcf = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf",
        bam = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Bam"]
    output:
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/snpfiltered_alignment.bam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_bam_filter_time"]]
    threads: config["popscle"]["popscle_bam_filter_threads"]
    params:
        out_dir = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        tag_group = config["popscle"]["popscle_tag_group"]
    log: config["outputs"]["output_dir"] + "log/popscle_bam_filter.{pool}.log"
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
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/snpfiltered_alignment.bam",
        sm_list = config["inputs"]["individual_list_dir"] + "{pool}.txt",
        vcf = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf",
        barcodes = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"],
    output:
        pileup = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup.var.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_pileup_time"]]
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
    log: config["outputs"]["output_dir"] + "log/popscle_pileup.{pool}.log"
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
        vcf = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf",
        group_list = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"],
        sm_list = config["inputs"]["individual_list_dir"] + "{pool}.txt"
    output:
        out = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_demuxlet_time"]]
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
    log: config["outputs"]["output_dir"] + "log/popscle_demuxlet.{pool}.log"
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
        bam = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Bam"],
        barcodes = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "Barcodes"],
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        vcf = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf"
    output:
        ambient_rna = config["outputs"]["output_dir"] + "{pool}/souporcell/ambient_rna.txt",
        clustering = config["outputs"]["output_dir"] + "{pool}/souporcell/clustering.done",
        genotypes = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf",
        clusters = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv",
        concensus = config["outputs"]["output_dir"] + "{pool}/souporcell/consensus.done",
        troublet = config["outputs"]["output_dir"] + "{pool}/souporcell/troublet.done"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_time"]]
    threads: config["souporcell"]["souporcell_threads"]
    params:
        out = config["outputs"]["output_dir"] + "{pool}/souporcell/",
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        clusters = lambda wildcards: SAMPLES_DF.loc[wildcards.pool, "N"],
        min_alt = config["souporcell_extra"]["min_alt"],
        min_ref = config["souporcell_extra"]["min_ref"],
        max_loci = config["souporcell_extra"]["max_loci"]
    log: config["outputs"]["output_dir"] + "log/souporcell.{pool}.log"
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_summary_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_summary_time"]]
    threads: config["souporcell"]["souporcell_summary_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/souporcell_summary.{pool}.log"
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
        vcf = config["outputs"]["output_dir"] + "vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf",
        cluster_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf"
    output:
        filtered_refs_temp = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf",
        filtered_refs = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_pool_vcf_time"]]
    threads: config["souporcell"]["souporcell_pool_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        individuals = config["inputs"]["individual_list_dir"] + "/{pool}.txt"
    log: config["outputs"]["output_dir"] + "log/souporcell_pool_vcf.{pool}.log"
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
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_correlations_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_correlations_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_correlations_time"]]
    threads: config["souporcell"]["souporcell_correlations_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        out = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations"
    log: config["outputs"]["output_dir"] + "log/souporcell_correlate_genotypes.{pool}.log"
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
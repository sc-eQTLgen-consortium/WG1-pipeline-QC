##############################################
############# POPSCLE BAM FILTER #############
##############################################
rule popscle_bam_filter:
    input:
        barcodes = lambda wildcards: scrnaseq_libs_df["BarcodeFile"][wildcards.pool],
        vcf = input_dict["vcf"],
        bam = lambda wildcards: scrnaseq_libs_df["BamFile"][wildcards.pool],
    output:
        bam = output_dict["output_dir"] + "/{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_bam_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_bam_filter_memory"]
    threads: popscle_dict["popscle_bam_filter_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"]
    log: output_dict["output_dir"] + "/logs/popscle_bam_filter.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools merge \
            -i {input.vcf} \
                singularity exec --bind {params.bind} {params.sif} samtools view \
                    --target-file \
                    - \
                    --tag-file CB:<(zcat {input.barcodes}) \
                    --output {output.bam} \
                    --write-index \
                    --threads {threads} \
                    {input.bam}
        """


########################################
############ POPSCLE PILEUP ############
########################################
rule popscle_pileup:
    input:
        bam = output_dict["output_dir"] + "/{pool}/popscle/bam_filter/{pool}_snpfiltered_alignment.bam",
        individuals = lambda wildcards: scrnaseq_libs_df["IndividualFile"][wildcards.pool],
        vcf = input_dict["vcf"],
        barcodes = lambda wildcards: scrnaseq_libs_df["BarcodeFile"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_pileup_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_pileup_memory"]
    threads: popscle_dict["popscle_pileup_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        tag_group = popscle_dict["popscle_tag_group"],
        tag_UMI = popscle_dict["popscle_tag_UMI"],
        cap_bq = popscle_extra_dict["cap_bq"],
        min_bq = popscle_extra_dict["min_bq"],
        min_mq = popscle_extra_dict["min_mq"],
        min_td = popscle_extra_dict["min_td"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"],
        out = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup"
    log: output_dict["output_dir"] + "/logs/popscle_pileup.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle dsc-pileup \
            --sam {input.bam} \
            --tag-group {params.tag_group} \
            --tag-UMI {params.tag_UMI} \
            --sm-list {input.individuals} \
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

##########################################
############ POPSCLE DEMUXLET ############
##########################################
rule popscle_demuxlet:
    input:
        pileup = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz",
        vcf = input_dict["vcf"],
        group_list = lambda wildcards: scrnaseq_libs_df["BarcodeFile"][wildcards.pool],
        sm_list = lambda wildcards: scrnaseq_libs_df["IndividualFile"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["popscle_demuxlet_memory"]
    threads: popscle_dict["popscle_demuxlet_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        plp = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup",
        field = popscle_dict["popscle_genotype_field"],
        geno_error_offset = popscle_extra_dict["geno_error_offset"],
        geno_error_coeff = popscle_extra_dict["geno_error_coeff"],
        r2_info = popscle_extra_dict["r2_info"],
        min_mac = popscle_extra_dict["min_mac"],
        min_callrate = popscle_extra_dict["min_callrate"],
        doublet_prior = popscle_extra_dict["doublet_prior"],
        cap_bq = popscle_extra_dict["cap_bq"],
        min_bq = popscle_extra_dict["min_bq"],
        min_mq = popscle_extra_dict["min_mq"],
        min_td = popscle_extra_dict["min_td"],
        excl_flag = popscle_extra_dict["excl_flag"],
        min_total = popscle_extra_dict["min_total"],
        min_snp = popscle_extra_dict["min_snp"],
        out = output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT"
    log: output_dict["output_dir"] + "/logs/popscle_demuxlet.{pool}.log"
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
        [[ -s {output} ]]
        echo $?
        """
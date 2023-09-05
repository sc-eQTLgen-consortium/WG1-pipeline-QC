###################################
############# POPSCLE #############
###################################
###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool],
        vcf = input_dict["snp_genotypes_filepath"],
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * popscle_dict["pileup_memory"]
    threads: popscle_dict["pileup_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        tag_group = popscle_dict["tag_group"],
        tag_UMI = popscle_dict["tag_UMI"],
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

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup.var.gz",
        vcf = input_dict["snp_genotypes_filepath"],
        group_list = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool],
        sm_list = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    output:
        output_dict["output_dir"] + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * popscle_dict["demuxlet_threads"]
    threads: popscle_dict["demuxlet_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        plp = output_dict["output_dir"] + "/{pool}/popscle/pileup/pileup",
        field = popscle_dict["genotype_field"],
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
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"]
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
#!/usr/bin/env python

# Input: WGS genotype data in VCF or VCF per chromosome
# Output: filtered PLINK binary with --max-alleles 2, --mind and --set-all-var-ids applied

######################################################
############ IF THE INPUT IS A SINGLE VCF ############
######################################################

rule split_input_vcf_by_chr:
    input:
        vcf = config["inputs"]["genotype_path"] + ".vcf.gz",
        index = config["inputs"]["genotype_path"] + ".vcf.gz.tbi"
    output:
        vcf = temp(config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz"),
        index = temp(config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz.csi")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["split_by_chr_time"]]
    threads: config["generic"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/split_input_vcf_by_chr.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


#####################################################


rule bcftools_norm:
    input:
        vcf = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz" if INPUT_OPTION == "VCF per chromosome" else config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz",
        index = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz.tbi" if INPUT_OPTION == "VCF per chromosome" else config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz.tbi",
    output:
        vcf = temp(config["outputs"]["output_dir"] + "bcftools_norm_by_chr/chr_{chr}_normalised.vcf.gz")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["bcftools_norm_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["bcftools_norm_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["wgs_filter"]["bcftools_norm_time"]]
    threads: config["wgs_filter"]["bcftools_norm_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/bcftools_norm.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools norm -m -any {input.vcf} -o {output.vcf}
        """


#########################################################
############ FILTER AND CONCAT THE WGS INPUT ############
#########################################################


rule wgs_filter:
    input:
        vcf = config["outputs"]["output_dir"] + "bcftools_norm_by_chr/chr_{chr}_normalised.vcf.gz"
    output:
        vcf = temp(config["outputs"]["output_dir"] + "wgs_filter_by_chr/chr_{chr}_normalised_filtered.vcf.gz"),
        log = config["outputs"]["output_dir"] + "wgs_filter_by_chr/chr_{chr}_normalised_filtered.log.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["wgs_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["wgs_filter_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["wgs_filter"]["wgs_filter_time"]]
    threads: config["wgs_filter"]["wgs_filter_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Imputation/scripts/wgs_vcf_filter.py",
        out_dir = config["outputs"]["output_dir"] + "wgs_filter",
        sex = lambda wildcards: "--sex {}".format(config["inputs"]["psam"]) if wildcards.chr in ["X", "Y"] else "",
        genotype_quality = config["wgs_filter_extra"]["genotype_quality"],
        allelic_balance_lower = config["wgs_filter_extra"]["allelic_balance_lower"],
        allelic_balance_upper = config["wgs_filter_extra"]["allelic_balance_upper"],
        inbreeding_coefficient = config["wgs_filter_extra"]["inbreeding_coefficient"],
        no_snv_vqsr_check = "--no_snv_vqsr_check" if config["wgs_filter_extra"]["no_snv_vqsr_check"] else "",
        vqsr_snv = config["wgs_filter_extra"]["vqsr_snv"],
        no_indel_vqsr_check = "--no_indel_vqsr_check" if config["wgs_filter_extra"]["no_indel_vqsr_check"] else "",
        vqsr_indel = config["wgs_filter_extra"]["vqsr_indel"],
        keep_multialleic = "--keep_multialleic" if config["wgs_filter_extra"]["keep_multialleic"] else "",
        keep_non_pass_snv = "--keep_non_pass_snv" if config["wgs_filter_extra"]["keep_non_pass_snv"] else "",
        keep_non_pass_indel = "--keep_non_pass_indel" if config["wgs_filter_extra"]["keep_non_pass_indel"] else "",
        keep_low_complexity = "--keep_low_complexity" if config["wgs_filter_extra"]["keep_low_complexity"] else "",
        minor_allele_frequency = config["wgs_filter_extra"]["minor_allele_frequency"],
        call_rate = config["wgs_filter_extra"]["call_rate"],
        hardy_weinberg_equilibrium = config["wgs_filter_extra"]["hardy_weinberg_equilibrium"],
        filtered_depth = config["wgs_filter_extra"]["filtered_depth"],
        keep_info_column = "--keep_info_column" if config["wgs_filter_extra"]["keep_info_column"]  else ""
    log: config["outputs"]["output_dir"] + "log/wgs_filter.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input {input.vcf} \
            --output {output.vcf} \
            --log {output.log} \
            {params.sex} \
            --genotype_quality {params.genotype_quality} \
            --allelic_balance_lower {params.allelic_balance_lower} \
            --allelic_balance_upper {params.allelic_balance_upper} \
            --inbreeding_coefficient {params.inbreeding_coefficient} \
            {params.no_snv_vqsr_check} \
            --vqsr_snv {params.vqsr_snv} \
            {params.no_indel_vqsr_check} \
            --vqsr_indel {params.vqsr_indel} \
            {params.keep_multialleic} \
            {params.keep_non_pass_snv} \
            {params.keep_non_pass_indel} \
            {params.keep_low_complexity} \
            --minor_allele_frequency {params.minor_allele_frequency} \
            --call_rate {params.call_rate} \
            --hardy_weinberg_equilibrium {params.hardy_weinberg_equilibrium} \
            --filtered_depth {params.filtered_depth} \
            {params.keep_info_column}
        """


rule wgs_filtered_vcf_to_pgen:
    input:
        vcf = config["outputs"]["output_dir"] + "wgs_filter_by_chr/chr_{chr}_normalised_filtered.vcf.gz"
    output:
        pgen = temp(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.pgen"),
        pvar = temp(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.pvar"),
        psam = temp(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.psam"),
        log = config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.log"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["vcf_to_plink_memory"] * config["generic"]["vcf_to_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["vcf_to_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["vcf_to_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["vcf_to_plink_time"]]
    threads: config["generic"]["vcf_to_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        psam = config["inputs"]["psam"],
        split_par_flag = lambda wildcards: "--split-par " + config["inputs"]["genome_build"] if wildcards.chr == "X" else "",
        max_allele_len = config["pre_processing_extra"]["max_allele_len"],
        out = config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered"
    log: config["outputs"]["output_dir"] + "log/wgs_filtered_vcf_to_pgen.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --vcf {input.vcf} \
            --psam {params.psam} \
            {params.split_par_flag} \
            --max-alleles 2 \
            --new-id-max-allele-len {params.max_allele_len} \
            --set-all-var-ids @:#:\$r_\$a \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """


rule wgs_filter_stats:
    input:
        vcf = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter_by_chr/chr_{chr}_normalised_filtered.log.gz", chr=INPUT_CHROMOSOMES)
    output:
        stats = config["outputs"]["output_dir"] + "wgs_filtered/wgs_filter_stats.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["wgs_filter_stats_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["wgs_filter"]["wgs_filter_stats_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["wgs_filter"]["wgs_filter_stats_time"]]
    threads: config["wgs_filter"]["wgs_filter_stats_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "Imputation/scripts/wgs_vcf_filter_stats.py",
        minor_allele_frequency = config["wgs_filter_extra"]["minor_allele_frequency"],
        call_rate = config["wgs_filter_extra"]["call_rate"],
        hardy_weinberg_equilibrium = config["wgs_filter_extra"]["hardy_weinberg_equilibrium"],
    log: config["outputs"]["output_dir"] + "log/wgs_filter_stats.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input {input.vcf} \
            --output {output.stats} \
            --minor_allele_frequency {params.minor_allele_frequency} \
            --call_rate {params.call_rate} \
            --hardy_weinberg_equilibrium {params.hardy_weinberg_equilibrium}
        """


# This also includes the indiv_missingness rule.
# IMPORTANT: --pmerge-list reorders to psam file so '--indiv-sort f' is very important here.
rule combine_wgs_filtered_pgens:
    input:
        pgen = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.pgen", chr=INPUT_CHROMOSOMES),
        pvar = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.pvar", chr=INPUT_CHROMOSOMES),
        psam = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered.psam", chr=INPUT_CHROMOSOMES),
    output:
        merge_pgen = temp(config["outputs"]["output_dir"] + "wgs_filtered/normalised_filtered.pgen"),
        merge_pvar = temp(config["outputs"]["output_dir"] + "wgs_filtered/normalised_filtered.pvar"),
        merge_psam = temp(config["outputs"]["output_dir"] + "wgs_filtered/normalised_filtered.psam"),
        pgen = config["outputs"]["output_dir"] + "wgs_filtered/data.pgen",
        pvar = config["outputs"]["output_dir"] + "wgs_filtered/data.pvar",
        psam = config["outputs"]["output_dir"] + "wgs_filtered/data.psam",
        log = config["outputs"]["output_dir"] + "wgs_filtered/data.log"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["process_plink_memory"] * config["generic"]["process_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["process_plink_time"]]
    threads: config["generic"]["process_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter_by_chr_pgen/chr_{chr}_normalised_filtered", chr=INPUT_CHROMOSOMES),
        pmerge_list = config["outputs"]["output_dir"] + "wgs_filtered/pmerge_list.txt",
        mind = config["pre_processing_extra"]["mind"],
        merge_out = config["outputs"]["output_dir"] + "wgs_filtered/normalised_filtered",
        psam = config["inputs"]["psam"],
        out = config["outputs"]["output_dir"] + "wgs_filtered/data"
    log: config["outputs"]["output_dir"] + "log/combine_wgs_filtered_pgens.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.infiles} | sed 's/ /\\n/g' > {params.pmerge_list}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pmerge-list {params.pmerge_list} pfile \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.merge_out}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {output.merge_pgen} \
            --pvar {output.merge_pvar} \
            --psam {output.merge_psam} \
            --indiv-sort f {params.psam} \
            --mind {params.mind} \
            --rm-dup 'force-first' \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

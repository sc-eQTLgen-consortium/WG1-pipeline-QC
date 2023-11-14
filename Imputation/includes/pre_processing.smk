#!/usr/bin/env python

# Input: microarray genotype data in VCF /PLINK binary or VCF / PLINK binary per chromosome
# Output: PLINK binary with --max-alleles 2, --mind and --set-all-var-ids applied

###############################################################
############ PREPARE FOR QC: OPTION 1 = SINGLE VCF ############
###############################################################

# This also includes the indiv_missingness rule.
rule input_vcf_to_pgen:
    input:
        vcf = config["inputs"]["genotype_path"] + ".vcf.gz",
        index = config["inputs"]["genotype_path"] + ".vcf.gz.tbi"
    output:
        pgen = config["outputs"]["output_dir"] + "pre_processed/data.pgen",
        pvar = config["outputs"]["output_dir"] + "pre_processed/data.pvar",
        psam = config["outputs"]["output_dir"] + "pre_processed/data.psam",
        log = config["outputs"]["output_dir"] + "pre_processed/data.log",
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
        genome_build = config["inputs"]["genome_build"],
        max_allele_len = config["pre_processing_extra"]["max_allele_len"],
        mind = config["pre_processing_extra"]["mind"],
        chr = INPUT_CHROMOSOMES,
        out = config["outputs"]["output_dir"] + "pre_processed/data"
    log: config["outputs"]["output_dir"] + "log/input_vcf_to_pgen.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --vcf {input.vcf} \
            --psam {params.psam} \
            --split-par {params.genome_build} \
            --max-alleles 2 \
            --new-id-max-allele-len {params.max_allele_len} \
            --mind {params.mind} \
            --chr {params.chr} \
            --rm-dup 'force-first' \
            --set-all-var-ids @:#:\$r_\$a \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

#################################################################
############ PREPARE FOR QC: OPTION 2 = SINGLE PLINK ############
#################################################################

# This also includes the indiv_missingness rule.
# Note that the PSAM is already been validated for being in the same order as the 'config["inputs"]["genotype_path"].psam'
# in the SnakeFile preprocessing checks.
rule input_pgen:
    input:
        pgen = config["inputs"]["genotype_path"] + ".pgen",
        pvar = config["inputs"]["genotype_path"] + ".pvar",
        psam = config["inputs"]["genotype_path"] + ".pvar"
    output:
        pgen = config["outputs"]["output_dir"] + "pre_processed/data.pgen",
        pvar = config["outputs"]["output_dir"] + "pre_processed/data.pvar",
        psam = config["outputs"]["output_dir"] + "pre_processed/data.psam",
        log = config["outputs"]["output_dir"] + "pre_processed/data.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["process_plink_memory"] * config["generic"]["process_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards,attempt: attempt * config["generic"]["process_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["process_plink_time"]]
    threads: config["generic"]["process_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        max_allele_len = config["pre_processing_extra"]["max_allele_len"],
        out = config["outputs"]["output_dir"] + "pre_processed/data",
        mind = config["pre_processing_extra"]["mind"],
        chr = INPUT_CHROMOSOMES
    log: config["outputs"]["output_dir"] + "log/input_pgen.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --max-alleles 2 \
            --new-id-max-allele-len {params.max_allele_len} \
            --mind {params.mind} \
            --chr {params.chr} \
            --rm-dup 'force-first' \
            --set-all-var-ids @:#:\$r_\$a \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """


#######################################################################
############ PREPARE FOR QC: OPTION 3 = VCF PER CHROMOSOME ############
#######################################################################

rule input_vcf_to_pgen_per_chr:
    input:
        vcf = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz" if INPUT_OPTION == "VCF per chromosome" else "",
        index = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz.tbi" if INPUT_OPTION == "VCF per chromosome" else ""
    output:
        pgen = config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.pgen",
        pvar = config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.pvar",
        psam = config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.psam",
        log = config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.log"
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
        out = config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data"
    log: config["outputs"]["output_dir"] + "log/input_vcf_to_pgen_per_chr.chr_{chr}.log"
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

# This also includes the indiv_missingness rule.
# IMPORTANT: --pmerge-list reorders to psam file so '--indiv-sort f' is very important here.
rule merge_vcf_to_pgens:
    input:
        pgen = lambda wildcards: expand(config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.pgen", chr=INPUT_CHROMOSOMES),
        pvar = lambda wildcards: expand(config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.pvar", chr=INPUT_CHROMOSOMES),
        psam = lambda wildcards: expand(config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}.psam", chr=INPUT_CHROMOSOMES),
    output:
        merge_pgen = temp(config["outputs"]["output_dir"] + "pre_processed/merged.pgen"),
        merge_pvar = temp(config["outputs"]["output_dir"] + "pre_processed/merged.pvar"),
        merge_psam = temp(config["outputs"]["output_dir"] + "pre_processed/merged.psam"),
        pgen = config["outputs"]["output_dir"] + "pre_processed/data.pgen",
        pvar = config["outputs"]["output_dir"] + "pre_processed/data.pvar",
        psam = config["outputs"]["output_dir"] + "pre_processed/data.psam",
        log = config["outputs"]["output_dir"] + "pre_processed/data.log"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["process_plink_memory"] * config["generic"]["process_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["process_plink_time"]]
    threads: config["generic"]["process_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["outputs"]["output_dir"] + "input_vcf_to_pgen_per_chr/data_{chr}", chr=INPUT_CHROMOSOMES),
        pmerge_list = config["outputs"]["output_dir"] + "pre_processed/pmerge_list.txt",
        mind = config["pre_processing_extra"]["mind"],
        merge_out = config["outputs"]["output_dir"] + "pre_processed/merged",
        psam = config["inputs"]["psam"],
        out = config["outputs"]["output_dir"] + "pre_processed/data"
    log: config["outputs"]["output_dir"] + "log/merge_vcf_to_pgens.log"
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


#########################################################################
############ PREPARE FOR QC: OPTION 4 = PLINK PER CHROMOSOME ############
#########################################################################

# This also includes the indiv_missingness rule.
# IMPORTANT: --pmerge-list reorders to psam file so '--indiv-sort f' is very important here.
rule merge_input_pgens:
    input:
        pgen = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".pgen", chr=INPUT_CHROMOSOMES) if INPUT_OPTION == "PLINK per chromosome" else "",
        pvar = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".pvar", chr=INPUT_CHROMOSOMES) if INPUT_OPTION == "PLINK per chromosome" else "",
        psam = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".psam", chr=INPUT_CHROMOSOMES) if INPUT_OPTION == "PLINK per chromosome" else "",
    output:
        merge_pgen = temp(config["outputs"]["output_dir"] + "pre_processed/merged.pgen"),
        merge_pvar = temp(config["outputs"]["output_dir"] + "pre_processed/merged.pvar"),
        merge_psam = temp(config["outputs"]["output_dir"] + "pre_processed/merged.psam"),
        pgen = config["outputs"]["output_dir"] + "pre_processed/data.pgen",
        pvar = config["outputs"]["output_dir"] + "pre_processed/data.pvar",
        psam = config["outputs"]["output_dir"] + "pre_processed/data.psam",
        log = config["outputs"]["output_dir"] + "pre_processed/data.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["process_plink_memory"] * config["generic"]["process_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["process_plink_time"]]
    threads: config["generic"]["process_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", wildcards.chr), chr=INPUT_CHROMOSOMES),
        pmerge_list = config["outputs"]["output_dir"] + "pre_processed/pmerge_list.txt",
        max_allele_len = config["pre_processing_extra"]["max_allele_len"],
        mind = config["pre_processing_extra"]["mind"],
        merge_out = config["outputs"]["output_dir"] + "pre_processed/merged",
        psam = config["inputs"]["psam"],
        out = config["outputs"]["output_dir"] + "pre_processed/data",
    log: config["outputs"]["output_dir"] + "log/merge_input_pgens.log"
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
            --max-alleles 2 \
            --new-id-max-allele-len {params.max_allele_len} \
            --mind {params.mind} \
            --rm-dup 'force-first' \
            --set-all-var-ids @:#:\$r_\$a \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

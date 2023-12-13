#!/usr/bin/env python

# Input: PLINK binary on genome build 38 with updated sex, ancestry and split per ancestry
# Output: PLINK binary imputed (possibly per split per dataset)
# IMPORTANT: be careful here, we can no longer assume the same order PSAM as the one we initially validated in the Snakefile.

rule harmonize:
    input:
        bed = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bed",
        bim = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bim",
        fam = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize/{ancestry}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize/{ancestry}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize/{ancestry}.fam",
        log = config["outputs"]["output_dir"] + "harmonize/{ancestry}.log",
        id_updates = config["outputs"]["output_dir"] + "harmonize/{ancestry}_idUpdates.txt",
        snp_log = config["outputs"]["output_dir"] + "harmonize/{ancestry}_snpLog.log",
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"] * config["imputation"]["harmonize_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["harmonize_time"]]
    threads: config["imputation"]["harmonize_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset",
        out = config["outputs"]["output_dir"] + "harmonize/{ancestry}",
    log: config["outputs"]["output_dir"] + "log/harmonize.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem_gb}g -jar {params.jar}\
            --input {params.infile} \
            --inputType PLINK_BED \
            --ref {input.ref_vcf} \
            --refType VCF \
            --update-id \
            --output {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.log}
        singularity exec --bind {params.bind} {params.sif} touch {output.id_updates}
        singularity exec --bind {params.bind} {params.sif} touch {output.snp_log}
        """


rule harmonized_bed_to_vcf:
    input:
        bed = config["outputs"]["output_dir"] + "harmonize/{ancestry}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize/{ancestry}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize/{ancestry}.fam"
    output:
        vcf = config["outputs"]["output_dir"] + "harmonized_bed_to_vcf/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonized_bed_to_vcf/{ancestry}_harmonised_hg38.vcf.gz.csi",
        log = config["outputs"]["output_dir"] + "harmonized_bed_to_vcf/{ancestry}_harmonised_hg38.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["plink_to_vcf_memory"] * config["generic"]["plink_to_vcf_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["plink_to_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["plink_to_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["plink_to_vcf_time"]]
    threads: config["generic"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "harmonized_bed_to_vcf/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "log/harmonized_bed_to_vcf.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --recode vcf id-paste=iid \
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

####################################################################
############ WGS VCF ARE TOO BIG SO NEED TO BE SPLITTED ############
####################################################################

rule split_by_chr_for_harmonize:
    input:
        bed = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bed",
        bim = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bim",
        fam = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.fam"
    output:
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.fam",
        log = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["split_by_chr_memory"] * config["generic"]["split_by_chr_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["split_by_chr_time"]]
    threads: config["generic"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset"
    log: config["outputs"]["output_dir"] + "log/split_by_chr_for_harmonize.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --chr {wildcards.chr} \
            --make-bed \
            --out {params.out}
        """


rule harmonize_per_chr:
    input:
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.fam",
        log = config["outputs"]["output_dir"] + "harmonize/{ancestry}_chr_{chr}.log",
        id_updates = config["outputs"]["output_dir"] + "harmonize/{ancestry}_chr_{chr}_idUpdates.txt",
        snp_log = config["outputs"]["output_dir"] + "harmonize/{ancestry}_chr_{chr}_snpLog.log",
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"] * config["imputation"]["harmonize_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["harmonize_time"]]
    threads: config["imputation"]["harmonize_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_subset",
        out = config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}",
    log: config["outputs"]["output_dir"] + "log/harmonize_per_chr.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -jar {params.jar} \
            --input {params.infile} \
            --inputType PLINK_BED \
            --ref {input.ref_vcf} \
            --refType VCF \
            --update-id \
            --output {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.log}
        singularity exec --bind {params.bind} {params.sif} touch {output.id_updates}
        singularity exec --bind {params.bind} {params.sif} touch {output.snp_log}
        """


rule harmonized_bed_per_chr_to_vcf:
    input:
        bed = expand(config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.bed", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        bim = expand(config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.bim", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        fam = expand(config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}.fam", ancestry=ANCESTRIES, chr=CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "harmonized_bed_per_chr_to_vcf/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonized_bed_per_chr_to_vcf/{ancestry}_harmonised_hg38.vcf.gz.csi",
        log = config["outputs"]["output_dir"] + "harmonized_bed_per_chr_to_vcf/{ancestry}_harmonised_hg38.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["plink_to_vcf_memory"] * config["generic"]["plink_to_vcf_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["plink_to_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["plink_to_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["plink_to_vcf_time"]]
    threads: config["generic"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_per_chr/{ancestry}_chr_{chr}", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        mergelist = config["outputs"]["output_dir"] + "harmonized_bed_per_chr_to_vcf/{ancestry}_mergelist.txt",
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "harmonized_bed_per_chr_to_vcf/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "log/harmonized_bed_per_chr_to_vcf.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.infiles} | sed 's/ /\\n/g' > {params.mergelist}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --pmerge-list {params.mergelist} bfile \
            --recode vcf id-paste=iid \
            --chr {params.chr_list} \
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

####################################################################

rule fixref:
    input:
        ref_fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi",
        vcf = config["outputs"]["output_dir"] + ("harmonized_bed_per_chr_to_vcf/" if config["settings"]["is_wgs"] else "harmonized_bed_to_vcf/") + "{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + ("harmonized_bed_per_chr_to_vcf/" if config["settings"]["is_wgs"] else "harmonized_bed_to_vcf/") + "{ancestry}_harmonised_hg38.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "fixref/{ancestry}_fixref_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "fixref/{ancestry}_fixref_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["fixref_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["fixref_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["fixref_time"]]
    threads: config["imputation"]["fixref_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/fixref.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fixref {input.vcf} -- -f {input.ref_fasta} -i {input.ref_vcf} | \
            singularity exec --bind {params.bind} {params.sif} bcftools norm \
                --check-ref x \
                --fasta-ref {input.ref_fasta} \
                -Oz \
                --output {output.vcf} > {log}

        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = config["outputs"]["output_dir"] + "fixref/{ancestry}_fixref_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "fixref/{ancestry}_fixref_hg38.vcf.gz.csi"
    output:
        tagged_vcf = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["filter_preimpute_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["filter_preimpute_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["filter_preimpute_vcf_time"]]
    threads: config["imputation"]["filter_preimpute_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        hwe = config["imputation"]["filter_preimpute_vcf_snp_hwe"],
        missing = config["imputation"]["filter_preimpute_vcf_snp_missing_pct"],
        maf = lambda wildcards: ANCESTRY_MAF_DICT[wildcards.ancestry]
    log: config["outputs"]["output_dir"] + "log/filter_preimpute_vcf.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        singularity exec --bind {params.bind} {params.sif} bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} |\
            singularity exec --bind {params.bind} {params.sif} bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
                singularity exec --bind {params.bind} {params.sif} bcftools filter -e "ALT='.'" |\
                    singularity exec --bind {params.bind} {params.sif} bcftools norm -d all |\
                        singularity exec --bind {params.bind} {params.sif} bcftools norm -m+any |\
                            singularity exec --bind {params.bind} {params.sif} bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        singularity exec --bind {params.bind} {params.sif} bcftools index {output.filtered_vcf}
        """


# Note, downgrading VCF version in the header to prevent the following error:
# Error: VCF version must be v4.0, v4.1 or v4.2:
# You are using version VCFv4.3
rule het:
    input:
        vcf = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(config["outputs"]["output_dir"] + "het/{ancestry}_filtered_temp.vcf.gz"),
        het = config["outputs"]["output_dir"] + "het/{ancestry}_het.het",
        failed_inds = config["outputs"]["output_dir"] + "het/{ancestry}_het_failed.inds",
        passed_inds = config["outputs"]["output_dir"] + "het/{ancestry}_het_passed.inds",
        passed_list = config["outputs"]["output_dir"] + "het/{ancestry}_het_passed.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["het_time"]]
    threads: config["imputation"]["het_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "het/{ancestry}_het",
        script = config["inputs"]["repo_dir"] + "Imputation/report_captions/filter_het.R"
    log: config["outputs"]["output_dir"] + "log/het.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | gzip -c > {output.tmp_vcf}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --gzvcf {output.tmp_vcf} \
            --het \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --input {output.het} \
            --out {params.out}
        """


rule het_filter:
    input:
        passed_list = config["outputs"]["output_dir"] + "het/{ancestry}_het_passed.txt",
        vcf = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["het_filter_time"]]
    threads: config["imputation"]["het_filter_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/het_filter.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


# Note, downgrading VCF version in the header to prevent the following error:
# Error: VCF version must be v4.0, v4.1 or v4.2:
# You are using version VCFv4.3
rule calculate_missingness:
    input:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(config["outputs"]["output_dir"] + "calculate_missingness/{ancestry}_het_filtered_temp.vcf.gz"),
        miss = config["outputs"]["output_dir"] + "calculate_missingness/{ancestry}_genotypes.imiss"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["calculate_missingness_time"]]
    threads: config["imputation"]["calculate_missingness_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "calculate_missingness/{ancestry}_genotypes"
    log: config["outputs"]["output_dir"] + "log/calculate_missingness.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | gzip -c > {output.tmp_vcf}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --gzvcf {output.tmp_vcf} \
            --missing-indv \
            --out {params.out}
        """


rule kinship:
    input:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz"
    output:
        subset_pgen = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset.pgen",
        subset_pvar = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset.pvar",
        subset_psam = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset.psam",
        subset_log = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset.log",
        prune_in = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruning.prune.in",
        prune_out = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruning.prune.out",
        prune_log = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruning.log",
        pruned_log = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruned.log",
        king = config["outputs"]["output_dir" ] + "kinship/{ancestry}_subset_pruned.king",
        king_id = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruned.king.id",
        kinship = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruned.kinship"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["imputation"]["kinship_memory"] * config["imputation"]["kinship_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["kinship_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["kinship_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["kinship_time"]]
    threads: config["imputation"]["kinship_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        genome_build = config["inputs"]["genome_build"],
        maf = config["imputation_extra"]["kinship_maf"],
        hwe = config["imputation_extra"]["kinship_hwe"],
        out_subset = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset",
        indep_window_size = config["imputation_extra"]["kinship_indep_window_size"],
        indep_step_size = config["imputation_extra"]["kinship_indep_step_size"],
        pairwise_r2_threshold = config["imputation_extra"]["kinship_pairwise_r2_threshold"],
        out_pruning = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruning",
        out_pruned = config["outputs"]["output_dir"] + "kinship/{ancestry}_subset_pruned",
        script = config["inputs"]["repo_dir"] + "Imputation/report_captions/kinship.R",
    log: config["outputs"]["output_dir"] + "log/kinship.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --vcf {input.vcf} \
            --split-par {params.genome_build} \
            --maf {params.maf} \
            --hwe {params.hwe} \
            --make-pgen \
            --out {params.out_subset}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --pgen {output.subset_pgen} \
            --pvar {output.subset_pvar} \
            --psam {output.subset_psam} \
            --indep-pairwise {params.indep_window_size} {params.indep_step_size} {params.pairwise_r2_threshold} \
            --bad-ld \
            --out {params.out_pruning}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --pgen {output.subset_pgen} \
            --pvar {output.subset_pvar} \
            --psam {output.subset_psam} \
            --extract {output.prune_in} \
            --make-king square \
            --out {params.out_pruned}
            
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --king {output.king} \
            --king_id {output.king_id} \
            --out {output.kinship}
        """


rule split_by_chr_for_prephasing:
    input:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "split_by_chr_for_prephasing/{ancestry}_chr_{chr}.vcf.gz",
        index = config["outputs"]["output_dir"] + "split_by_chr_for_prephasing/{ancestry}_chr_{chr}.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["split_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["split_by_chr_time"]]
    threads: config["generic"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/split_by_chr_for_prephasing.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = config["outputs"]["output_dir"] + "split_by_chr_for_prephasing/{ancestry}_chr_{chr}.vcf.gz",
        map_file = config["refs"]["ref_dir"] + config["refs_extra"]["relative_map_path"],
        phasing_bcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "chr{chr}.bcf",
        phasing_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "chr{chr}.bcf.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "eagle_prephasing_by_chr/{ancestry}_chr_{chr}_phased.vcf.gz",
        index = config["outputs"]["output_dir"] + "eagle_prephasing_by_chr/{ancestry}_chr_{chr}_phased.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["eagle_prephasing_time"]]
    threads: config["imputation"]["eagle_prephasing_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "eagle_prephasing_by_chr/{ancestry}_chr_{chr}_phased"
    log: config["outputs"]["output_dir"] + "log/eagle_prephasing.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} eagle \
            --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_bcf} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads} > {log}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule merge_eagle_prephasing_stats:
    input:
        logs = expand(config["outputs"]["output_dir"] + "log/eagle_prephasing.{ancestry}.chr_{chr}.log", ancestry=ANCESTRIES, chr=CHROMOSOMES)
    output:
        stats = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_phase_confidence.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["merge_eagle_prephasing_stats_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["merge_eagle_prephasing_stats_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["merge_eagle_prephasing_stats_time"]]
    threads: config["imputation"]["merge_eagle_prephasing_stats_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        script = config["inputs"]["repo_dir"] + "Imputation/report_captions/merge_eagle_prephasing_stats.py",
    log: config["outputs"]["output_dir"] + "log/merge_eagle_prephasing_stats.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --logs {input.logs} \
            --out {output.stats}
        """


# Results are slightly different compared to v1.0.2
# ID column is rsID in NEW but chr:pos:ref:alt in OLD
# NEW INFO column contains AVG_CS + values are slightly different
# NEW FORMAT is reordered as GT:GP:DS
# temp-prefix parameter fixes issue I had running it on our cluster
rule minimac_imputation:
    input:
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["relative_imputation_dir"] + "chr{chr}.msav",
        target = config["outputs"]["output_dir"] + "eagle_prephasing_by_chr/{ancestry}_chr_{chr}_phased.vcf.gz",
        target_index = config["outputs"]["output_dir"] + "eagle_prephasing_by_chr/{ancestry}_chr_{chr}_phased.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "minimac_imputation_by_chr/{ancestry}_chr_{chr}.dose.vcf.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["minimac_imputation_time"]]
    threads: config["imputation"]["minimac_imputation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        temp_prefix = config["outputs"]["output_dir"] + "minimac_imputation_by_chr/{ancestry}_chr_{chr}_m4_",
        chunk = config["imputation"]["minimac_chunk"],
        overlap = config["imputation"]["minimac_overlap"]
    log: config["outputs"]["output_dir"] + "log/minimac_imputation.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} minimac4 \
            {input.reference} \
            {input.target} \
            --chunk {params.chunk} \
            --format GT,DS,GP \
            --temp-prefix {params.temp_prefix} \
            --output {output.vcf} \
            --output-format vcf.gz \
            --threads {threads} \
            --overlap {params.overlap}
        """


rule merge_vcfs:
    input:
        vcfs = expand(config["outputs"]["output_dir"] + "minimac_imputation_by_chr/{ancestry}_chr_{chr}.dose.vcf.gz", ancestry=ANCESTRIES, chr=CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["combine_vcfs_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["combine_vcfs_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["combine_vcfs_time"]]
    threads: config["generic"]["combine_vcfs_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"]
    log: config["outputs"]["output_dir"] + "log/merge_vcfs.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools concat -Oz {input.vcfs} > {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


# --force-samples required since samples may have been removed in the het_filter
rule split_by_dataset:
    input:
        vcf = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_imputed_hg38.vcf.gz.csi",
        samples = lambda wildcards: config["inputs"]["dataset_samples"][wildcards.dataset]
    output:
        vcf = config["outputs"]["output_dir"] + "minimac_imputed_by_datatset/{dataset}_{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "minimac_imputed_by_datatset/{dataset}_{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_dataset_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_dataset_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["split_by_dataset_time"]]
    threads: config["imputation"]["split_by_dataset_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"]
    log: config["outputs"]["output_dir"] + "log/split_by_dataset.{dataset}.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.samples} {input.vcf} --force-samples -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

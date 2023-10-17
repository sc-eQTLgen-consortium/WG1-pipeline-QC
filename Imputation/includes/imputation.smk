#!/usr/bin/env python

##############################################
############ hg18, b36, hg19, b37 ############
##############################################

# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.pgen",
        pvar = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.pvar",
        psam = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.psam",
    output:
        bed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_input.bed",
        outbed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_output.bed",
        excluded_ids = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_excluded_ids.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["crossmap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["crossmap_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["crossmap_time"]]
    threads: config["imputation"]["crossmap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink",
        chain_file = config["refs"]["ref_dir"] + (config["refs_extra"]["relative_grch37_to_grch38_chain_path"] if config["inputs"]["genome_build"] in ["hg19", "GRCh37"] else config["refs_extra"]["relative_hg18_to_hg38_chain_path"])
    log: config["outputs"]["output_dir"] + "log/crossmap.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        singularity exec --bind {params.bind} {params.sif} CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        singularity exec --bind {params.bind} {params.sif} awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --exclude {output.excluded_ids} \
            --make-bed \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """


rule sort_bed:
    input:
        bed = config["outputs"]["output_dir"] + ("crossmapped/{ancestry}_crossmapped_plink.bed" if config["inputs"]["genome_build"] in ["hg18", "b36", "hg19", "b37"] else "subset_ancestry/{ancestry}_subset.pgen"),
        bim = config["outputs"]["output_dir"] + ("crossmapped/{ancestry}_crossmapped_plink.bim" if config["inputs"]["genome_build"] in ["hg18", "b36", "hg19", "b37"] else "subset_ancestry/{ancestry}_subset.pvar"),
        fam = config["outputs"]["output_dir"] + ("crossmapped/{ancestry}_crossmapped_plink.fam" if config["inputs"]["genome_build"] in ["hg18", "b36", "hg19", "b37"] else "subset_ancestry/{ancestry}_subset.psam")
    output:
        bed = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bed",
        bim = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bim",
        fam = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.fam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["sort_bed_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["sort_bed_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["sort_bed_time"]]
    threads: config["imputation"]["sort_bed_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        input = ("--bfile " + config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink" if config["inputs"]["genome_build"] in ["hg18", "b36", "hg19", "b37"] else "--pfile " + config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset"),
        out = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted"
    log: config["outputs"]["output_dir"] + "log/sort_bed_after_crossmap.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            {input} \
            --make-bed \
            --max-alleles 2 \
            --out {params.out}
        """

####################################################


rule harmonize_hg38:
    input:
        bed = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bed",
        bim = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bim",
        fam = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["harmonize_hg38_time"]]
    threads: config["imputation"]["harmonize_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted",
        out = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}",
    log: config["outputs"]["output_dir"] + "log/harmonize_hg38.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -jar {params.jar}\
            --input {params.infile} \
            --inputType PLINK_BED \
            --ref {input.ref_vcf} \
            --refType VCF \
            --update-id \
            --output {params.out}
        """


rule plink_to_vcf:
    input:
        bed = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.fam"
    output:
        vcf = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["plink_to_vcf_time"]]
    threads: config["imputation"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "log/plink_to_vcf.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --recode vcf id-paste=iid \
            --chr {params.chr_list} \
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

####################################################################
############ WGS VCF ARE TOO BIG SO NEED TO BE SPLITTED ############
####################################################################

rule split_by_chr_for_harmonize:
    input:
        bed = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bed",
        bim = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.bim",
        fam = config["outputs"]["output_dir"] + "sorted/{ancestry}_sorted.fam"
    output:
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.fam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["split_by_chr_time"]]
    threads: config["imputation"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted"
    log: config["outputs"]["output_dir"] + "log/split_by_chr_for_harmonize.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --chr {wildcards.chr} \
            --make-bed \
            --out {params.out}
        """


rule harmonize_hg38_per_chr:
    input:
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["harmonize_hg38_time"]]
    threads: config["imputation"]["harmonize_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.27-SNAPSHOT/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_sorted",
        out = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}",
    log: config["outputs"]["output_dir"] + "log/harmonize_hg38_per_chr.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem}g -jar {params.jar} \
            --input {params.infile} \
            --inputType PLINK_BED \
            --ref {input.ref_vcf} \
            --refType VCF \
            --update-id \
            --output {params.out}
        """


rule plink_per_chr_to_vcf:
    input:
        bed = expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bed", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        bim = expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bim", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        fam = expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.fam", ancestry=ANCESTRIES, chr=CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_harmonised_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["plink_to_vcf_time"]]
    threads: config["imputation"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}", ancestry=ANCESTRIES, chr=CHROMOSOMES),
        mergelist = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_mergelist.txt",
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "log/plink_to_vcf.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.infiles} | sed 's/ /\\n/g' > {params.mergelist}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --pmerge-list {params.mergelist} bfile \
            --recode vcf id-paste=iid \
            --chr {params.chr_list} \
            --out {params.out}

        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

####################################################################

rule vcf_fixref_hg38:
    input:
        ref_fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi",
        vcf = config["outputs"]["output_dir"] + ("harmonize_hg38_per_chr" if config["settings"]["is_wgs"] else "harmonize_hg38") + "/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + ("harmonize_hg38_per_chr" if config["settings"]["is_wgs"] else "harmonize_hg38") + "/{ancestry}_harmonised_hg38.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["vcf_fixref_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["vcf_fixref_hg38_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["vcf_fixref_hg38_time"]]
    threads: config["imputation"]["vcf_fixref_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/vcf_fixref_hg38.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fixref {input.vcf} -- -f {input.ref_fasta} -i {input.ref_vcf} | \
            singularity exec --bind {params.bind} {params.sif} bcftools norm \
                --check-ref x \
                --fasta-ref {input.ref_fasta} \
                -Oz \
                --output {output.vcf}

        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


# Add tags
# Filter rare and non-HWE variants and those with abnormal alleles and duplicates
rule filter_preimpute_vcf:
    input:
        vcf = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
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


rule het:
    input:
        vcf = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(config["outputs"]["output_dir"] + "het/{ancestry}_filtered_temp.vcf"),
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
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R"
    log: config["outputs"]["output_dir"] + "log/het.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --vcf {output.tmp_vcf} \
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


rule calculate_missingness:
    input:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals = config["outputs"]["output_dir"] + "genotype_donor_annotation/{ancestry}_individuals.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["calculate_missingness_time"]]
    threads: config["imputation"]["calculate_missingness_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "filter_preimpute_vcf/{ancestry}_genotypes"
    log: config["outputs"]["output_dir"] + "log/calculate_missingness.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}
        singularity exec --bind {params.bind} {params.sif} vcftools --gzvcf {output.tmp_vcf} --missing-indv --out {params.out}
        singularity exec --bind {params.bind} {params.sif} bcftools query -l {input.vcf} >> {output.individuals}
        """


rule genotype_donor_annotation:
    input:
        individuals = expand(config["outputs"]["output_dir"] + "genotype_donor_annotation/{ancestry}_individuals.tsv", ancestry=ANCESTRIES),
    output:
        out_temp = temp(config["outputs"]["output_dir"] + "genotype_donor_annotation/genotype_donor_annotation_temp.tsv"),
        combined_individuals = config["outputs"]["output_dir"] + "genotype_donor_annotation/combined_individuals.tsv",
        final = config["outputs"]["output_dir"] + "genotype_donor_annotation/genotype_donor_annotation.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["genotype_donor_annotation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["genotype_donor_annotation_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["genotype_donor_annotation_time"]]
    threads: config["imputation"]["genotype_donor_annotation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        psam = config["outputs"]["output_dir"] + "pca_sex_checks/updated.psam"
    log: config["outputs"]["output_dir"] + "log/genotype_donor_annotation.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cut -f2,5- {params.psam} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.out_temp}
        singularity exec --bind {params.bind} {params.sif} cat {input.individuals} >> {output.combined_individuals}
        singularity exec --bind {params.bind} {params.sif} sed -i '1 i\IID' {output.combined_individuals}
        singularity exec --bind {params.bind} {params.sif} awk -F"\t" 'NR==FNR {{a[$1]; next}} $1 in a' {output.combined_individuals} {output.out_temp} | awk 'BEGIN{{FS=OFS="\t"}}{{sub("1","M",$2);print}}' | awk 'BEGIN{{FS=OFS="\t"}}{{sub("2","F",$2);print}}' > {output.final}
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^IID\tSEX\tProvided_Ancestry/donor_id\tsex\tethnicity_super_population/g' {output.final}
        singularity exec --bind {params.bind} {params.sif} sed -i 's/$/\tsceQTL-Gen_hg38_imputation_pipeline\t1000g_30x_GRCh38_ref/' {output.final}
        singularity exec --bind {params.bind} {params.sif} sed -i '1s/SangerImputationServer/imputation_method/g' {output.final}
        singularity exec --bind {params.bind} {params.sif} sed -i '1s/1000g_30x_GRCh38_ref/imputation_reference/g' {output.final}
        """


rule split_by_chr:
    input:
        vcf = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        index = config["outputs"]["output_dir"] + "split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["split_by_chr_time"]]
    threads:
        config["imputation"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/split_by_chr.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = config["outputs"]["output_dir"] + "split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = config["refs"]["ref_dir"] + config["refs_extra"]["relative_map_path"],
        phasing_bcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "chr{chr}.bcf",
        phasing_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "chr{chr}.bcf.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr_{chr}_phased.vcf.gz",
        index = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr_{chr}_phased.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["eagle_prephasing_time"]]
    threads: config["imputation"]["eagle_prephasing_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr_{chr}_phased"
    log: config["outputs"]["output_dir"] + "log/eagle_prephasing.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} eagle \
            --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_bcf} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


# Results are slightly different compared to v1.0.2
# ID column is rsID in NEW but chr:pos:ref:alt in OLD
# NEW INFO column contains AVG_CS + values are slightly different
# NEW FORMAT is reordered as GT:GP:DS
# temp-prefix parameter fixes issue I had running it on our cluster
rule minimac_imputation:
    input:
        reference = config["refs"]["ref_dir"] + config["refs_extra"]["relative_imputation_dir"] + "chr{chr}.msav",
        target = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr_{chr}_phased.vcf.gz",
        target_index = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr_{chr}_phased.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr_{chr}.dose.vcf.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["minimac_imputation_time"]]
    threads: config["imputation"]["minimac_imputation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "minimac_imputed",
        temp_prefix = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr_{chr}_m4_",
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


rule combine_vcfs_ancestry:
    input:
        vcfs = expand(config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr_{chr}.dose.vcf.gz", ancestry=ANCESTRIES, chr=CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["combine_vcfs_ancestry_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["combine_vcfs_ancestry_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["combine_vcfs_ancestry_time"]]
    threads: config["imputation"]["combine_vcfs_ancestry_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        files_begin = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    log: config["outputs"]["output_dir"] + "log/combine_vcfs_ancestry.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools concat -Oz {params.files_begin} > {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule split_per_dataset:
    input:
        vcf = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi",
        samples = lambda wildcards: config["inputs"]["dataset_samples"][wildcards.dataset]
    output:
        vcf = config["outputs"]["output_dir"] + "split_per_dataset/{dataset}_{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "split_per_dataset/{dataset}_{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_per_dataset_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_per_dataset_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["imputation"]["split_per_dataset_time"]]
    threads: config["imputation"]["split_per_dataset_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"]
    log: config["outputs"]["output_dir"] + "log/split_per_dataset.{dataset}.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.samples} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

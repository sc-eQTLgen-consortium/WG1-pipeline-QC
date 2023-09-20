#!/usr/bin/env python


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
        inbed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmap_input.bed",
        outbed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmap_output.bed",
        excluded_ids = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_excluded_ids.txt"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["crossmap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["crossmap_memory"]
    threads: config["imputation"]["crossmap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink",
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    log: config["outputs"]["output_dir"] + "logs/crossmap.{ancestry}.log"
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
            --output-chr MT \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """


rule sort_bed:
    input:
        bed = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = config["outputs"]["output_dir"] + "crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["sort_bed_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["sort_bed_memory"]
    threads: config["imputation"]["sort_bed_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted"
    log: config["outputs"]["output_dir"] + "logs/sort_bed.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --max-alleles 2 \
            --output-chr MT \
            --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"]
    threads: config["imputation"]["harmonize_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}",
    log: config["outputs"]["output_dir"] + "logs/harmonize_hg38.{ancestry}.log"
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"]
    threads: config["imputation"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "logs/plink_to_vcf.{ancestry}.log"
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
        bed = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = config["outputs"]["output_dir"] + "crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    output:
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.fam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"]
    threads: config["imputation"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted"
    log: config["outputs"]["output_dir"] + "logs/split_by_chr_for_harmonize.{ancestry}.chr_{chr}.log"
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
        bed = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.bed",
        bim = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.bim",
        fam = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted.fam",
        ref_vcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"],
        ref_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_vcf_path"] + ".tbi"
    output:
        bed = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bed",
        bim = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bim",
        fam = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_memory"]
    threads: config["imputation"]["harmonize_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar",
        infile = config["outputs"]["output_dir"] + "split_by_chr_for_harmonize/{ancestry}_chr_{chr}_crossmapped_sorted",
        out = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}",
    log: config["outputs"]["output_dir"] + "logs/harmonize_hg38_per_chr.{ancestry}.chr_{chr}.log"
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
        bed = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bed", chr = CHROMOSOMES),
        bim = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.bim", chr = CHROMOSOMES),
        fam = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}.fam", chr = CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["plink_to_vcf_memory"]
    threads: config["imputation"]["plink_to_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_chr_{chr}", chr = CHROMOSOMES),
        mergelist = config["outputs"]["output_dir"] + "harmonize_hg38_per_chr/{ancestry}_mergelist.txt",
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38"
    log: config["outputs"]["output_dir"] + "logs/plink_to_vcf.{ancestry}.log"
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
        vcf = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["vcf_fixref_hg38_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["vcf_fixref_hg38_memory"]
    threads: config["imputation"]["vcf_fixref_hg38_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "logs/vcf_fixref_hg38.{ancestry}.log"
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["filter_preimpute_vcf_memory"]
    threads: config["imputation"]["filter_preimpute_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        hwe = config["imputation"]["filter_preimpute_vcf_snp_hwe"],
        missing = config["imputation"]["filter_preimpute_vcf_snp_missing_pct"],
        maf = lambda wildcards: ANCESTRY_MAF_DICT[wildcards.ancestry]
    log: config["outputs"]["output_dir"] + "logs/filter_preimpute_vcf.{ancestry}.log"
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_memory"]
    threads: config["imputation"]["het_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "het/{ancestry}_het",
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R"
    log: config["outputs"]["output_dir"] + "logs/het.{ancestry}.log"
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"]
    threads: config["imputation"]["het_filter_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "logs/het_filter.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule genotype_donor_annotation:
    input:
        vcf = expand(config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz", ancestry = ANCESTRIES),
        index = expand(config["outputs"]["output_dir"] + "het_filter/{ancestry}_het_filtered.vcf.gz.csi", ancestry = ANCESTRIES),
        psam = config["outputs"]["output_dir"] + "pca_sex_checks/sex_ancestry_updated.psam"
    output:
        individuals = config["outputs"]["output_dir"] + "genotype_donor_annotation/{ancestry}_individuals.tsv",
        out_temp = temp(config["outputs"]["output_dir"] + "genotype_donor_annotation/genotype_donor_annotation_temp.tsv"),
        combined_individuals = config["outputs"]["output_dir"] + "genotype_donor_annotation/combined_individuals.tsv",
        final = config["outputs"]["output_dir"] + "genotype_donor_annotation/genotype_donor_annotation.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["genotype_donor_annotation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["genotype_donor_annotation_memory"]
    threads: config["imputation"]["genotype_donor_annotation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
    log: config["outputs"]["output_dir"] + "logs/genotype_donor_annotation.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools query -l {input.vcf} >> {output.individuals}
        singularity exec --bind {params.bind} {params.sif} cut -f2,5- {input.psam} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1,1"}}' > {output.out_temp}
        singularity exec --bind {params.bind} {params.sif} cat {output.individuals} >> {output.combined_individuals}
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
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["split_by_chr_memory"]
    threads:
        config["imputation"]["split_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "logs/split_by_chr.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = config["outputs"]["output_dir"] + "split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = config["refs"]["ref_dir"] + config["refs_extra"]["relative_map_path"],
        phasing_bcf = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "ls",
        phasing_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_phasing_dir"] + "chr{chr}.bcf.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["eagle_prephasing_memory"]
    threads: config["imputation"]["eagle_prephasing_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr{chr}_phased"
    log: config["outputs"]["output_dir"] + "logs/eagle_prephasing.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} eagle \
            --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_bcf} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        """


rule minimac_imputation:
    input:
        impute_file = config["refs"]["ref_dir"] + config["refs_extra"]["relative_imputation_dir"] + "chr{chr}.m3vcf.gz",
        vcf = config["outputs"]["output_dir"] + "eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    output:
        dose = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["minimac_imputation_memory"]
    threads: config["imputation"]["minimac_imputation_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        minimac4 = "/opt/bin/minimac4",
        out = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr{chr}",
        chunk_length = config["imputation"]["minimac_chunk_length"]
    log: config["outputs"]["output_dir"] + "logs/minimac_imputation.{ancestry}.chr_{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.minimac4} \
            --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix {params.out} \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads} \
            --ChunkLengthMb {params.chunk_length}
        """


rule combine_vcfs_ancestry:
    input:
        vcfs = expand(config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", chr = CHROMOSOMES, ancestry = ANCESTRIES)
    output:
        vcf = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["combine_vcfs_ancestry_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["imputation"]["combine_vcfs_ancestry_memory"]
    threads: config["imputation"]["combine_vcfs_ancestry_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        files_begin = config["outputs"]["output_dir"] + "minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    log: config["outputs"]["output_dir"] + "logs/combine_vcfs_ancestry.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools concat -Oz {params.files_begin} > {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """

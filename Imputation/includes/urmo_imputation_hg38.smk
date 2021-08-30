
#!/usr/bin/env python
shell.executable('bash')



# Converts BIM to BED and converts the BED file via CrossMap. 
# Finds excluded SNPs and removes them from the original plink file. 
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pgen",
        psam = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.psam",
        pvar = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pvar",
        chain_file = "/opt/chain/GRCh37_to_GRCh38.chain"
    output:
        bed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmap_input.bed",
        outbed = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmap_output.bed",
        ids = output_dict["output_dir"] + "/crossmapped/{ancestry}_input_ids.txt",
        out_ids = output_dict["output_dir"] + "/crossmapped/{ancestry}_output_ids.txt",
        excluded_ids = output_dict["output_dir"] + "/crossmapped/{ancestry}_excluded_ids.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["crossmap_memory"]
    threads:
        imputation_dict["crossmap_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/crossmapped/{ancestry}_subset",
        out = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink"
    shell: 
        """
        singularity exec --bind {params.bind} {params.sif} awk '{{print $1,$4,$4+1,$2,$5,$6,$2 "___" $5 "___" \$6}}' {input.pvar} > {output.inbed}
        singularity exec --bind {params.bind} {params.sif} CrossMap.py bed {input.chain_file} {output.inbed} {output.outbed}
        singularity exec --bind {params.bind} {params.sif} awk '{{print $7}}' {output.inbed} | \
            singularity exec --bind {params.bind} {params.sif} sort > {output.ids}
        singularity exec --bind {params.bind} {params.sif} awk '{{print $7}}' {output.outbed} | \
            singularity exec --bind {params.bind} {params.sif} sort > {output.out_ids}
        singularity exec --bind {params.bind} {params.sif} comm -23 {output.ids} {output.out_ids} | \
            singularity exec --bind {params.bind} {params.sif} awk '{{split($0,a,"___"); print a[1]}'} > {output.excluded_ids}
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile ${study_name_bed.simpleName} --exclude {output.excluded_ids} --make-bed --output-chr MT --out {params.out}
        singularity exec --bind {params.bind} {params.sif} awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$5,$6}}' {output.outbed} > {output.bim}
        """


rule sort_bed:
    input:
        pgen = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.pgen",
        psam = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.psam",
        pvar = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink.pvar"
    output:
        bed = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["sort_bed_memory"]
    threads:
        imputation_dict["sort_bed_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/crossmapped/{ancestry}_crossmapped_plink",
        out = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.infile} --make-bed --output-chr MT --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi"
        # tuple file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes_hg38.collect()
    output:
        bed = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bed",
        bim = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bim",
        fam = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["harmonize_hg38_memory"]
    threads:
        imputation_dict["harmonize_hg38_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx25g -jar /usr/bin/GenotypeHarmonizer.jar\
            --input {params.infile}\
            --inputType PLINK_BED\
            --ref {input.vcf}\
            --refType VCF\
            --update-id\
            --output {params.out}
        """


rule plink_to_vcf:
    input:
        bed = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bed",
        bim = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.bim",
        fam = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}.fam"
    output:
        output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["plink_to_vcf_memory"]
    threads:
        imputation_dict["plink_to_vcf_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}",
        out = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.infile} --recode vcf-iid --chr 1-22 --out {params.out}
        """


rule vcf_fixref_hg38:
    input:
        data_vcf = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf",
        fasta = fasta,
        vcf = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = vcf_dir + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi"
    output:
        data_vcf_gz = output_dict["output_dir"] + "/harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        vcf = output_dict["output_dir"] + "/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["vcf_fixref_hg38_memory"]
    threads:
        imputation_dict["vcf_fixref_hg38_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    script:
        """
        singularity exec --bind {params.bind} {params.sif} bgzip {input.data_vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.data_vcf_gz}
        
        singularity exec --bind {params.bind} {params.sif} bcftools +fixref {output.data_vcf_gz} -- -f {input.fasta} -i ${input.vcf} | \
        singularity exec --bind {params.bind} {params.sif} bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = output_dict["output_dir"] + "/vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["filter_preimpute_vcf_memory"]
    threads:
        imputation_dict["filter_preimpute_vcf_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],

    script:
        """
        #Index
        singularity exec --bind {params.bind} {params.sif} bcftools index {input.vcf}

        #Add tags
        singularity exec --bind {params.bind} {params.sif} bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        singularity exec --bind {params.bind} {params.sif} bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' {output.tagged_vcf} |\
        singularity exec --bind {params.bind} {params.sif} bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        singularity exec --bind {params.bind} {params.sif} bcftools filter -e "ALT='.'" |\
        singularity exec --bind {params.bind} {params.sif} bcftools norm -d all |\
        singularity exec --bind {params.bind} {params.sif} bcftools norm -m+any |\
        singularity exec --bind {params.bind} {params.sif} bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        #Index the output file
        bcftools index filtered.vcf.gz
        """

rule het:
    input:
        vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        inds = output_dict["output_dir"] + "/het/{ancestry}_het_failed.inds",
        het = output_dict["output_dir"] + "/het/{ancestry}_het.het",
        passed = output_dict["output_dir"] + "/het/{ancestry}_het_passed.inds",
        passed_list = output_dict["output_dir"] + "/het/{ancestry}_het_passed.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_memory"]
    threads: imputation_dict["het_threads"]
    params:
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R",
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het/{ancestry}_het",
        sif = input_dict["singularity_image"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} vcftools --gzvcf {input.vcf} --het --out {output.het}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list = output_dict["output_dir"] + "/het/{ancestry}_het_passed.txt",
        vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz",
        index = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["het_filter_memory"]
    threads: imputation_dict["het_filter_threads"]
    params:
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}

        #Index the output file
        bcftools index {output.vcf}
        """


rule calculate_missingness:
    input:
        filtered_vcf = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_genotypes.imiss"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["calculate_missingness_memory"]
    threads:
        imputation_dict["calculate_missingness_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_genotypes"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} vcftools --gzvcf {input.tagged_vcf} --missing-indv --out {params.out}
        """


rule split_by_chr:
    input:
        filtered_vcf = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = output_dict["output_dir"] + "/filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    output:
        output_dict["output_dir"] + "/split_by_chr/{ancestry}_chr_{chr}.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["split_by_chr_memory"]
    threads:
        imputation_dict["split_by_chr_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"]
    script:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output}
        """

 
rule eagle_prephasing:
    input:
        vcf = output_dict["output_dir"] + "/split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = genetic_map,
        phasing_file = phasing_dir + "/chr{chr}.bcf"
    output:
        output_dict["output_dir"] + "{ancestry}_chr{chr}_phased.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["eagle_prephasing_memory"]
    threads: imputation_dict["eagle_prephasing_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "{ancestry}_chr{chr}_phased"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools index {input.vcf}
        singularity exec --bind {params.bind} {params.sif} eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        """


rule minimac_imputation:
    input:
        vcf = output_dict["output_dir"] + "{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file = impute_dir + "chr{chr}.m3vcf.gz"
    output:
        output_dict["output_dir"] + "{ancestry}_chr{chr}.dose.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["minimac_imputation_memory"]
    threads:
        imputation_dict["minimac_imputation_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "{ancestry}_chr{chr}"
    script:
        """
        singularity exec --bind {params.bind} {params.sif} minimac4 --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix {params.out} \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads}
        """




rule genotype_donor_annotation:
    input:
        het_psams = lambda wildcards: expand(output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.psam", ancestry = ancestry_file["unique_ancestry"]),
        updated_psam = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam"
    output:
        output_dict["output_dir"] + "/genotype_donor_annotation.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["genotype_donor_annotation_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * imputation_dict["genotype_donor_annotation_memory"]
    threads: imputation_dict["genotype_donor_annotation_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out_temp = output_dict["output_dir"] + "/genotype_donor_annotation_temp.tsv",
        out_temp2 = output_dict["output_dir"] + "/genotype_donor_annotation_temp2.tsv"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} head -n 1 {input.updated_psam}  > {params.out_temp}
        singularity exec --bind {params.bind} {params.sif} cat {input.het_psams} | grep -v "#" >> {params.out_temp}
        singularity exec --bind {params.bind} {params.sif} sed -r 's/\t/_/' {params.out_temp} | \
        singularity exec --bind {params.bind} {params.sif} sed 's/\t1/\tM/g' | \
        singularity exec --bind {params.bind} {params.sif} sed 's/\t2/\tF/g' | \
        singularity exec --bind {params.bind} {params.sif} sed 's/#FID_IID/donor_id/g' > {params.out_temp2}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN {{FS="\t"; OFS="\t"}}{{$2=$3="";gsub("\t+","\t",$0)}}1' {params.out_temp2} > {output}
        singularity exec --bind {params.bind} {params.sif} sed -i 's/$/\tsceQTL-Gen_hg38_imputation_pipeline/' {output}
        singularity exec --bind {params.bind} {params.sif} sed -i '1s/SangerImputationServer/imputation_server/g' {output}
        singularity exec --bind {params.bind} {params.sif} rm {params.out_temp}
        singularity exec --bind {params.bind} {params.sif} rm {params.out_temp2}
        """
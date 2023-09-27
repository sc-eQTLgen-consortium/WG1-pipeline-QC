#!/usr/bin/env python


######################################################
############ IF THE INPUT IS A SINGLE VCF ############
######################################################

rule split_input_vcf_by_chr:
    input:
        vcf = config["inputs"]["genotype_path"] + ".vcf.gz",
        index = config["inputs"]["genotype_path"] + ".vcf.gz.tbi"
    output:
        vcf = config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz",
        index = config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["split_vcf_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["split_vcf_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["split_vcf_by_chr_time"]]
    threads: config["pre_processing"]["split_vcf_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/split_input_vcf_by_chr.{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools view -r {wildcards.chr} {input.vcf} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule normalise_vcf_by_chr:
    input:
        vcf = config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz",
        index = config["outputs"]["output_dir"] + "split_input_vcf_by_chr/chr_{chr}.vcf.gz.csi"
    output:
        vcf = config["outputs"]["output_dir"] + "normalise/chr_{chr}_normalised.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["normalise_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["normalise_by_chr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["normalise_by_chr_time"]]
    threads: config["pre_processing"]["normalise_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/normalise_vcf_by_chr.{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools norm -m any {input.vcf} -o {output.vcf}
        """


#####################################################
############ IF THE INPUT IS VCF PER CHR ############
#####################################################


rule normalise_input_vcf_per_chr:
    input:
        vcf = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz",
        index = lambda wildcards: config["inputs"]["genotype_path"].replace("CHR", wildcards.chr) + ".vcf.gz.tbi",
    output:
        vcf = config["outputs"]["output_dir"] + "normalise/chr_{chr}_normalised.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["normalise_by_chr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["normalise_by_chr_memory"],
        time = lambda wildcards,attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["normalise_by_chr_time"]]
    threads: config["pre_processing"]["normalise_by_chr_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/normalise_input_vcf_per_chr.{chr}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools norm -m any {input.vcf} -o {output.vcf}
        """


#########################################################
############ FILTER AND CONCAT THE WGS INPUT ############
#########################################################


rule wgs_filter:
    input:
        vcf = config["outputs"]["output_dir"] + "normalise/chr_{chr}_normalised.vcf.gz"
    output:
        vcf = config["outputs"]["output_dir"] + "wgs_filter/chr_{chr}_normalised-filtered.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["wgs_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["wgs_filter_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["wgs_filter_time"]]
    threads: config["pre_processing"]["wgs_filter_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/custom_vcf_filter.py",
        out_dir = config["outputs"]["output_dir"] + "wgs_filter",
        sex = lambda wildcards: "--sex {}".format(None) if wildcards.chr in ["X", "Y"] else "", # TODO sex file does not exist,
        genotype_quality = config["wgs_filter"]["genotype_quality"],
        allelic_balance_lower = config["wgs_filter"]["allelic_balance_lower"],
        allelic_balance_upper = config["wgs_filter"]["allelic_balance_upper"],
        inbreeding_coefficient = config["wgs_filter"]["inbreeding_coefficient"],
        no_snv_vqsr_check = config["wgs_filter"]["no_snv_vqsr_check"],
        vqsr_snv = config["wgs_filter"]["vqsr_snv"],
        no_indel_vqsr_check = config["wgs_filter"]["no_indel_vqsr_check"],
        vqsr_indel = config["wgs_filter"]["vqsr_indel"],
        keep_multialleic = config["wgs_filter"]["keep_multialleic"],
        keep_non_pass_snv = config["wgs_filter"]["keep_non_pass_snv"],
        keep_non_pass_indel = config["wgs_filter"]["keep_non_pass_indel"],
        keep_low_complexity = config["wgs_filter"]["keep_low_complexity"],
        minor_allele_frequency = config["wgs_filter"]["minor_allele_frequency"],
        call_rate = config["wgs_filter"]["call_rate"],
        hardy_weinberg_equilibrium = config["wgs_filter"]["hardy_weinberg_equilibrium"],
        filtered_depth = config["wgs_filter"]["filtered_depth"],
        keep_info_column = config["wgs_filter"]["keep_info_column"]
    log: config["outputs"]["output_dir"] + "log/wgs_filter.{chr}.log"
    shell:
        """
        mkdir -p {params.out_dir} && \
            singularity exec --bind {params.bind} {params.sif} python {params.script} \
                --input {input.vcf} \
                --output {output.vcf} \
                --log {log} \
                {params.sex} \
                --genotype_quality {params.genotype_quality} \
                --allelic_balance_lower {params.allelic_balance_lower} \
                --allelic_balance_upper {params.allelic_balance_upper} \
                --inbreeding_coefficient {params.inbreeding_coefficient} \
                --no_snv_vqsr_check {params.no_snv_vqsr_check} \
                --vqsr_snv {params.vqsr_snv} \
                --no_indel_vqsr_check {params.no_indel_vqsr_check} \
                --vqsr_indel {params.vqsr_indel} \
                --keep_multialleic {params.keep_multialleic} \
                --keep_non_pass_snv {params.keep_non_pass_snv} \
                --keep_non_pass_indel {params.keep_non_pass_indel} \
                --keep_low_complexity {params.keep_low_complexity} \
                --minor_allele_frequency {params.minor_allele_frequency} \
                --call_rate {params.call_rate} \
                --hardy_weinberg_equilibrium {params.hardy_weinberg_equilibrium} \
                --filtered_depth {params.filtered_depth} \
                --keep_info_column {params.keep_info_column}
        """


rule combine_wgs_filtered_vcfs:
    input:
        vcf = lambda wildcards: expand(config["outputs"]["output_dir"] + "wgs_filter/chr_{chr}_normalised-filtered.vcf.gz", chr=CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "combine_wgs_filtered_vcfs/normalised-filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "combine_wgs_filtered_vcfs/normalised-filtered.vcf.gz.csi",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_vcfs_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_vcfs_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["combine_vcfs_time"]]
    threads: config["pre_processing"]["combine_vcfs_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        file_list = config["outputs"]["output_dir"] + "combine_wgs_filtered_vcfs/file_list.txt"
    log: config["outputs"]["output_dir"] + "log/combine_wgs_filtered_vcfs.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {input.vcf} | sed 's/ /\\n/g' > {params.file_list}
        singularity exec --bind {params.bind} {params.sif} bcftools concat --file-list {params.file_list} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
    """


rule wgs_filtered_vcf_to_pgen:
    input:
        vcf = config["outputs"]["output_dir"] + "combine_wgs_filtered_vcfs/normalised-filtered.vcf.gz",
        index = config["outputs"]["output_dir"] + "combine_wgs_filtered_vcfs/normalised-filtered.vcf.gz.csi",
    output:
        pgen = config["outputs"]["output_dir"] + "wgs_filtered_vcf_to_pgen/data.pgen",
        pvar = config["outputs"]["output_dir"] + "wgs_filtered_vcf_to_pgen/data.pvar",
        psam = config["outputs"]["output_dir"] + "wgs_filtered_vcf_to_pgen/data.psam",
        tmp_psam = config["outputs"]["output_dir"] + "wgs_filtered_vcf_to_pgen/tmp.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["vcf_to_pgen_time"]]
    threads: config["pre_processing"]["vcf_to_pgen_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data"
    log: config["outputs"]["output_dir"] + "log/wgs_filtered_vcf_to_pgen.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --vcf {input.vcf} \
            --make-pgen \
            --threads {threads} \
            --out {params.out}
        
        singularity exec --bind {params.bind} {params.sif} cp {output.psam} {output.tmp_psam}
        singularity exec --bind {params.bind} {params.sif} echo "#FID\tIID\tPAT\tMAT\tSEX\tProvided_Ancestry\tgenotyping_platform\tarray_available\twgs_available\twes_available\tage\tage_range\tStudy\tsmoking_status\thormonal_contraception_use_currently\tmenopause\tpregnancy_status" > {output.psam}
        singularity exec --bind {params.bind} {params.sif} cut -f1 {output.tmp_psam} | sed -e 's/^/\t/' | tail -n +2 >> {output.psam}
        """

###############################################################
############ PREPARE FOR QC: OPTION 1 = SINGLE VCF ############
###############################################################

rule input_vcf_to_pgen:
    input:
        vcf = config["inputs"]["genotype_path"] + ".vcf.gz",
        index = config["inputs"]["genotype_path"] + ".vcf.gz.tbi"
    output:
        pgen = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pgen",
        pvar = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pvar",
        psam = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.psam",
        tmp_psam = config["outputs"]["output_dir"] + "wgs_filtered_vcf_to_pgen/tmp.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["vcf_to_pgen_time"]]
    threads: config["pre_processing"]["vcf_to_pgen_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data"
    log: config["outputs"]["output_dir"] + "log/input_vcf_to_pgen.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --vcf {input.vcf} \
            --make-pgen \
            --out {params.out}
            
        singularity exec --bind {params.bind} {params.sif} cp {output.psam} {output.tmp_psam}
        singularity exec --bind {params.bind} {params.sif} echo -e "#FID\tIID\tPAT\tMAT\tSEX\tProvided_Ancestry\tgenotyping_platform\tarray_available\twgs_available\twes_available\tage\tage_range\tStudy\tsmoking_status\thormonal_contraception_use_currently\tmenopause\tpregnancy_status" > {output.psam}
        singularity exec --bind {params.bind} {params.sif} cut -f1 {output.tmp_psam} | sed -e 's/^/\t/' | tail -n +2 >> {output.psam}
        """


#######################################################################
############ PREPARE FOR QC: OPTION 3 = VCF PER CHROMOSOME ############
#######################################################################


rule combine_input_vcfs:
    input:
        vcf = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}") + ".vcf.gz", chr = CHROMOSOMES),
        indices = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}") + ".vcf.gz.tbi", chr = CHROMOSOMES)
    output:
        vcf = config["outputs"]["output_dir"] + "combine_input_vcfs/data.vcf.gz",
        index = config["outputs"]["output_dir"] + "combine_input_vcfs/data.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_vcfs_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_vcfs_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["combine_vcfs_time"]]
    threads: config["pre_processing"]["combine_vcfs_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        file_list = config["outputs"]["output_dir"] + "combine_input_vcfs/file_list.txt"
    log: config["outputs"]["output_dir"] + "log/combine_input_vcfs.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {input.vcf} | sed 's/ /\\n/g' > {params.file_list}
        singularity exec --bind {params.bind} {params.sif} bcftools concat --file-list {params.file_list} -Oz -o {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


rule combined_input_vcf_to_pgen:
    input:
        vcf = config["outputs"]["output_dir"] + "combine_input_vcfs/data.vcf.gz",
        index = config["outputs"]["output_dir"] + "combine_input_vcfs/data.vcf.gz.csi"
    output:
        pgen = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pgen",
        pvar = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pvar",
        psam = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["vcf_to_pgen_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["vcf_to_pgen_time"]]
    threads: config["pre_processing"]["vcf_to_pgen_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data"
    log: config["outputs"]["output_dir"] + "log/combined_input_vcf_to_pgen.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --vcf {input.vcf} \
            --make-pgen \
            --out {params.out}
        """


#########################################################################
############ PREPARE FOR QC: OPTION 4 = PLINK PER CHROMOSOME ############
#########################################################################


rule combine_input_plink:
    input:
        pgen = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}") + ".pgen", chr = CHROMOSOMES),
        pvar = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}") + ".pvar", chr = CHROMOSOMES),
        psam = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}") + ".psam", chr = CHROMOSOMES)
    output:
        pgen = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pgen",
        pvar = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.pvar",
        psam = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["pre_processing"]["combine_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["pre_processing"]["combine_plink_time"]]
    threads: config["pre_processing"]["combine_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infiles = lambda wildcards: expand(config["inputs"]["genotype_path"].replace("CHR", "{chr}"), chr = CHROMOSOMES),
        pmerge_list = config["outputs"]["output_dir"] + "harmonize_hg38/pmerge_list.txt",
        out = config["outputs"]["output_dir"] + "plink_gender_ancestry_input/data"
    log: config["outputs"]["output_dir"] + "log/combine_input_plink.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.infiles} | sed 's/ /\\n/g' > {params.pmerge_list}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pmerge-list {params.pmerge_list} bfile \
            --make-pgen \
            --out {params.out}
        """
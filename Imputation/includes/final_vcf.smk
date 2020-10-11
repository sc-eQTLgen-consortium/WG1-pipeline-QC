#!/usr/local/envs/py36/bin python3


if os.path.exists(output_dict["output_dir"] + "/gtm_projection/european.inds"): ### Path to European file
    europeans = pd.read_csv(output_dict["output_dir"] + "/gtm_projection/european.inds", sep = "\t")
    if (len(europeans)) > 0:
        rule bed2vcf_european:
            input:
                european = output_dict["output_dir"] + "/gtm_projection/european.inds",
                bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}.bed"
            output:
                output_dict["output_dir"] + "/vcf/european_QC_filtered_chr{chr}.vcf"
            resources:
                mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2vcf_european_memory"],
                disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2vcf_european_memory"]
            threads: final_vcf_dict["bed2vcf_european_threads"]
            params:
                sif = input_dict["singularity_image"],
                bfile = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}
            shell:
                """
                singularity exec {params.sif} plink --bfile {input.bfile} --recode --keep {input.european} vcf --out {output}
                """

        rule european_vcf_sort:
            input:
                output_dict["output_dir"] + "/vcf/european_QC_filtered_chr{chr}.vcf"
            output:
                output_dict["output_dir"] + "/vcf/european_QC_filtered_sorted_chr{chr}.vcf.gz"
            resources:
                mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["european_vcf_sort_memory"],
                disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["european_vcf_sort_memory"]
            threads: final_vcf_dict["european_vcf_sort_threads"]
            params:
                sif = input_dict["singularity_image"]
            shell:
                """
                singularity exec {params.sif} bcftools sort {input} -Oz -o {output}
                """


if os.path.exists(output_dict["output_dir"] + "/gtm_projection/non_european.inds"): ### Path to European file
    non_europeans = pd.read_csv(output_dict["output_dir"] + "/gtm_projection/non_european.inds", sep = "\t")
    if (len(non_europeans)) > 0:
        rule bed2vcf_non_european:
            input:
                non_european = output_dict["output_dir"] + "/gtm_projection/non_european.inds",
                bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}.bed"
            output:
                output_dict["output_dir"] + "/vcf/non_european_QC_filtered_chr{chr}.vcf"
            resources:
                mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2vcf_non_european_memory"],
                disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2vcf_non_european_memory"]
            threads: final_vcf_dict["bed2vcf_non_european_threads"]
            params:
                sif = input_dict["singularity_image"],
                bfile = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}
            shell:
                """
                singularity exec {params.sif} plink --bfile {input.bfile} --recode --keep {input.non_european} vcf --out {output}
                """

        rule non_european_vcf_sort:
            input:
                output_dict["output_dir"] + "/vcf/non_european_QC_filtered_chr{chr}.vcf"
            output:
                output_dict["output_dir"] + "/vcf/non_european_QC_filtered_sorted_chr{chr}.vcf"
            resources:
                mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["european_vcf_sort_memory"],
                disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["european_vcf_sort_memory"]
            threads: final_vcf_dict["bed2vcf_non_european_memory"]
            params:
                sif = input_dict["european_vcf_sort_threads"]
            shell:
                """
                singularity exec {params.sif} bcftools sort {input} -Oz -o {output}
                """


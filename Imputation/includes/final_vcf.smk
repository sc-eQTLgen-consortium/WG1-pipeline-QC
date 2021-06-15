#!/usr/local/envs/py36/bin python3

rule bed2pgen:
    input:
        output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand-updated-chr23.vcf"
    output:
        output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pgen"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2pgen_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["bed2pgen_memory"]
    threads: final_vcf_dict["bed2pgen_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        bfile = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand-updated-chr{chr}",
        out = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --bfile {params.bfile} --make-pgen 'psam-cols='fid,parents,sex,phenos 'id-delim='_ 'id-paste='fid,iid --out {params.out}
        """

rule fix_pvar:
    input:
        pgen = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pgen",
        pvar = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pvar"
    output:
        old = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.old",
        temp = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.temp",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["fix_pvar_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["fix_pvar_memory"]
    threads: final_vcf_dict["fix_pvar_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        pfile = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cp {params.pfile}.pvar {output.old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.old} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$3,$4,$5)}}' > {output.temp}
        singularity exec --bind {params.bind} {params.sif} grep "#" {input.pvar} | sed 's/, /,/g' > {params.pfile}.pvar
        singularity exec --bind {params.bind} {params.sif} awk -F"\t" 'BEGIN{{OFS=FS = "\t"}} NR==FNR{{a[$1 FS $2 FS $3] = $0; next}} {{ind = $1 FS $2 FS $3}} ind in a {{print a[ind], $6, $7}}' {output.temp} {input.pvar} >> {params.pfile}.pvar
        """


rule pgen2vcf:
    input:
        output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}.pvar.temp"
    output:
        # vcf = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf",
        gvcf = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["pgen2vcf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["pgen2vcf_memory"]
    threads: final_vcf_dict["pgen2vcf_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        pfile = output_dict["output_dir"] + "/convert_files/{ancestry}/{ancestry}_strand_check-updated-chr{chr}",
        out = output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.pfile} --recode vcf --out {params.out}
        singularity exec --bind {params.bind} {params.sif} bgzip {params.out}.vcf
        singularity exec --bind {params.bind} {params.sif} tabix -p vcf {output.gvcf}
        """

# if os.path.exists(output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv"):
    # ancestry_file = pd.read_csv(output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv", sep = "\t", header = None, names = ["Ancestry"])

rule combine_ancestries:
    input:
        vcfs = lambda wildcards: expand(output_dict["output_dir"] + "/vcf/{ancestry}/{ancestry}_QC_filtered_chr{chr}.vcf.gz", chr = wildcards.chr, ancestry = ancestry_file["unique_ancestry"])
    output:
        output_dict["output_dir"] + "/vcf/combined/QC_filtered_chr{chr}.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_ancestries_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_ancestries_memory"]
    threads: final_vcf_dict["combine_ancestries_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        bfile = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} vcf-merge {input.vcfs} > {output}
        """


rule vcf_sort:
    input:
        output_dict["output_dir"] + "/vcf/combined/QC_filtered_chr{chr}.vcf"
    output:
        output_dict["output_dir"] + "/vcf/combined_sorted/QC_filtered_sorted_chr{chr}.vcf.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["vcf_sort_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["vcf_sort_memory"]
    threads: final_vcf_dict["vcf_sort_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools sort {input} -Oz -o {output}
        """


rule combine_vcfs:
    input:
        vcfs = expand(output_dict["output_dir"] + "/vcf/combined_sorted/QC_filtered_sorted_chr{chr}.vcf.gz", chr = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]),
        sexes = output_dict["output_dir"] + "/vcf/files4submission/samples.txt"
    output:
        combined = output_dict["output_dir"] + "/vcf/combined_sorted/QC_filtered_sorted_chr.vcf.gz",
        updated = output_dict["output_dir"] + "/vcf/combined_sorted/QC_filtered_sorted_chr_updated.vcf.gz",
        final = output_dict["output_dir"] + "/vcf/files4submission/QC_filtered_sorted_chr_updated_final.vcf.gz",
        ind = output_dict["output_dir"] + "/vcf/files4submission/QC_filtered_sorted_chr_updated_final.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_vcfs_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["combine_vcfs_memory"]
    threads: final_vcf_dict["combine_vcfs_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        conversion = "/opt/WG1-pipeline-QC/Imputation/chr_conversions.tsv"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools concat -Oz {input.vcfs} > {output.combined}
        singularity exec --bind {params.bind} {params.sif} bcftools annotate -Oz --rename-chrs {params.conversion} {output.combined} > {output.updated}
        singularity exec --bind {params.bind} {params.sif} bcftools +fixploidy {output.updated} -Oz -o {output.final} -- -s {input.sexes}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.final}
        """

rule sex4imputation:
    input:
        output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam"
    output:
        output_dict["output_dir"] + "/vcf/files4submission/samples.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["sex4imputation_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["sex4imputation_memory"]
    threads: final_vcf_dict["sex4imputation_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS="\t"}}{{OFS=""}}{{print($1,"_",$2,"\t",$5)}}' {input} | \
            singularity exec --bind {params.bind} {params.sif} sed '1d' | \
            singularity exec --bind {params.bind} {params.sif} sed 's/\t1/\tM/g' | \
            singularity exec --bind {params.bind} {params.sif} sed 's/\t2/\tF/g' > {output}
        """



rule genotype_donor_annotation:
    input:
        het_psams = lambda wildcards: expand(output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.psam", ancestry = ancestry_file["unique_ancestry"]),
        updated_psam = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam"
    output:
        output_dict["output_dir"] + "/genotype_donor_annotation.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["sex4imputation_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * final_vcf_dict["sex4imputation_memory"]
    threads: final_vcf_dict["sex4imputation_threads"]
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
        singularity exec --bind {params.bind} {params.sif} sed -i 's/$/\tSangerImputationServer/' {output}
        singularity exec --bind {params.bind} {params.sif} sed -i '1s/SangerImputationServer/imputation_server/g' {output}
        singularity exec --bind {params.bind} {params.sif} rm {params.out_temp}
        singularity exec --bind {params.bind} {params.sif} rm {params.out_temp2}
        """
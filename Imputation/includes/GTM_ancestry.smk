#!/usr/local/envs/py36/bin python3


rule gtm_preprocess:
    input:
        mhc = input_dict["pipeline_dir"] + "/MHC_location.txt",
        bed = output_dict["output_dir"] + "/hrc_check/grm_subset-updated-chr{chr}.bed"
    output:
        output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}.prune.in"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_preprocess_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_preprocess_memory"]
    threads: GTM_check_dict["gtm_preprocess_threads"]
    params:
        bind = input_dict["bind_paths"],
        bed = output_dict["output_dir"] + "/hrc_check/grm_subset-updated-chr{chr}",
        out = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}",
        hrc_check = output_dict["output_dir"] + "/hrc_check",
        sif = input_dict["singularity_image"],
        out_base = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess"
    shell:
        """
        if [[ {wildcards.chr} == 6 ]]
        then
            singularity exec --bind {params.bind} {params.sif} plink --exclude {input.mhc} --bfile {params.bed} --indep-pairwise 1000 10 0.02 --maf 0.05 --out {params.out} --make-bed
        else
            singularity exec --bind {params.bind} {params.sif} plink --bfile {params.bed} --indep-pairwise 1000 10 0.02 --maf 0.05 --out {params.out} --make-bed
        fi
        """

rule gtm_prune:
    input:
        bed = output_dict["output_dir"] + "/hrc_check/grm_subset-updated-chr{chr}.bed",
        prune = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}.prune.in"
    output:
        bed = output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * GTM_check_dict["gtm_prune_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * GTM_check_dict["gtm_prune_memory"]
    threads: GTM_check_dict["gtm_prune_threads"]
    params:
        bind = input_dict["bind_paths"],
        bed = output_dict["output_dir"] + "/hrc_check/grm_subset-updated-chr{chr}",
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.bed} --extract {input.prune} --out {params.out} --make-bed
        """


rule gtm_merge:
    input:
        expand(output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}.bed", chr = range(1, 24))
    output:
        merge = output_dict["output_dir"] + "/gtm_merge/tomerge.txt",
        out = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.raw",
        mat = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.mat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_merge_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_merge_memory"]
    threads: GTM_check_dict["gtm_merge_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr",
        in_file = output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr"
    shell:
        """
        echo "starting"
        for chrom in {{1..23}}
        do
            echo $chrom
            singularity exec --bind {params.bind} {params.sif} echo {params.in_file}$chrom.bed | singularity exec --bind {params.bind} {params.sif} sed 's/.bed//g' >> {output.merge}
        done
            singularity exec --bind {params.bind} {params.sif} plink --merge-list {output.merge} --recodeA --out {params.out}
            singularity exec --bind {params.bind} {params.sif} tail -n +2 {output.out} > {output.mat}
        """

rule gtm_projection:
    input:
        mat = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.mat"
    output:
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_indiv_predictions.csv",
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_indiv_probabilities.csv",
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_group_probabilities.csv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_projection_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_projection_memory"]
    threads: GTM_check_dict["gtm_projection_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out_base = output_dict["output_dir"] + "/gtm_projection/gtm_projection"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python /opt/ancestry_viz/runGTM.py \
            --model GTM \
            --test /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-pipeline-QC/Imputation/recoded_1000G.noadmixed_test.mat \
            --data {input.mat} \
            --labels /opt/ancestry_viz/recoded_1000G.raw.noadmixed.lbls3_3 \
            --labeltype discrete \
            --out out_base \
            --pca \
            --n_components 10 \
            --missing \
            --missing_strategy median \
            --random_state 8 
        """

            # --test {input.mat} \
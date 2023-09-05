#!/usr/bin/env python
import pandas as pd

# Get list of pools to process
samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")
samples.columns = ["Pool", "N"]
samples["Pool"] = samples["Pool"].astype(str)

####################################
############ SOUPORCELL ############
####################################
rule souporcell_unzip_barcodes:
    input:
        barcodes = lambda wildcards: scrnaseq_libs_df["Barcode_Files"][wildcards.pool]
    threads: 1
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 5
    output:
        output_dict["output_dir"] + "/{pool}/souporcell/barcodes.tsv"
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"]
    log: output_dict["output_dir"] + "/logs/souporcell_unzip_barcodes.{pool}.log"
    shell:
        """
        if [[ {input.barcodes} == *".gz"* ]]
        then
            singularity exec --bind {params.bind} {params.sif} gunzip < {input.barcodes} > {output}
        else 
            singularity exec --bind {params.bind} {params.sif} cp {input.barcodes} {output}
        fi
        """

rule souporcell:
    input:
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
        barcodes = output_dict["output_dir"] + "/{pool}/souporcell/barcodes.tsv",
        fasta = fasta,
        common_variants = input_dict["snp_genotypes_filepath"]
    threads: souporcell_dict["souporcell_threads"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"]
    output:
        clusters = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv",
        genotypes = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf",
        troublet = output_dict["output_dir"] + "/{pool}/souporcell/troublet.done",
        concensus = output_dict["output_dir"] + "/{pool}/souporcell/consensus.done",
        clustering = output_dict["output_dir"] + "/{pool}/souporcell/clustering.done"
    params:
        out = output_dict["output_dir"] + "/{pool}/souporcell/",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        clusters = lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0],
        min_alt = souporcell_extra_dict["min_alt"],
        min_ref = souporcell_extra_dict["min_ref"],
        max_loci = souporcell_extra_dict["max_loci"] 
    log: output_dict["output_dir"] + "/logs/souporcell.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} souporcell_pipeline.py \
            --bam {input.bam} \
            --barcodes {input.barcodes} \
            --fasta {input.fasta} \
            --threads {threads} \
            --out_dir {params.out} \
            --clusters {params.clusters} \
            --common_variants {input.common_variants} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} 2> {log}
        [[ -s {output.genotypes} ]]
        echo $?
        """

#####################################################
############ REFORMAT SOUPORCELL RESULTS ############
#####################################################
rule souporcell_results_temp:
    input:
        souporcell = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv"
    output:
        output_dict["output_dir"] + "/{pool}/CombinedResults/souporcell_results.txt"
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"]
    log: output_dict["output_dir"] + "/logs/souporcell_results_temp.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} \
            awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,$4,$5}}' {input.souporcell} | \
            awk 'BEGIN{{FS=OFS="\t"}} $2=="doublet" {{$3="doublet"}}1' | \
            awk 'BEGIN{{FS=OFS="\t"}} $2=="unassigned" {{$4="unassigned"}}1' | \
            sed "s/status/DropletType/g" | sed "s/assignment/Assignment/g" | \
            sed "s/log_prob_singleton/LogProbSinglet/g" | \
            sed "s/log_prob_doublet/LogProbDoublet/g" | \
            sed "s/barcode/Barcode/g" | \
            sed "1s/\t/\tsouporcell_/g" | \
            awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output} 2> {log}
        """

#####################################
############ SUBSET VCFS ############
#####################################
rule souporcell_pool_vcf:
    input:
        genotypes = input_dict["snp_genotypes_filepath"], 
        cluster_geno = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf"
    output:
        filtered_refs = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        filtered_refs_temp = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        individuals = lambda wildcards: scrnaseq_libs_df["Individuals_Files"][wildcards.pool]
    log: output_dict["output_dir"] + "/logs/souporcell_pool_vcf.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools intersect -a {input.genotypes} -b {input.cluster_geno} -f 1.0 -r -wa -header > {output.filtered_refs_temp} 2> {log}
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {params.individuals} -Oz -o {output.filtered_refs} {output.filtered_refs_temp} 2>> {log}
        """


###############################################################################
############ CORRELATE INDIVIDUAL GENOTYPES WITH CLUSTER GENOTYPES ############
###############################################################################
## To take in souporcell_pool_vcf output and souporcell cluster vcf output
rule souporcell_correlate_genotypes:
    input:
        genotypes = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        results_file = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    output:
        assignments = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv",
        correlation = report(output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations/pearson_correlation.png", category = "Souporcell Genotype Correlations", subcategory = "{pool}", caption = "../report_captions/souporcell.rst")
    resources:
        mem_per_thread_gb = souporcell_dict["souporcell_correlations_memory"],
        disk_per_thread_gb = souporcell_dict["souporcell_correlations_memory"]
    threads: souporcell_dict["souporcell_correlations_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        basedir = output_dict["output_dir"],
        cor_thresh = souporcell_dict["souporcell_genotype_correlation_threshold"]
    log: output_dict["output_dir"] + "/logs/souporcell_correlate_genotypes.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --basedir {params.out} \
            --pool {wildcards.pool} \
            --result_file {input.results_file} \
            --correlation_limit {params.cor_thresh} \
            2> {log}
        [[ -s {output.assignments} ]]
        echo $?
        """

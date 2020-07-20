#!/usr/bin/env python
import os
import pandas as pd
from glob import glob

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["refs"]
souporcell_dict = config["souporcell"]
souporcell_extra_dict = config["souporcell_extra"]


####################################
############ SOUPORCELL ############
####################################
rule souporcell:
    input:
        bam = ,
        barcodes = ,
        fasta = ref_dict["fasta_filepath"]
        snps = input_dict["snp_genotypes_filepath"]
    threads: souporcell_dict["souporcell_threads"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"]
    output:
        clusters = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv"
        genotypes = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf"
    params:
        out = output_dict["output_dir"] + "/{pool}/souporcell/",
        sif = input_dict["singularity_image"],
        N = lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
        min_alt = souporcell_extra_dict["min_alt"] 
        min_ref = souporcell_extra_dict["min_ref"] 
        max_loci = souporcell_extra_dict["max_loci"] 
    shell:
        """
        singularity exec {params.sif} souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            -t {threads} \s
            -o {params.out} \
            -k {params.N} \
            --common_variants {input.snps} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} \
        [[ -s {output.genotypes} ]]
        [[ -s {output.clusters} ]]
        echo $?
        """

#####################################################
############ REFORMAT SOUPORCELL RESULTS ############
#####################################################
rule souporcell_results_temp:
    input:
        souporcell = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv"
    output:
        souporcell_temp = temp(output_dict["output_dir"] + "/{pool}/CombinedResults/souporcell_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,$4,$5}}' {input.souporcell} | \
            singularity exec {params.sif} awk 'BEGIN{{FS=OFS="\t"}} $2=="doublet" {{$3="doublet"}}1' | \
            singularity exec {params.sif} awk 'BEGIN{{FS=OFS="\t"}} $2=="unassigned" {{$4="unassigned"}}1' | \
            singularity exec {params.sif} sed "s/status/DropletType/g" | sed "s/assignment/Assignment/g" | \
            singularity exec {params.sif} sed "s/log_prob_singleton/LogProbSinglet/g" | \
            singularity exec {params.sif} sed "s/log_prob_doublet/LogProbDoublet/g" | \
            singularity exec {params.sif} sed "s/barcode/Barcode/g" | \
            singularity exec {params.sif} sed "1s/\t/\tsouporcell_/g" | \
            singularity exec {params.sif} awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output.souporcell_temp}
        """

#####################################
############ SUBSET VCFS ############
#####################################
rule souporcell_pool_vcf:
    input:
        genotypes = input_dict["snp_genotypes_filepath"], 
        cluster_geno = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf"
    output:
        output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        sif = input_dict["singularity_image"],
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    shell:
        "singularity exec {params.sif} bcftools view -s {params.individuals} -R {input.cluster_geno} -Oz -o {output} {input.genotypes}"


###############################################################################
############ CORRELATE INDIVIDUAL GENOTYPES WITH CLUSTER GENOTYPES ############
###############################################################################
## To take in souporcell_pool_vcf output and souporcell cluster vcf output
rule souporcell_correlate_genotypes:
    input:
        genotypes = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        assignments = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    output:
        assignments = output_dict["output_dir"] + "/{pool}/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"
        variables = temp(output_dict["output_dir"] + "/{pool}/souporcell/souporcel_genotypes_variables")
    resorces:
        mem_per_thread_gb = souporcell_dict["souporcell_correlations_memory"],
        disk_per_thread_gb = souporcell_dict["souporcell_correlations_memory"]
    threads: souporcell_dict["souporcell_correlations_threads"]
    params:
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"],
        script = input_dict["pipeline_dir"] + "/scripts/Assign_Indiv_by_Geno.R"
    shell:
        """
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {wildcards.pool} >> {output.variables}
        singularity exec {params.sif} echo {input.assignments} >> {output.variables}
        singularity exec {params.sif} Rscript {params.script} {output.variables}
        [[ -s {output.assignments} ]]
        echo $?
        """

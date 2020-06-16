#!/usr/local/envs/py36/bin python3

import os 
import pandas as pd

#### The commented section here is what we would like to implement - the ability to provide additional global environments from the command line instead of editing the file
# envvars:
#     "DATADIR",
#     "OUTDIR",
#     "SAMPLE_FILE"

# datadir = os.environ["DATADIR"]
# outdir = os.environ["OUTDIR"]

samples = pd.read_csv("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Pool_Name_N_Individuals_temp.txt", sep = "\t")
datadir = "/directflow/SCCGGroupShare/projects/data/experimental_data/CLEAN/OneK1K_scRNA/OneK1K_scRNA_V1"
outdir = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium"
FASTA="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa"
FAI="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa.fai"
SNP_GENOTYPES="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Imputed/Merged_MAF0.01.dose_GeneFiltered_hg38_nochr_NoAdditionalchr.vcf" #GENES only
SNP_GENOTYPES_LIST="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Imputed/ImputedGenotypeLocations.tsv"
SNVs_list="/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh38SNPvcfs1000genomes/MergedAutosomesFilteredGenes.recode.MAF0.01.vcf"

T = 8

rule all:
    input:
        expand(outdir + "/{pool}/popscle/demuxlet/demuxletOUT.best",  pool=samples.Pool),
        expand(outdir + "/{pool}/souporcell/cluster_genotypes.vcf", pool=samples.Pool),
        expand(outdir + "/{pool}/souporcell/clusters.tsv", pool=samples.Pool),
        expand(outdir +  "/{pool}/CombinedResults/CombinedDropletAssignments.tsv", pool=samples.Pool),
        expand(outdir + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz", pool=samples.Pool)


###################################
############# POPSCLE #############
###################################

###### popscle Preprocessing ######
rule popscle_pileup:
    input:
        vcf=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    output:
        directory(outdir + "/{pool}/popscle/pileup/")
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 10
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        singularity exec {params.sif} popscle dsc-pileup --sam {input.bam} --vcf {input.vcf} --group-list {input.barcodes} --out {output}pileup
        [[ -s {output}pileup.var.gz ]]
        echo $?
        """

##### Popscle Demuxlet Individual File Generation #####
rule popscle_demuxlet_ind_files:
    input:
        outdir + "/{pool}/popscle/pileup/"
    output:
        outdir + "/{pool}/popscle/Individuals.txt"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} echo {params.individuals} | tr "," "\n" > {output}
        [[ -s {output} ]]
        echo $?
        """

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        snps=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        individuals=outdir + "/{pool}/popscle/Individuals.txt"
    output:
        outdir + "/{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 25,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 25
    threads: 4
    params:
        out=outdir + "/{pool}/popscle/demuxlet/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        field="GP"
    shell:
        """
        singularity exec {params.sif} popscle demuxlet --plp {input.pileup}pileup --vcf {input.snps} --field {params.field} --group-list {input.barcodes} --out {params.out}demuxletOUT --sm-list {input.individuals}
        [[ -s {output} ]]
        echo $?
        """

####################################
############ SOUPORCELL ############
####################################
rule souporcell:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        fasta=FASTA,
        snps=SNP_GENOTYPES
    threads: T
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    output:
        clusters=outdir + "/{pool}/souporcell/clusters.tsv"
    params:
        out=outdir + "/{pool}/souporcell/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} souporcell_pipeline.py -i {input.bam} -b {input.barcodes} -f {input.fasta} -t {threads} -o {params.out} -k {params.N} --common_variants {input.snps}
        [[ -s {output.clusters} ]]
        echo $?
        """

##################################
############ SCRUBLET ############
##################################
rule scrublet:
    input:
        matrix=outdir + "/{pool}/matrix_out/matrix.mtx",
        genes=outdir + "/{pool}/matrix_out/features.tsv"
    output:
        outdir + "/{pool}/scrublet_{pctl}/predicted_doublet_mask.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
    threads: 1
    params:
        pctl="{pctl}",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/DoubletDetection.sif",
        out=outdir + "/{pool}/scrublet_{pctl}/"
    benchmark:
        outdir + "benchmarks/{pool}.scrublet_{pctl}.txt"
    shell:
        """
        singularity exec {params.sif} python /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/scripts/ONEK1K/hg38/refSNVs/sceQTL-Gen-Demultiplex/scrublet_pipeline.py --counts_matrix {input.matrix} --genes {input.genes} --min_gene_variability_pctl {params.pctl} -o {params.out}
        [[ -s {output} ]]
        echo $?
        """

##############################
############ SCDS ############
##############################
rule scds:
    input:
        script="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/scripts/ONEK1K/hg38/TestTxnDoublets/scds.R",
        genes=outdir + "/{pool}/matrix_out/genes.tsv"
    output: 
        doublets=outdir + "/{pool}/scds/scds_doublets.txt",
        variables=outdir + "/{pool}/scds/scds_variables.txt",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 1
    params:
        matrix_dir=outdir + "/{pool}/matrix_out/",
        out=outdir + "/{pool}/scds/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/DoubletDetection.sif"
    benchmark:
        outdir + "benchmarks/{pool}.scds.txt"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} Rscript {input.script} {output.variables}
        [[ -s {output.doublets} ]]
        echo $?
        """


###########################################
############ DOUBLET DETECTION ############
###########################################
rule DoubletDetection:
    input:
        script="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/scripts/ONEK1K/hg38/TestTxnDoublets/DoubletDetection.py",
        matrix=outdir + "/{pool}/matrix_out/matrix.mtx.gz"
    output:
        doublets = outdir + "/{pool}/DoubletDetection/DoubletDetection_predicted_doublet.txt",
        variables=outdir + "/{pool}/DoubletDetection/DoubletDetection_variables.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 2
    params:
        matrix_dir=outdir + "/{pool}/matrix_out/",
        out=outdir + "/{pool}/DoubletDetection/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/DoubletDetection.sif"
    benchmark:
        outdir + "benchmarks/{pool}.DoubletDetection.txt"
    shell:
        """
        singularity exec {params.sif} echo {wildcards.pool} > {output.variables}
        singularity exec {params.sif} echo {params.out} >> {output.variables}
        singularity exec {params.sif} echo {params.matrix_dir} >> {output.variables}
        singularity exec {params.sif} python {input.script} --counts_matrix {input.matrix} -o {params.out}
        [[ -s {output.doublets} ]]
        echo $?
        """

#################################
######## COMBINE RESULTS ########
#################################
rule demuxlet_results_temp:
    input:
        demuxlet=outdir + "/{pool}/popscle/demuxlet/demuxletOUT_impute_vars.best"
    output:
        demuxlet_temp=temp(outdir + "/{pool}/CombinedResults/demuxlet_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$6,$14,$19,$20}}' {input.demuxlet} | sed "s/SNG/singlet/g" | sed "s/DBL/doublet/g" | awk 'BEGIN{{FS=OFS="\t"}} $3=="doublet" {{$4="doublet"}}1' | sed -E "s/,[0-9]+_[0-9]+,[0-9].[0-9]+\t/\t/g" | sed "s/NUM.SNPS/nSNP/g" | sed "s/DROPLET.TYPE/DropletType/g" | sed "s/BEST.GUESS/Assignment/g" | sed "s/singlet.BEST.LLK/SingletLLK/g" | sed "s/doublet.BEST.LLK/DoulbetLLK/g" | sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | sed "1s/\t/\tdemuxlet_/g" | sed "s/BARCODE/Barcode/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}'  > {output.demuxlet_temp}
        """

rule souporcell_results_temp:
    input:
        souporcell=outdir + "/{pool}/souporcell/clusters.tsv"
    output:
        souporcell_temp=temp(outdir + "/{pool}/CombinedResults/souporcell_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,$4,$5}}' {input.souporcell} | awk 'BEGIN{{FS=OFS="\t"}} $2=="doublet" {{$3="doublet"}}1' | awk 'BEGIN{{FS=OFS="\t"}} $2=="unassigned" {{$4="unassigned"}}1' | sed "s/status/DropletType/g" | sed "s/assignment/Assignment/g" | sed "s/log_prob_singleton/LogProbSinglet/g" | sed "s/log_prob_doublet/LogProbDoublet/g" | sed "s/barcode/Barcode/g" | sed "1s/\t/\tsouporcell_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output.souporcell_temp}
        """

rule join_results:
    input:
        demuxlet=outdir + "/{pool}/CombinedResults/demuxlet_temp.txt",
        souporcell=outdir + "/{pool}/CombinedResults/souporcell_temp.txt",
    output:
        outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    shell:
        """
         join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7" {input.demuxlet} {input.freemuxlet} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3" - {input.scSplit} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2,2.3,2.4,2.5" - {input.souporcell} | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' | \
            join -a1 -a2 -1 1 -2 1 -t "\t" -e"NA" -o "0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.2,2.3,2.4,2.5,2.6" - {input.vireo} > {output}
        """

#####################################################################
############ SCRIPT FOR CALLING FINAL BARCODE ASSIGNMENT ############
#####################################################################

#####################################
############ SUBSET VCFS ############
#####################################
rule souporcell_pool_vcf:
    input:
        genotypes=SNP_GENOTYPES + ".gz",
        cluster_geno=outdir + "/{pool}/souporcell/cluster_genotypes.vcf"
    output:
        outdir + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
        mem_per_thread_gb=5,
        disk_per_thread_gb=5
    threads: 1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    shell:
        "singularity exec {params.sif} bcftools view -s {params.individuals} -R {input.cluster_geno} -Oz -o {output} {input.genotypes}"

###############################################################################
############ CORRELATE INDIVIDUAL GENOTYPES WITH CLUSTER GENOTYPES ############
###############################################################################
## To take in souporcell_pool_vcf output and souporcell cluster vcf output



###############################################################################
############ SCRIPT TO PRODUCE QC PLOTS ############
###############################################################################
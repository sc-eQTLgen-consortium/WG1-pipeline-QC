#!/usr/local/envs/py36/bin python3

import os 
import pandas as pd

# envvars:
#     "DATADIR",
#     "OUTDIR",
#     "SAMPLE_FILE"

# datadir = os.environ["DATADIR"]
# outdir = os.environ["OUTDIR"]

samples = pd.read_csv("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Pool_Name_N_Individuals_temp.txt", sep = "\t")
# samples = pd.read_csv("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Pool_Name_N_Individuals_test.txt", sep = "\t")
datadir = "/directflow/SCCGGroupShare/projects/data/experimental_data/CLEAN/OneK1K_scRNA/OneK1K_scRNA_V1"
outdir = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium"
FASTA="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa"
FAI="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa.fai"
SNP_GENOTYPES="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Imputed/Merged_MAF0.01.dose_GeneFiltered_hg38_nochr.vcf" #GENES only
SNVs_list="/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh38SNPvcfs1000genomes/MergedAutosomesFilteredGenes.recode.MAF0.01.vcf"

T = 8

rule all:
    input:
        expand(outdir + "/{pool}/popscle/demuxlet/",  pool=samples.Pool),
        expand(outdir + "/{pool}/popscle/freemuxlet/", pool=samples.Pool),
        expand(outdir + "/{pool}/souporcell/", pool=samples.Pool),
        expand(outdir + "/{pool}/vireo/results/", pool=samples.Pool),
        expand(outdir +  "/{pool}/scSplit/genotypes/", pool=samples.Pool)

rule scSplit_sam_header:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    threads: T
    output:
        temp(outdir + "/{pool}/scSplit/SAM_header")
    resources:
        mem_per_thread_gb=1
    shell:
        "samtools view -@ {threads} -H {input.bam} > {output}"

rule scSplit_sam_body:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    threads: T
    resources:
        mem_per_thread_gb=1
    output:
        temp(outdir + "/{pool}/scSplit/filtered_SAM_body")
    shell:
        "samtools view -@ {threads} -S -q 10 -F 3844 {input.bam} | LC_ALL=C grep -F -f {input.barcodes} > {output}"

rule scSplit_sam_combine:
    input:
        header=outdir + "/{pool}/scSplit/SAM_header",
        body=outdir + "/{pool}/scSplit/filtered_SAM_body"
    threads: T
    resources:
        mem_per_thread_gb=1
    output:
        temp(outdir + "/{pool}/scSplit/filtered.bam")
    shell:
        "cat {input.header} {input.body} | samtools view -@ {threads} -b - > {output}"

rule scSplit_rmdupe:
    input:
        bam=outdir + "/{pool}/scSplit/filtered.bam"
    output:
        temp(outdir + "/{pool}/scSplit/dedup_filtered.bam")
    resources:
        mem_gb=20
    shell:
        "samtools rmdup {input.bam} {output}"

rule scSplit_sort:
    input:
        outdir + "/{pool}/scSplit/dedup_filtered.bam"
    threads: T
    resources:
        mem_per_thread_gb=15
    output:
        outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        samtools index {output}
        """

##### scSplit Allele Counting #####
rule scSplit_allele_matrices:
    input:
        snvs=SNVs_list,
        vcf=SNP_GENOTYPES,
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    output:
        directory(outdir + "/{pool}/scSplit/allele_matrices/")
    params:
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    resources:
        mem_gb=25
    shell:
        "scSplit count -c {input.snvs} -v {input.vcf} -i {input.bam} -b {params.barcodes} -r {output}ref_filtered.csv -a {output}alt_filtered.csv -o {output}"

##### scSplit Demultiplexing #####
rule scSplit_demultiplex:
    input:
        outdir + "/{pool}/scSplit/allele_matrices/"
    output:
        directory(outdir + "/{pool}/scSplit/demultiplex/")
    params:
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    resources:
        mem_gb=25
    shell:
        "scSplit run -r {input}ref_filtered.csv -a ${input}alt_filtered.csv -n {params.N} -o {output}"

##### scSplit Get Genotypes #####
rule scSplit_genotypes:
    input:
        matrices=outdir + "/{pool}/scSplit/allele_matrices/",
        demultiplex=outdir + "/{pool}/scSplit/demultiplex/"
    output:
        directory(outdir + "/{pool}/scSplit/genotypes/")
    resources:
        mem_gb=25
    shell:
        "scSplit genotype -r {input.matrices}ref_filtered.csv -a {input.matrices}alt_filtered.csv -p {input.demultiplex}scSplit_P_s_c.csv -o {output}"

##### Popscle Pileup #####
rule popscle_pileup:
    input:
        vcf=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    output:
        directory(outdir + "/{pool}/popscle/pileup/")
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 40
    shell:
        "popscle dsc-pileup --sam {input.bam} --vcf {input.vcf} --group-list {input.barcodes} --out {output}pileup"

##### Popscle Freemuxlet Demultiplexing #####
rule popscle_freemuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    output:
        directory(outdir + "/{pool}/popscle/freemuxlet/")
    params:
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    resources:
        mem_gb=ambda wildcards, attempt: attempt * 20
    shell:
        "popscle freemuxlet --plp {input.pileup}pileup --out {output}freemuxletOUT --group-list {input.barcodes} --nsample {params.N}"

##### Popscle Demuxlet Individual File Generation #####
rule popscle_demuxlet_ind_files:
    input:
        expand(outdir + "/{pool}/popscle/pileup/",  pool=samples.Pool)
    output:
        temp(outdir + "/{pool}/popscle/Individuals.txt")
    params:
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    resources:
        mem_gb=5
    shell:
        "echo {params.individuals} | sed 's/,/\n/g' > {output}; ls {input}"

##### Popscle Demuxlet Demultiplexing #####
rule popscle_demuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        snps=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        individuals=outdir + "/{pool}/popscle/Individuals.txt"
    output:
        directory(outdir + "/{pool}/popscle/demuxlet/")
    params:
        field="GP"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 50
    shell:
        "popscle demuxlet --plp {input.pileup}pileup --vcf {input.snps} --field {params.field} --group-list {input.barcodes} --out {output}demuxletOUT --sm-list {input.individuals}"

##### cellSNP Pileup #####
rule cellSNP:
    input:
        vcf=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    output:
        outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz"
    params:
        p=20,
        maf=0.1,
        count=20
    resources:
        mem_gb=40
    shell:
        "cellSNP -s {input.bam} -b {input.barcodes} -o {output} -R {input.vcf} -p {params.p} --minMAF {params.maf} --minCOUNT {params.count}"

##### Subset the imputed genotype files by the individuals in the pools #####
rule subset_vcf:
    input:
        pileup=outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz",
        snps=SNP_GENOTYPES,
    output:
        outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf"
    params:
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    resources:
        mem_gb=10
    shell:
        "bcftools view -R {input.pileup} -s {params.individuals} -Oz -o {output} {input.snps}"

##### Vireo demultiplexing #####
rule vireo:
    input:
        pileup=outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz",
        snps=outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf"
    output:
        outdir + "/{pool}/vireo/results/"
    params:
        field="GP"
    resources:
        mem_gb=55
    shell:
        "vireo -c {input.pileup} -d {input.snps} -o {outdir} -t {params.field}"


##### Run the souporcell pipeline #####
rule souporcell:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        fasta=FASTA,
        snps=SNP_GENOTYPES
    threads: T
    resources:
        mem_gb=40
    output:
        directory(outdir + "/{pool}/souporcell/")
    params:
        individuals=lambda wildcards: str(samples.Individuals[samples.Pool == wildcards.pool].iloc[0]).replace(",", " "),
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        "souporcell_pipeline.py -i {input.bam} -b {input.barcodes} -f {input.fasta} -t {threads} -o {output} -k {params.N} --known_genotypes {input.snps} --known_genotypes_sample_names {params.individuals}"


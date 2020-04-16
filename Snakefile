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
SNP_GENOTYPES="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Imputed/Merged_MAF0.01.dose_GeneFiltered_hg38_nochr_NoAdditionalchr.vcf" #GENES only
SNVs_list="/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh38SNPvcfs1000genomes/MergedAutosomesFilteredGenes.recode.MAF0.01.vcf"

T = 8

rule all:
    input:
        expand(outdir + "/{pool}/popscle/demuxlet/demuxletOUT.best",  pool=samples.Pool),
        expand(outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz", pool=samples.Pool),
        expand(outdir + "/{pool}/souporcell/cluster_genotypes.vcf", pool=samples.Pool),
        expand(outdir + "/{pool}/vireo/results/donor_ids.tsv", pool=samples.Pool),
        expand(outdir +  "/{pool}/scSplit/scSplit.vcf", pool=samples.Pool)

rule scSplit_sam_header:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    threads: T
    output:
        temp(outdir + "/{pool}/scSplit/SAM_header")
    resources:
        mem_per_thread_gb=1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "scSplit"
    shell:
        "singularity exec {params.sif} samtools view -@ {threads} -H {input.bam} > {output}"

rule scSplit_sam_body:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    threads: T
    resources:
        mem_per_thread_gb=1
    output:
        temp(outdir + "/{pool}/scSplit/filtered_SAM_body")
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "scSplit"
    shell:
        "singularity exec {params.sif} samtools view -@ {threads} -S -q 10 -F 3844 {input.bam} | LC_ALL=C grep -F -f {input.barcodes} > {output}"

rule scSplit_sam_combine:
    input:
        header=outdir + "/{pool}/scSplit/SAM_header",
        body=outdir + "/{pool}/scSplit/filtered_SAM_body"
    threads: T
    resources:
        mem_per_thread_gb=1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "scSplit"
    output:
        temp(outdir + "/{pool}/scSplit/filtered.bam")
    shell:
        "singularity exec {params.sif} cat {input.header} {input.body} | singularity exec {params.sif} samtools view -@ {threads} -b - > {output}"

rule scSplit_rmdupe:
    input:
        bam=outdir + "/{pool}/scSplit/filtered.bam"
    output:
        temp(outdir + "/{pool}/scSplit/dedup_filtered.bam")
    resources:
        mem_gb=20
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "scSplit"
    shell:
        "singularity exec {params.sif} samtools rmdup {input.bam} {output}"

rule scSplit_sort:
    input:
        outdir + "/{pool}/scSplit/dedup_filtered.bam"
    threads: T
    resources:
        mem_per_thread_gb=15
    output:
        outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "scSplit"
    shell:
        """
        singularity exec {params.sif} samtools sort -@ {threads} -o {output} {input}
        singularity exec {params.sif} samtools index {output}
        """

rule scSplit_regions:
    input:
        fai=FAI,
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    output:
        temp(outdir + "/{pool}/scSplit/regions_file")
    resources:
        mem_gb=5
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    group: "freebayes"
    shell:
        """
            singularity exec {params.sif} fasta_generate_regions.py {input.fai} 100000 > {output}
        """

rule scSplit_freebayes:
    input:
        fasta=FASTA,
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam",
        regions=outdir + "/{pool}/scSplit/regions_file"
    output:
        outdir + "/{pool}/scSplit/freebayes_var.vcf"
    group: "freebayes"
    threads: 40
    resources:
        mem_per_thread_gb=2
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
            export TMPDIR=/tmp
            singularity exec {params.sif} freebayes-parallel {input.regions} {threads} -f {input.fasta} -iXu -C 2 -q 1 {input.bam} > {output}
        """
 
rule scSplit_vcf_qual_filt:
    input:
        vcf=outdir + "/{pool}/scSplit/freebayes_var.vcf"
    output:
        outdir + "/{pool}/scSplit/frebayes_var_qual30.vcf"
    group: "freebayes"
    resources:
        mem_gb=10
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
            singularity exec {params.sif} vcftools --gzvcf {input.vcf} --minQ 30 --recode --recode-INFO-all --out {output}
        """      

##### scSplit Allele Counting #####
rule scSplit_allele_matrices:
    input:
        snvs=SNVs_list,
        vcf=outdir + "/{pool}/scSplit/frebayes_var_qual30.vcf",
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    output:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv"
    resources:
        mem_gb=120
    params:
        out=outdir + "/{pool}/scSplit/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} scSplit count -c {input.snvs} -v {input.vcf} -i {input.bam} -b {params.barcodes} -r {output.ref} -a {output.alt} -o {params.out}
        [[ -s {output.alt} ]]
        echo $?
        """

##### scSplit Demultiplexing #####
rule scSplit_demultiplex:
    input:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/alt_filtered.csv"
    output:
        outdir + "/{pool}/scSplit/scSplit_P_s_c.csv"
    resources:
        mem_gb=25
    params:
        out=outdir + "/{pool}/scSplit/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} scSplit run -r {input.ref} -a {input.alt} -n {params.N} -o {params.out}
        [[ -s {output} ]]
        echo $?
        """

##### scSplit Get Genotypes #####
rule scSplit_genotypes:
    input:
        matrices=outdir + "/{pool}/scSplit/",
        demultiplex=outdir + "/{pool}/scSplit/scSplit_P_s_c.csv"
    output:
        outdir + "/{pool}/scSplit/scSplit.vcf"
    resources:
        mem_gb=25
    params:
        out=outdir + "/{pool}/scSplit/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        scSplit genotype -r {input.matrices}ref_filtered.csv -a {input.matrices}alt_filtered.csv -p {input.demultiplex} -o {params.out}
        [[ -s {output} ]] 
        echo $?
        """

##### Popscle Pileup #####
rule popscle_pileup:
    input:
        vcf=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    output:
        directory(outdir + "/{pool}/popscle/pileup/")
    resources:
        mem_gb=250
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        singularity exec {params.sif} popscle dsc-pileup --sam {input.bam} --vcf {input.vcf} --group-list {input.barcodes} --out {output}pileup
        [[ -s {output}pileup.var.gz ]]
        echo $?
        """

##### Popscle Freemuxlet Demultiplexing #####
rule popscle_freemuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    output:
        outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz"
    resources:
        mem_gb=80
    params:
        out=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} popscle freemuxlet --plp {input.pileup}pileup --out {params.out} --group-list {input.barcodes} --nsample {params.N}
        [[ -s {output} ]]
        echo $?
        """

##### Popscle Demuxlet Individual File Generation #####
rule popscle_demuxlet_ind_files:
    input:
        outdir + "/{pool}/popscle/pileup/"
    output:
        temp(outdir + "/{pool}/popscle/Individuals.txt")
    resources:
        mem_gb=5
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
        mem_gb=500
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

##### cellSNP Pileup #####
rule cellSNP:
    input:
        vcf=SNP_GENOTYPES,
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    output:
        outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz"
    resources:
        mem_gb=40
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        p=20,
        maf=0.1,
        count=20
    shell:
        "singularity exec {params.sif} cellSNP -s {input.bam} -b {input.barcodes} -o {output} -R {input.vcf} -p {params.p} --minMAF {params.maf} --minCOUNT {params.count}"

##### Subset the imputed genotype files by the individuals in the pools #####
rule subset_vcf:
    input:
        pileup=outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz",
        snps=SNP_GENOTYPES + ".gz"
    output:
        outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf.gz"
    resources:
        mem_gb=10
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        individuals=lambda wildcards: samples.Individuals[samples.Pool == wildcards.pool].iloc[0]
    shell:
        "singularity exec {params.sif} bcftools view -R {input.pileup} -s {params.individuals} -Oz -o {output} {input.snps}"

##### Vireo demultiplexing #####
rule vireo:
    input:
        pileup=outdir + "/{pool}/vireo/cellSNPpileup.vcf.gz",
        snps=outdir + "/{pool}/vireo/Merged_MAF0.01.dose_GeneFiltered_hg38_individualSubset.vcf.gz"
    output:
        outdir + "/{pool}/vireo/results/donor_ids.tsv"
    resources:
        mem_gb=120
    params:
        out=outdir + "/{pool}/vireo/results/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        field="GP",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} vireo -c {input.pileup} -d {input.snps} -o {params.out} -t {params.field} -N {params.N}
        [[ -s {output} ]]
        echo $?
        """


##### Run the souporcell pipeline #####
rule souporcell:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0],
        fasta=FASTA,
        snps=SNP_GENOTYPES
    threads: T
    resources:
        mem_gb=200
    output:
        outdir + "/{pool}/souporcell/cluster_genotypes.vcf"
    params:
        out=outdir + "/{pool}/souporcell/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        # individuals=lambda wildcards: str(samples.Individuals[samples.Pool == wildcards.pool].iloc[0]).replace(",", " "),
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} souporcell_pipeline.py -i {input.bam} -b {input.barcodes} -f {input.fasta} -t {threads} -o {params.out} -k {params.N}
        [[ -s {output} ]]
        echo $?
        """


rule join_results:
    input:
        demuxlet=outdir + "/{pool}/popscle/demuxlet/demuxletOUT.best",
        freemuxlet=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz",
        scSplit=outdir + "/{pool}/scSplit/scSplit_results.csv",
        souporcell=outdir + "/{pool}/souporcell/clusters.tsv",
        vireo=outdir + "/{pool}/souporcell/donor_ids.tsv"
    output:
        outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    params:
        out=outdir + "/{pool}/CombinedResults/"
    shell:
        """
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$2,\$3,\$5,\$6,\$14,\$19,\$20)}}" {input.demuxlet} | sed "1s/\t/\tdemulxet_/" > {params.out}demuxlet_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$2,\$3,\$5,\$6,\$14,\$19,\$20)}}" {input.freemuxlet} | sed "1s/\t/\tfreemulxet_/" > {params.out}freemuxlet_temp.txt
        sed "1s/\t/\tscSplit_/" {input.scSplit} > {params.out}scSplit_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$1,\$2,\$3,\$4,\$5)}}" {input.souporcell} | sed "1s/\t/\tsouporcell_/" > {params.out}souporcell_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$1,\$2,\$3,\$4,\$5)}}" {input.vireo} | sed "1s/\t/\tvireo_/" > {params.out}vireo_temp.txt

        join -a1 -a2 -1 1 -2 1 {params.out}demuxlet_temp.txt {params.out}freemuxlet_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}scSplit_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}souporcell_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}vireo_temp.txt > {output}
        """


rule join_snps:
    input:
        popscle=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz",
        scSplit=outdir + "/{pool}/scSplit/scSplit_results.csv",
        souporcell=outdir + "/{pool}/souporcell/clusters.tsv",
        vireo=outdir + "/{pool}/souporcell/donor_ids.tsv"
    output:
        outdir + "/{pool}/CombinedResults/CombinedDropletAssignments.tsv"
    params:
        out=outdir + "/{pool}/CombinedResults/"
    shell:
        """
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$2,\$3,\$5,\$6,\$14,\$19,\$20)}}" {input.demuxlet} | sed "1s/\t/\tdemulxet_/" > {params.out}demuxlet_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$2,\$3,\$5,\$6,\$14,\$19,\$20)}}" {input.freemuxlet} | sed "1s/\t/\tfreemulxet_/" > {params.out}freemuxlet_temp.txt
        sed "1s/\t/\tscSplit_/" {input.scSplit} > {params.out}scSplit_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$1,\$2,\$3,\$4,\$5)}}" {input.souporcell} | sed "1s/\t/\tsouporcell_/" > {params.out}souporcell_temp.txt
        awk "BEGIN{{OFS=FS="\t"}}{{print(\$1,\$2,\$3,\$4,\$5)}}" {input.vireo} | sed "1s/\t/\tvireo_/" > {params.out}vireo_temp.txt

        join -a1 -a2 -1 1 -2 1 {params.out}demuxlet_temp.txt {params.out}freemuxlet_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}scSplit_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}souporcell_temp.txt | \
            join - -a1 -a2 -1 1 -2 1 {params.out}vireo_temp.txt > {output}
        """
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
        expand(outdir +  "/{pool}/CombinedResults/CombinedDropletAssignments.tsv", pool=samples.Pool)


rule scSplit_sam_header:
    input:
        bam=datadir + "/{pool}_V1/outs/possorted_genome_bam.bam"
    threads: T
    output:
        temp(outdir + "/{pool}/scSplit/SAM_header")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
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
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
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
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
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
        mem_per_thread_gb=20,
        disk_per_thread_gb=20
    threads: 1
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
        mem_per_thread_gb=15,
        disk_per_thread_gb=15
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
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
    threads: 1
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
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 2,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 2
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
        outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf"
    group: "freebayes"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 1
    params:
        out=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        singularity exec {params.sif} vcftools --gzvcf {input.vcf} --minQ 30 --recode --recode-INFO-all --out {params.out}
        [[ -s {output} ]]
        echo $?
        """      

rule scSplit_bgzip:
    input:
        outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf"
    output:
        gz=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz",
        index=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz.tbi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        singularity exec {params.sif} bgzip  -c {input} > {output.gz}
        singularity exec {params.sif} tabix -p vcf {output.gz}
        [[ -s {output.index} ]]
        echo $?
        """

##### This is how it should be done but forgot before pileup #####
rule scSplit_subset_vcf:
    input:
        pileup=outdir + "/{pool}/scSplit/freebayes_var_qual30.vcf.recode.vcf.gz",
        snps=SNP_GENOTYPES + ".gz"
    output:
        outdir + "/{pool}/scSplit/frebayes_var_qual30_subset.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 1
    params:
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
    shell:
        "singularity exec {params.sif} bcftools view {input.pileup} -R {input.snps} -Ov -o {output}"


##### scSplit Allele Counting #####
rule scSplit_allele_matrices:
    input:
        snvs=SNVs_list,
        vcf=outdir + "/{pool}/scSplit/frebayes_var_qual30_subset.vcf",
        bam=outdir + "/{pool}/scSplit/possort_dedup_filtered.bam"
    output:
        alt=outdir + "/{pool}/scSplit/alt_filtered.csv",
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 50,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 50
    threads: 4
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
        ref=outdir + "/{pool}/scSplit/ref_filtered.csv"
    output:
        Psc=outdir + "/{pool}/scSplit/scSplit_P_s_c.csv",
        result=outdir + "/{pool}/scSplit/scSplit_result.csv"
    threads: 5
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 30,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 30
    params:
        out=outdir + "/{pool}/scSplit/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """  
        singularity exec {params.sif} scSplit run -r {input.ref} -a {input.alt} -n {params.N} -o {params.out}
        [[ -s {output.Psc} ]]
        [[ -s {output.result} ]]
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
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 1
    params:
        out=outdir + "/{pool}/scSplit/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif"
    shell:
        """
        singularity exec {params.sif} scSplit genotype -r {input.matrices}ref_filtered.csv -a {input.matrices}alt_filtered.csv -p {input.demultiplex} -o {params.out}
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

##### Popscle Freemuxlet Demultiplexing #####
rule popscle_freemuxlet:
    input:
        pileup=outdir + "/{pool}/popscle/pileup/",
        barcodes=lambda wildcards: samples.Barcodes[samples.Pool == wildcards.pool].iloc[0]
    output:
        outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 30,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 30
    threads: 1
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
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 5
    params:
        out=outdir + "/{pool}/popscle/demuxlet/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        field="GP"
    shell:
        """
        singularity exec {params.sif} popscle demuxlet --plp {input.pileup}pileup --vcf {input.snps} --field {params.field} --group-list {input.barcodes} --geno-error-coeff 1.0 --geno-error-offset 0.05 --out {params.out}demuxletOUT_impute_vars --sm-list {input.individuals}
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
        mem_per_thread_gb=40,
        disk_per_thread_gb=40
    threads: 1
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
        mem_per_thread_gb=10,
        disk_per_thread_gb=10
    threads: 1
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
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 40,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 40
    threads: 1
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
        vcf=SNP_GENOTYPES
    threads: T
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 10
    output:
        genotypes=outdir + "/{pool}/souporcell/cluster_genotypes.vcf",
        clusters=outdir + "/{pool}/souporcell/clusters.tsv"
    params:
        out=outdir + "/{pool}/souporcell/",
        sif="/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Singularity_Buckets/docker/AllSoftwares.sif",
        individuals=lambda wildcards: str(samples.Individuals[samples.Pool == wildcards.pool].iloc[0]).replace(",", " "),
        N=lambda wildcards: samples.N[samples.Pool == wildcards.pool].iloc[0]
    shell:
        """
        singularity exec {params.sif} souporcell_pipeline.py -i {input.bam} -b {input.barcodes} -f {input.fasta} -t {threads} -o {params.out} -k {params.N} --known_genotypes {input.vcf} --known_genotypes_sample_names {params.individuals} --skip_remap SKIP_REMAP
        [[ -s {output.clusters} ]]
        [[ -s {output.genotypes} ]]
        echo $?
        """

rule demuxlet_results_temp:
    input:
        demuxlet=outdir + "/{pool}/popscle/demuxlet/demuxletOUT.best"
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

rule freemuxlet_results_temp:
    input:
        freemuxlet=outdir + "/{pool}/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz"
    output:
        freemuxlet_temp=temp(outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    shell:
        """
        gunzip -c {input.freemuxlet} | awk 'BEGIN{{OFS=FS="\\t"}}{{print $2,$3,$5,$6,$14,$19,$20 }}' | sed "s/SNG/singlet/g" | sed "s/DBL/doublet/g" | awk 'BEGIN{{FS=OFS="\\t"}} $3=="doublet" {{$4="doublet"}}1' | sed -E "s/,[0-9]+\t/\t/g" | sed "s/NUM.SNPS/nSNP/g" | sed "s/DROPLET.TYPE/DropletType/g" | sed "s/BEST.GUESS/Assignment/g" | sed "s/singlet.BEST.LLK/SingletLLK/g" | sed "s/doublet.BEST.LLK/DoulbetLLK/g" | sed "s/DIFF.LLK.singlet.doublet/DiffLLK/g" | sed "s/BARCODE/Barcode/g" | sed "1s/\t/\tfreemuxlet_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output.freemuxlet_temp}
        """

rule scSplit_results_temp:
    input:
        scSplit=outdir + "/{pool}/scSplit/scSplit_result.csv"
    output:
        scSplit_temp=temp(outdir + "/{pool}/CombinedResults/scSplit_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    shell:
        """
        sed -E 's/\tDBL-[0-9]+/\tdoublet\tdoublet/g' {input.scSplit} | sed 's/SNG-/singlet\t/g' | sed 's/Cluster/DropletType\tAssignment/g' | sed "1s/\t/\tscSplit_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output.scSplit_temp}
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

rule vireo_results_temp:
    input:
        vireo=outdir + "/{pool}/vireo/results/donor_ids.tsv"
    output:
        vireo_temp=temp(outdir + "/{pool}/CombinedResults/vireo_temp.txt")
    resources:
        mem_per_thread_gb=1,
        disk_per_thread_gb=1
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$2,$3,$4,$5}}' {input.vireo} | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("donor[0-9]+","singlet",$3)}}1' | awk 'BEGIN{{FS=OFS="\\t"}}{{gsub("[0-9]+_[0-9]+","singlet",$3)}}1' | sed "s/donor_id\tdonor_id/Assignment\tDropletType/g" | sed "s/prob_max/ProbSinglet/g" | sed "s/prob_doublet/ProbDoublet/g" | sed "s/n_vars/nSNP/g" | sed "s/cell/Barcode/g" | sed "1s/\t/\tvireo_/g" | awk 'NR<2{{print $0;next}}{{print $0| "sort -k1"}}' > {output.vireo_temp}
        """

rule join_results:
    input:
        demuxlet=outdir + "/{pool}/CombinedResults/demuxlet_temp.txt",
        freemuxlet=outdir + "/{pool}/CombinedResults/freemuxlet_temp.txt",
        scSplit=outdir + "/{pool}/CombinedResults/scSplit_temp.txt",
        souporcell=outdir + "/{pool}/CombinedResults/souporcell_temp.txt",
        vireo=outdir + "/{pool}/CombinedResults/vireo_temp.txt"
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

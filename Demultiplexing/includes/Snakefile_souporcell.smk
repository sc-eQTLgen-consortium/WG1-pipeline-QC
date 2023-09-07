####################################
############ SOUPORCELL ############
####################################
rule souporcell:
    input:
        bam = lambda wildcards: scrnaseq_libs_df["BamFile"][wildcards.pool],
        barcodes = lambda wildcards: scrnaseq_libs_df["BarcodeFile"][wildcards.pool],
        fasta = ref_dict["ref_dir"] + "/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        vcf = input_dict["vcf"]
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_memory"]
    threads: souporcell_dict["souporcell_threads"]
    output:
        ambient_rna = output_dict["output_dir"] + "/{pool}/souporcell/ambient_rna.txt",
        clustering = output_dict["output_dir"] + "/{pool}/souporcell/clustering.done",
        genotypes = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf",
        clusters = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv",
        concensus = output_dict["output_dir"] + "/{pool}/souporcell/consensus.done",
        troublet = output_dict["output_dir"] + "/{pool}/souporcell/troublet.done"
    params:
        out = output_dict["output_dir"] + "/{pool}/souporcell/",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        clusters = lambda wildcards: souporcell_size_dict[wildcards.pool],
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
            --common_variants {input.vcf} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} 2> {log}
        [[ -s {output.genotypes} ]]
        echo $?
        """

#################################
############ SUMMARY ############
#################################
rule souporcell_summary:
    input:
        clusters = output_dict["output_dir"] + "/{pool}/souporcell/clusters.tsv"
    output:
        summary = output_dict["output_dir"] + "/{pool}/souporcell/souporcell_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_summary_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_summary_memory"]
    threads: souporcell_dict["souporcell_summary_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        individuals = lambda wildcards: scrnaseq_libs_df["IndividualFile"][wildcards.pool]
    log: output_dict["output_dir"] + "/logs/souporcell_summary.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} \
            awk 'BEGIN{FS=OFS="\t"} $2=="unassigned" {$3="unassigned"}1' $1 \
            | awk 'BEGIN{FS=OFS="\t"}{print $3}' \
            | sed -E 's|[0-9]+/[0-9]+|doublet|g' \
            | tail -n+2 \
            | sort \
            | uniq -c \
            | sed -E 's/^ +//g' \
            | sed 's/ /\t/g' \
            | sed '1 i\Assignment N\tClassification' \
            | awk 'BEGIN{FS=OFS="\t"}{print($2,$1)}' {input.clusters} > {output.out}
        """

#####################################
############ SUBSET VCFS ############
#####################################
rule souporcell_pool_vcf:
    input:
        vcf = input_dict["vcf"],
        cluster_vcf = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf"
    output:
        filtered_refs_temp = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf",
        filtered_refs = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_pool_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * souporcell_dict["souporcell_pool_vcf_memory"]
    threads: souporcell_dict["souporcell_pool_vcf_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        individuals = lambda wildcards: scrnaseq_libs_df["IndividualFile"][wildcards.pool]
    log: output_dict["output_dir"] + "/logs/souporcell_pool_vcf.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools intersect -a {input.vcf} -b {input.cluster_vcf} -f 1.0 -r -wa -header > {output.filtered_refs_temp} 2> {log}
        singularity exec --bind {params.bind} {params.sif} bcftools view -S {params.individuals} -Oz -o {output.filtered_refs} {output.filtered_refs_temp} 2>> {log}
        """


###############################################################################
############ CORRELATE INDIVIDUAL GENOTYPES WITH CLUSTER GENOTYPES ############
###############################################################################
rule souporcell_correlate_genotypes:
    input:
        reference_vcf = output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        cluster_vcf = output_dict["output_dir"] + "/{pool}/souporcell/cluster_genotypes.vcf",
    output:
        correlation_file = output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlations.tsv",
        correlation_img = output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlation.png",
        assignments = output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt"
    resources:
        mem_per_thread_gb = souporcell_dict["souporcell_correlations_memory"],
        disk_per_thread_gb = souporcell_dict["souporcell_correlations_memory"]
    threads: souporcell_dict["souporcell_correlations_threads"]
    params:
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        out = output_dict["output_dir"] + "/{pool}/souporcell/genotype_correlations"
    log: output_dict["output_dir"] + "/logs/souporcell_correlate_genotypes.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script}
            --reference_vcf {input.reference_vcf} \
            --cluster_vcf {input.cluster_vcf} \
            --out {params.out} \
            2> {log}
        [[ -s {output.assignments} ]]
        echo $?
        """

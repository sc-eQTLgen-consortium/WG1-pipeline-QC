#!/usr/bin/env python
import json

#########################################
############# PREPROCESSING #############
#########################################


# In case of multiple inputs
rule combine_vcfs_all:
    input:
        vcfs = config["inputs"]["vcf"],
    output:
        vcf = config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38.vcf.gz",
        index = config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38.vcf.gz.csi"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["combine_vcfs_all_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["combine_vcfs_all_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["combine_vcfs_all_time"]]
    threads: config["demultiplex_preprocessing"]["combine_vcfs_all_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
    log: config["outputs"]["output_dir"] + "log/combine_vcfs_all.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools merge -Oz {input.vcfs} > {output.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools index {output.vcf}
        """


# Add all the info fields
# Filter the Imputed SNP Genotype by Minor Allele Frequency (MAF) and INFO scores
# TODO: maybe this fails if the input VCF is not gzipped?
rule filter4demultiplexing:
    input:
        vcf = config["inputs"]["vcf"][0] if len(config["inputs"]["vcf"]) == 1 else config["outputs"]["output_dir"] + "vcf_all_merged/imputed_hg38.vcf.gz",
    output:
        info_filled = temp(config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38_info_filled.vcf.gz"),
        qc_filtered = temp(config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38_qc_filtered.vcf.gz"),
        location_filtered = temp(config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38_qc_filtered_exons.recode.vcf.gz"),
        complete_cases = config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38_qc_filtered_exons_complete_cases.recode.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["filter4demultiplexing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["filter4demultiplexing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["filter4demultiplexing_time"]]
    threads: config["demultiplex_preprocessing"]["filter4demultiplexing_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        maf = config["demultiplex_preprocessing_extra"]["filter4demultiplexing_maf"],
        r2 = config["demultiplex_preprocessing_extra"]["filter4demultiplexing_r2"],
        bed = config["refs"]["ref_dir"] + config["refs_extra"]["relative_hg38_exons_ucsc_bed_path"]
    log: config["outputs"]["output_dir"] + "log/filter4demultiplexing.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bcftools +fill-tags -Oz --output {output.info_filled} {input.vcf}
        singularity exec --bind {params.bind} {params.sif} bcftools filter --include 'MAF>={params.maf} & R2>={params.r2}' -Oz --output {output.qc_filtered} {output.info_filled}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --gzvcf {output.qc_filtered} \
            --max-alleles 2 \
            --remove-indels \
            --bed {params.bed} \
            --recode \
            --recode-INFO-all \
            --stdout | gzip -c > {output.location_filtered}
        singularity exec --bind {params.bind} {params.sif} vcftools \
            --recode \
            --recode-INFO-all \
            --gzvcf {output.location_filtered} \
            --max-missing 1 \
            --stdout | gzip -c > {output.complete_cases}
        """


rule sort4demultiplexing:
    input:
        complete_cases = config["outputs"]["output_dir"] + "genotypes/vcf_all_merged/imputed_hg38_qc_filtered_exons_complete_cases.recode.vcf.gz"
    output:
        complete_cases_sorted = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz"
    resources:
        java_mem_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_memory"] * config["demultiplex_preprocessing"]["sort4demultiplexing_threads"] - config["settings_extra"]["java_memory_buffer"],
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["demultiplex_preprocessing"]["sort4demultiplexing_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["demultiplex_preprocessing"]["sort4demultiplexing_time"]]
    threads: config["demultiplex_preprocessing"]["sort4demultiplexing_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        jar = "/opt/picard-3.1.0/build/libs/picard.jar"
    log: config["outputs"]["output_dir"] + "log/sort4demultiplexing.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.java_mem_gb}g -Xms{resources.java_mem_gb}g -jar {params.jar} SortVcf \
            I={input.complete_cases} \
            O={output.complete_cases_sorted}
        """


###################################
############# POPSCLE #############
###################################


# TODO, gives error if barcodes is gzipped?
rule popscle_bam_filter:
    input:
        vcf = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz",
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"],
        bam = lambda wildcards: POOL_DF.loc[wildcards.pool, "Bam"]
    output:
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/snpfiltered_alignment.bam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_bam_filter_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_bam_filter_time"]]
    threads: config["popscle"]["popscle_bam_filter_threads"]
    params:
        out_dir = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        tag_group = config["popscle"]["popscle_pileup_tag_group"]
    log: config["outputs"]["output_dir"] + "log/popscle_bam_filter.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools merge -i {input.vcf} | \
            singularity exec --bind {params.bind} {params.sif} samtools view \
                --target-file \
                - \
                --tag-file {params.tag_group}:{input.barcodes} \
                --output {output.bam} \
                --write-index \
                --threads {threads} \
                {input.bam}
        """


rule popscle_pileup:
    input:
        bam = config["outputs"]["output_dir"] + "{pool}/popscle/bam_filter/snpfiltered_alignment.bam",
        sm_list = lambda wildcards: POOL_DF.loc[wildcards.pool, "Individuals"],
        vcf = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz",
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"],
    output:
        pileup = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup.var.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_pileup_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_pileup_time"]]
    threads: config["popscle"]["popscle_pileup_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        tag_group = config["popscle"]["popscle_pileup_tag_group"],
        tag_UMI = config["popscle"]["popscle_pileup_tag_UMI"],
        exclude_flag = config["popscle_extra"]["exclude_flag"],
        out = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup",
        sam_verbose = config["popscle_extra"]["sam_verbose"],
        vcf_verbose = config["popscle_extra"]["vcf_verbose"],
        skip_umi = "--skip-umi" if config["popscle_extra"]["skip_umi"] else "",
        cap_bq = config["popscle_extra"]["cap_bq"],
        min_bq = config["popscle_extra"]["min_bq"],
        min_mq = config["popscle_extra"]["min_mq"],
        min_td = config["popscle_extra"]["min_td"],
        excl_flag = config["popscle_extra"]["excl_flag"],
        min_total = config["popscle_extra"]["min_total"],
        min_uniq = config["popscle_extra"]["min_uniq"],
        min_snp = config["popscle_extra"]["min_snp"],
    log: config["outputs"]["output_dir"] + "log/popscle_pileup.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle dsc-pileup \
            --sam {input.bam} \
            --tag-group {params.tag_group} \
            --tag-UMI {params.tag_UMI} \
            --exclude-flag {params.exclude_flag} \
            --vcf {input.vcf} \
            --sm-list {input.sm_list} \
            --out {params.out} \
            --sam-verbose {params.sam_verbose} \
            --vcf-verbose {params.vcf_verbose} \
            {params.skip_umi} \
            --cap-BQ {params.cap_bq} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --min-TD {params.min_td} \
            --excl-flag {params.excl_flag} \
            --group-list {input.barcodes} \
            --min-total {params.min_total} \
            --min-uniq {params.min_uniq} \
            --min-snp {params.min_snp}
        """


rule popscle_demuxlet:
    input:
        pileup = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup.var.gz",
        vcf = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz",
        group_list = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"],
        sm_list = lambda wildcards: POOL_DF.loc[wildcards.pool, "Individuals"]
    output:
        out = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT.best"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["popscle"]["popscle_demuxlet_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["popscle"]["popscle_demuxlet_time"]]
    threads: config["popscle"]["popscle_demuxlet_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        plp = config["outputs"]["output_dir"] + "{pool}/popscle/pileup/pileup",
        field = config["popscle"]["popscle_demuxlet_field"],
        geno_error_offset = config["popscle_extra"]["geno_error_offset"],
        geno_error_coeff = config["popscle_extra"]["geno_error_coeff"],
        r2_info = config["popscle_extra"]["r2_info"],
        min_mac = config["popscle_extra"]["min_mac"],
        min_callrate = config["popscle_extra"]["min_callrate"],
        out = config["outputs"]["output_dir"] + "{pool}/popscle/demuxlet/demuxletOUT",
        alpha = lambda wildcards: "" if config["popscle_extra"]["alpha"] is None else "--alpha " + str(config["popscle_extra"]["alpha"]),
        doublet_prior = config["popscle_extra"]["doublet_prior"],
        sam_verbose = config["popscle_extra"]["sam_verbose"],
        vcf_verbose = config["popscle_extra"]["vcf_verbose"],
        cap_bq = config["popscle_extra"]["cap_bq"],
        min_bq = config["popscle_extra"]["min_bq"],
        min_mq = config["popscle_extra"]["min_mq"],
        min_td = config["popscle_extra"]["min_td"],
        excl_flag = config["popscle_extra"]["excl_flag"],
        min_total = config["popscle_extra"]["min_total"],
        min_umi = config["popscle_extra"]["min_umi"],
        min_snp = config["popscle_extra"]["min_snp"]
    log: config["outputs"]["output_dir"] + "log/popscle_demuxlet.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} popscle demuxlet \
            --plp {params.plp} \
            --vcf {input.vcf} \
            --field {params.field} \
            --geno-error-offset {params.geno_error_offset} \
            --geno-error-coeff {params.geno_error_coeff} \
            --r2-info {params.r2_info} \
            --min-mac {params.min_mac} \
            --min-callrate {params.min_callrate} \
            --sm-list {input.sm_list} \
            --out {params.out} \
            {params.alpha} \
            --doublet-prior {params.doublet_prior} \
            --sam-verbose {params.sam_verbose} \
            --vcf-verbose {params.vcf_verbose} \
            --cap-BQ {params.cap_bq} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --min-TD {params.min_td} \
            --excl-flag {params.excl_flag} \
            --group-list {input.group_list} \
            --min-total {params.min_total} \
            --min-umi {params.min_umi} \
            --min-snp {params.min_snp}
        """

####################################
############ SOUPORCELL ############
####################################
#!/usr/bin/env python
# Adapted from: https://github.com/wheaton5/souporcell/blob/master/souporcell_pipeline.py
# Author: Martijn Vochteloo
# Note: some functionality from the original souporcell_pipeline is not implemented if it isn't used by us.

# To prevent AmbiguousRuleException between 'souporcell_freebayes' and 'souporcell_freebayes_combine'
wildcard_constraints:
    index="\d+"


rule souporcell_preflights:
    input:
        bam = lambda wildcards: POOL_DF.loc[wildcards.pool, "Bam"],
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"],
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        common_variants = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz"
    output:
        settings = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_settings.json"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_preflights_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_preflights_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_preflights_time"]]
    threads: config["souporcell"]["souporcell_preflights_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/souporcell_preflights.py",
        n_remap_splits = config["souporcell"]["souporcell_remap_splits"],
        n_freebayes_splits = config["souporcell"]["souporcell_freebayes_splits"],
        out_dir = config["outputs"]["output_dir"] + "{pool}/souporcell/",
        clusters = lambda wildcards: POOL_DF.loc[wildcards.pool, "N_Individuals"],
        ploidy = config["souporcell_extra"]["ploidy"],
        min_alt = config["souporcell_extra"]["min_alt"],
        min_ref = config["souporcell_extra"]["min_ref"],
        max_loci = config["souporcell_extra"]["max_loci"],
        restarts = config["souporcell_extra"]["restarts"],
        no_umi = config["souporcell_extra"]["no_umi"],
        umi_tag = config["souporcell_extra"]["umi_tag"],
        cell_tag = "CB"
    log: config["outputs"]["output_dir"] + "log/souporcell_define_remap_bam_regions.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --bam {input.bam} \
            --barcodes {input.barcodes} \
            --fasta {input.fasta} \
            --n_remap_splits {params.n_remap_splits} \
            --n_freebayes_splits {params.n_freebayes_splits} \
            --out_dir {params.out_dir} \
            --clusters {params.clusters} \
            --ploidy {params.ploidy} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --max_loci {params.max_loci} \
            --restarts {params.restarts} \
            --common_variants {input.common_variants} \
            --no_umi {params.no_umi} \
            --umi_tag {params.umi_tag} \
            --cell_tag {params.cell_tag}
        """


rule souporcell_define_remap_bam_regions:
    input:
        settings = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_settings.json",
        bam = lambda wildcards: POOL_DF.loc[wildcards.pool, "Bam"]
    output:
        bam_regions = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/bam_remap_regions.json")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_get_bam_regions_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_get_bam_regions_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_get_bam_regions_time"]]
    threads: config["souporcell"]["souporcell_get_bam_regions_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/get_bam_regions.py",
        n_splits = config["souporcell"]["souporcell_remap_splits"]
    log: config["outputs"]["output_dir"] + "log/souporcell_define_remap_bam_regions.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --bam {input.bam} \
            --n_splits {params.n_splits} \
            --out {output.bam_regions}
        """


def get_remap_bam_region_info(wildcards):
    fh = open(config["outputs"]["output_dir"] + "{pool}/souporcell/bam_remap_regions.json".format(pool=wildcards.pool))
    regions = json.load(fh)
    fh.close()
    combined_regions = []
    for region in regions[wildcards.index]:
        combined_regions.append('"{}:{}:{}"'.format(region["chr"], region["start"], region["stop"]))
    return " ".join(combined_regions)


rule souporcell_make_fastqs:
    input:
        bam = lambda wildcards: POOL_DF.loc[wildcards.pool, "Bam"],
        bam_index = lambda wildcards: POOL_DF.loc[wildcards.pool, "Bam"] + ".bai",
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool,"Barcodes"],
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        fasta_index = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"] + ".fai",
        region = config["outputs"]["output_dir"] + "{pool}/souporcell/bam_remap_regions.json"
    output:
        tmpfq = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/tmp_{index}.fq.gz")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_make_fastqs_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_make_fastqs_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_make_fastqs_time"]]
    threads: config["souporcell"]["souporcell_make_fastqs_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/renamer.py",
        region_array = get_remap_bam_region_info,
        no_umi = config["souporcell_extra"]["no_umi"],
        umi_tag = config["souporcell_extra"]["umi_tag"],
        cell_tag = "CB",
    log: config["outputs"]["output_dir"] + "log/souporcell_make_fastqs.{pool}.{index}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo -n "" > {output.tmpfq}
        for region in {params.region_array}
        do
            IFS=':' read -ra region_array <<< "$region"
            singularity exec --bind {params.bind} {params.sif} echo ${{region_array[0]}} ${{region_array[1]}} ${{region_array[2]}}
            singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
                --bam {input.bam} \
                --barcodes {input.barcodes} \
                --chrom ${{region_array[0]}} \
                --start ${{region_array[1]}} \
                --end ${{region_array[2]}} \
                --no_umi {params.no_umi} \
                --umi_tag {params.umi_tag} \
                --cell_tag {params.cell_tag} | gzip -c >> {output.tmpfq}
        done
        """

# Did not implement HISAT2 method here but it is not part of the image we use anyway.
# Arguments:
# -a = Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
# -x splice = Long-read spliced alignment (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). In the splice mode, 1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 2) long insertions are disabled; 3) deletion and insertion gap costs are different during chaining; 4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes.
# -t = Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences, and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently idle and takes little CPU time).
# -G 50k = Stop chain enlongation if there are no minimizers within NUM-bp [10k].
# -k 21 = Minimizer k-mer length [15]
# -w 11 = Minimizer window size [2/3 of k-mer length]. A minimizer is the smallest k-mer in a window of w consecutive k-mers.
# --sr = Enable short-read alignment heuristics. In the short-read mode, minimap2 applies a second round of chaining with a higher minimizer occurrence threshold if no good chain is found. In addition, minimap2 attempts to patch gaps between seeds with ungapped alignment.
# -A 2 = Matching score [2]
# -B 8 = Mismatching penalty [4]
# -O 12,32 = If query sequence name/length are identical to the target name/length, ignore diagonal anchors. This option also reduces DP-based extension along the diagonal.
# -E 2,1 = Gap extension penalty [2,1]. A gap of length k costs min{O1+k*E1,O2+k*E2}. In the splice mode, the second gap penalties are not used.
# -r 200 = Bandwidth for chaining and base alignment [500,20k]. NUM1 is used for initial chaining and alignment extension; NUM2 for RMQ-based re-chaining and closing gaps in alignments.
# -p .5 = Minimal secondary-to-primary score ratio to output secondary mappings [0.8]. Between two chains overlaping over half of the shorter chain (controlled by -M), the chain with a lower score is secondary to the chain with a higher score. If the ratio of the scores is below FLOAT, the secondary chain will not be outputted or extended with DP alignment later. This option has no effect when -X is applied.
# -N 20 = Output at most INT secondary alignments [5]. This option has no effect when -X is applied.
# -f 1000,5000 = If fraction, ignore top FLOAT fraction of most frequent minimizers [0.0002]. If integer, ignore minimizers occuring more than INT1 times. INT2 is only effective in the --sr or -xsr mode, which sets the threshold for a second round of seeding.
# -n 2 = Discard chains consisting of <INT number of minimizers [3]
# -m 20 = Discard chains with chaining score <INT [40]. Chaining score equals the approximate number of matching bases minus a concave gap penalty. It is computed with dynamic programming.
# -s 40 = Minimal peak DP alignment score to output [40]. The peak score is computed from the final CIGAR. It is the score of the max scoring segment in the alignment and may be different from the total alignment score.
# -g 2000 = Stop chain enlongation if there are no minimizers within NUM-bp [10k].
# -2 = Use two I/O threads during mapping. By default, minimap2 uses one I/O thread. When I/O is slow (e.g. piping to gzip, or reading from a slow pipe), the I/O thread may become the bottleneck. Apply this option to use one thread for input and another thread for output, at the cost of increased peak RAM.
# -K 50m = Number of bases loaded into memory to process in a mini-batch [500M]. Similar to option -I, K/M/G/k/m/g suffix is accepted. A large NUM helps load balancing in the multi-threading mode, at the cost of increased memory.
# --secondary=no = Whether to output secondary alignments [yes]
# -o = Output alignments to FILE [stdout].
# Note: I use minimap2 v2.26, the same version as the authors of the souporcell pipeline, while the old sc-eQTLgen
# used v2.7. The new version gives slightly different results but since this is the version that was intended for souporcell I assume
# it is fine.
rule souporcell_remap:
    input:
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        tmpfq = config["outputs"]["output_dir"] + "{pool}/souporcell/tmp_{index}.fq.gz",
    output:
        samfile = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tmp_{index}.sam")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_remap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_remap_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_remap_time"]]
    threads: config["souporcell"]["souporcell_remap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/souporcell_remap.{pool}.{index}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} minimap2 \
            -ax splice -t {threads} -G50k -k 21 -w 11 --sr -A2 -B8 -O12,32 \
            -E2,1 -r200 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g2000 -2K50m --secondary=no \
            {input.fasta} {input.tmpfq} -o {output.samfile}
        """


# Note: samtools sort gives slightly different results when using old vs new image, syntax is correct
rule souporcell_retag:
    input:
        minimap_tmp_files = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tmp_{index}.sam",
    output:
        retag_bam = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_retag_tmp_{index}.bam"),
        retag_sorted_bam = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_retag_sorted_tmp_{index}.bam")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_retag_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_retag_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_retag_time"]]
    threads: config["souporcell"]["souporcell_retag_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/retag.py",
        no_umi = config["souporcell_extra"]["no_umi"],
        umi_tag = config["souporcell_extra"]["umi_tag"],
        cell_tag = "CB",
    log: config["outputs"]["output_dir"] + "log/souporcell_retag.{pool}.{index}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --sam {input.minimap_tmp_files} \
            --out {output.retag_bam} \
            --no_umi {params.no_umi} \
            --umi_tag {params.umi_tag} \
            --cell_tag {params.cell_tag}

        singularity exec --bind {params.bind} {params.sif} samtools sort {output.retag_bam} -o {output.retag_sorted_bam}
        """


rule souporcell_samtools_merge:
    input:
        retag_sorted_bams = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_retag_sorted_tmp_{index}.bam", pool=wildcards.pool, index=range(config["souporcell"]["souporcell_remap_splits"]))
    output:
        final_bam = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam"),
        final_index = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam.bai")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_souporcell_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_souporcell_time"]]
    threads: config["souporcell"]["souporcell_souporcell_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/souporcell_samtools_merge.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} samtools merge {output.final_bam} {input.retag_sorted_bams}
        singularity exec --bind {params.bind} {params.sif} samtools index {output.final_bam}
        """


rule souporcell_define_bam_regions:
    input:
        bam = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam",
    output:
        bam_regions = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/bam_regions.json")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_get_bam_regions_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_get_bam_regions_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_get_bam_regions_time"]]
    threads: config["souporcell"]["souporcell_get_bam_regions_threads"]
    params:
        sif = config["inputs"]["singularity_image"],
        bind = config["inputs"]["bind_path"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/get_bam_regions.py",
        n_splits = config["souporcell"]["souporcell_freebayes_splits"]
    log: config["outputs"]["output_dir"] + "log/souporcell_define_bam_regions.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --bam {input.bam} \
            --n_splits {params.n_splits} \
            --out {output.bam_regions}
        """


def get_bam_region_info(wildcards):
    fh = open(config["outputs"]["output_dir"] + "{pool}/souporcell/bam_regions.json".format(pool=wildcards.pool))
    regions = json.load(fh)
    fh.close()
    combined_regions = []
    for region in regions[wildcards.index]:
        combined_regions.append("{}:{}-{}".format(region["chr"], region["start"], region["stop"]))
    return " ".join(combined_regions)


# Note, did not implement option where common_variants == None
# Touch the index file since I got some 'The index file is older than the data file:' errors.
# Also, gives slightly different results when using old vs new image, syntax is correct
# IMPORTANT: the old version of samtools depth defaults to a maximum coverage depth of 8000 while the newest version
# has no limit.
rule souporcell_freebayes:
    input:
        bam = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam",
        index = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam.bai",
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"],
        region = config["outputs"]["output_dir"] + "{pool}/souporcell/bam_regions.json"
    output:
        bed = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/depth_{index}.bed"),
        merged_bed = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/depth_{index}_merged.bed")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_freebayes_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_freebayes_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_freebayes_time"]]
    threads: config["souporcell"]["souporcell_freebayes_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        region_args = get_bam_region_info,
        min_cov = config["souporcell_extra"]["min_ref"] + config["souporcell_extra"]["min_alt"],
        max_cov = 100000
    log: config["outputs"]["output_dir"] + "log/souporcell_freebayes.{pool}.{index}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} touch -c {input.index}
        
        singularity exec --bind {params.bind} {params.sif} samtools view -hb {input.bam} {params.region_args} | singularity exec --bind {params.bind} {params.sif} samtools depth - | singularity exec --bind {params.bind} {params.sif} awk '{{ if ($3 >= {params.min_cov} && $3 < {params.max_cov}) {{ print $1 "\t" $2 "\t" $2+1 "\t" $3 }} }}' > {output.bed}
        
        singularity exec --bind {params.bind} {params.sif} bedtools merge -i {output.bed} > {output.merged_bed}
        """


rule souporcell_freebayes_combine:
    input:
        bed = lambda wildcards: expand(config["outputs"]["output_dir"] + "{pool}/souporcell/depth_{index}_merged.bed", pool=POOLS, index=range(config["souporcell"]["souporcell_freebayes_splits"]))
    output:
        bed = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/depth_merged.bed")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_freebayes_combine_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_freebayes_combine_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_freebayes_combine_time"]]
    threads: config["souporcell"]["souporcell_freebayes_combine_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
    log: config["outputs"]["output_dir"] + "log/souporcell_freebayes_combine.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cat {input.bed} > {output.bed}
        """


rule souporcell_common_variants:
    input:
        bed = config["outputs"]["output_dir"] + "{pool}/souporcell/depth_merged.bed",
        common_variants = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz"
    output:
        final_vcf = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/common_variants_covered.vcf.gz")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_common_variants_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_common_variants_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_common_variants_time"]]
    threads: config["souporcell"]["souporcell_common_variants_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
    log: config["outputs"]["output_dir"] + "log/souporcell_common_variants.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} zcat {input.common_variants} | grep "#" | gzip -c > {output.final_vcf}
        singularity exec --bind {params.bind} {params.sif} bedtools intersect \
            -wa \
            -a {input.common_variants} \
            -b {input.bed} | gzip -c >> {output.final_vcf}
        """


# NOTE: cell_tag: DOES NOT WORK, vartrix doesnt support this! https://github.com/wheaton5/souporcell/commit/6872d8803eebd5fd85d16370036aeb2a69942b22
# Skipped parameters:
# --out_variants
# --out_barcodes
# --primary_alignments
# --no_duplicates
rule souporcell_vartrix:
    input:
        final_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/common_variants_covered.vcf.gz",
        final_bam = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam",
        final_index = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_minimap_tagged_sorted.bam.bai",
        barcodes= lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"],
        fasta = config["refs"]["ref_dir"] + config["refs_extra"]["relative_fasta_path"]
    output:
        ref_mtx = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/ref.mtx"),
        alt_mtx = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/alt.mtx")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_vartrix_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_vartrix_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_vartrix_time"]]
    threads: config["souporcell"]["souporcell_vartrix_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        padding = 100,
        scoring_method = "coverage",
        log_level = "error",
        mapq = 30, # Default: 0
        umi = "--umi" if not config["souporcell_extra"]["no_umi"] and config["souporcell_extra"]["umi_tag"] == "UB" else "",
        bam_tag = "CB",
        valid_chars = "ATGCatgc"
    log: config["outputs"]["output_dir"] + "log/souporcell_vartrix.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} vartrix \
            --vcf {input.final_vcf} \
            --bam {input.final_bam} \
            --fasta {input.fasta} \
            --cell-barcodes {input.barcodes} \
            --out-matrix {output.alt_mtx} \
            --padding {params.padding} \
            --scoring-method {params.scoring_method} \
            --ref-matrix {output.ref_mtx} \
            --log-level {params.log_level} \
            --threads {threads} \
            --mapq {params.mapq} \
            {params.umi} \
            --bam-tag {params.bam_tag} \
            --valid-chars {params.valid_chars}
        """


rule souporcell_souporcell:
    input:
        ref_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/ref.mtx",
        alt_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/alt.mtx",
        barcodes = lambda wildcards: POOL_DF.loc[wildcards.pool, "Barcodes"]
    output:
        cluster_file = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/clusters_tmp.tsv.gz"),
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_souporcell_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_souporcell_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_souporcell_time"]]
    threads: config["souporcell"]["souporcell_souporcell_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        clusters = lambda wildcards: POOL_DF.loc[wildcards.pool, "N_Individuals"],
        max_loci = config["souporcell_extra"]["max_loci"],
        min_alt = config["souporcell_extra"]["min_alt"],
        min_ref = config["souporcell_extra"]["min_ref"],
        output_dir = config["outputs"]["output_dir"] + "{pool}/souporcell",
        restarts = config["souporcell_extra"]["restarts"],
        cluster_file = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters_tmp.tsv",
    log: config["outputs"]["output_dir"] + "log/souporcell_souporcell.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} souporcell \
            --alt_matrix {input.alt_mtx} \
            --ref_matrix {input.ref_mtx} \
            --barcodes {input.barcodes} \
            --num_clusters {params.clusters} \
            --min_alt {params.min_alt} \
            --min_ref {params.min_ref} \
            --threads {threads} \
            --restarts {params.restarts} 2> {log} 1> {params.cluster_file}
        singularity exec --bind {params.bind} {params.sif} gzip {params.cluster_file}
        """


# Does not accept gzipped cluster_files.
rule souporcell_doublets:
    input:
        ref_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/ref.mtx",
        alt_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/alt.mtx",
        cluster_file = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters_tmp.tsv.gz"
    output:
        cluster_file = temp(config["outputs"]["output_dir"] + "{pool}/souporcell/clusters_tmp.tsv"),
        doublet_file = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_doublets_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_doublets_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_doublets_time"]]
    threads: config["souporcell"]["souporcell_doublets_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        doublet_prior = 0.5,
        doublet_threshold = 0.9,
        singlet_threshold = 0.9,
        doublet_file = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv",
    log: config["outputs"]["output_dir"] + "log/souporcell_doublets.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} gunzip -c {input.cluster_file} > {output.cluster_file}
        singularity exec --bind {params.bind} {params.sif} troublet \
            --alts {input.alt_mtx} \
            --refs {input.ref_mtx} \
            --clusters {output.cluster_file} \
            --doublet_prior {params.doublet_prior} \
            --doublet_threshold {params.doublet_threshold} \
            --singlet_threshold {params.singlet_threshold} 2> {log} 1> {params.doublet_file}
        singularity exec --bind {params.bind} {params.sif} gzip {params.doublet_file}
        """


rule souporcell_consensus:
    input:
        ref_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/ref.mtx",
        alt_mtx = config["outputs"]["output_dir"] + "{pool}/souporcell/alt.mtx",
        doublet_file = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz",
        final_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/common_variants_covered.vcf.gz"
    output:
        soup_out = config["outputs"]["output_dir"] + "{pool}/souporcell/ambient_rna.txt",
        vcf_out = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf.gz",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_consensus_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_consensus_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_consensus_time"]]
    threads: config["souporcell"]["souporcell_consensus_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        python = "/opt/conda/envs/souporcell/bin/python3.6",
        script = "/opt/souporcell/consensus.py",
        ploidy = config["souporcell_extra"]["ploidy"],
        output_dir = config["outputs"]["output_dir"] + "{pool}/souporcell",
    log: config["outputs"]["output_dir"] + "log/souporcell_consensus.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} {params.python} {params.script} \
            --clusters {input.doublet_file} \
            --alt_matrix {input.alt_mtx} \
            --ref_matrix {input.ref_mtx} \
            --ploidy {params.ploidy} \
            --soup_out {output.soup_out} \
            --vcf_out {output.vcf_out} \
            --output_dir {params.output_dir} \
            --vcf {input.final_vcf}
        """


rule souporcell_summary:
    input:
        clusters = config["outputs"]["output_dir"] + "{pool}/souporcell/clusters.tsv.gz"
    output:
        summary = config["outputs"]["output_dir"] + "{pool}/souporcell/souporcell_summary.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_summary_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_summary_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_summary_time"]]
    threads: config["souporcell"]["souporcell_summary_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"]
    log: config["outputs"]["output_dir"] + "log/souporcell_summary.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $3}}' <(gzip -dc {input.clusters}) \
            | sed -E 's|[0-9]+/[0-9]+|doublet|g' \
            | tail -n+2 \
            | sort \
            | uniq -c \
            | sed -E 's/^ +//g' \
            | sed 's/ /\t/g' \
            | sed '1 i\Assignment N\tClassification' \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print($2,$1)}}' > {output.summary}
        """


rule souporcell_pool_vcf:
    input:
        vcf = config["outputs"]["output_dir"] + "genotypes/vcf_4_demultiplex/imputed_hg38_qc_filtered_exons_sorted.vcf.gz",
        cluster_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf.gz"
    output:
        filtered_refs = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_pool_vcf_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_pool_vcf_time"]]
    threads: config["souporcell"]["souporcell_pool_vcf_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        min_overlap_of_a = 1.0,
        individuals = lambda wildcards: POOL_DF.loc[wildcards.pool, "Individuals"]
    log: config["outputs"]["output_dir"] + "log/souporcell_pool_vcf.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} bedtools intersect \
            -a {input.vcf} \
            -b {input.cluster_vcf} \
            -f {params.min_overlap_of_a} \
            -r \
            -wa \
            -header | \
                singularity exec --bind {params.bind} {params.sif} bcftools view \
                    -S {params.individuals} \
                    -Oz \
                    -o {output.filtered_refs} - 2> {log}
        """


rule souporcell_correlate_genotypes:
    input:
        reference_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/Individual_genotypes_subset.vcf.gz",
        cluster_vcf = config["outputs"]["output_dir"] + "{pool}/souporcell/cluster_genotypes.vcf.gz",
    output:
        correlation_file = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlations.tsv.gz",
        correlation_img = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/ref_clust_pearson_correlation.png",
        assignments = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations/Genotype_ID_key.txt.gz"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_correlations_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["souporcell"]["souporcell_correlations_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["souporcell"]["souporcell_correlations_time"]]
    threads: config["souporcell"]["souporcell_correlations_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Demultiplexing/scripts/Assign_Indiv_by_Geno.R",
        out = config["outputs"]["output_dir"] + "{pool}/souporcell/genotype_correlations"
    log: config["outputs"]["output_dir"] + "log/souporcell_correlate_genotypes.{pool}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --reference_vcf {input.reference_vcf} \
            --cluster_vcf {input.cluster_vcf} \
            --out {params.out}
        """

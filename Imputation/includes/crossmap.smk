#!/usr/bin/env python

# Input: PLINK binary with --max-alleles 2, --mind and --set-all-var-ids applied
# Output: PLINK binary on genome build 38

########################################################
############ CROSSMAP: hg18, b36, hg19, b37 ############
########################################################

def get_chain_file():
    if config["inputs"]["genome_build"] in ["hg18", "b36"]:
        return config["refs_extra"]["relative_hg18_to_hg38_chain_path"]
    elif config["inputs"]["genome_build"] in ["hg19", "b37"]:
        return config["refs_extra"]["relative_grch37_to_grch38_chain_path"]
    else:
        return ""

# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the PVAR with CrossMap's output.
rule crossmap:
    input:
        pgen = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.pgen",
        pvar = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.pvar",
        psam = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.psam",
    output:
        inbed = config["outputs"]["output_dir"] + "crossmapped/input.bed",
        outbed = config["outputs"]["output_dir"] + "crossmapped/output.bed",
        outbed_unmap = config["outputs"]["output_dir"] + "crossmapped/output.bed.unmap",
        excluded_ids = config["outputs"]["output_dir"] + "crossmapped/excluded_ids.txt",
        pgen = config["outputs"]["output_dir"] + "crossmapped/data.pgen",
        pvar = config["outputs"]["output_dir"] + "crossmapped/data.pvar",
        psam = config["outputs"]["output_dir"] + "crossmapped/data.psam",
        original_pvar = config["outputs"]["output_dir"] + "crossmapped/data_original.pvar",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["crossmap"]["crossmap_time"]]
    threads: config["crossmap"]["crossmap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        chain_file = config["refs"]["ref_dir"] + get_chain_file(),
        max_allele_len = config["pre_processing_extra"]["max_allele_len"],
        out = config["outputs"]["output_dir"] + "crossmapped/data",
    log: config["outputs"]["output_dir"] + "log/crossmap.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        singularity exec --bind {params.bind} {params.sif} CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        singularity exec --bind {params.bind} {params.sif} awk '{{print $4}}' {output.outbed_unmap} > {output.excluded_ids}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --exclude {output.excluded_ids} \
            --make-pgen \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {output.pvar} {output.original_pvar}
        singularity exec --bind {params.bind} {params.sif} echo -e "#CHROM\tPOS\tID\tREF\tALT" > {output.pvar}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$1":"$2":"$5"_"$6,$5,$6}}' {output.outbed} >> {output.pvar}
        """

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
# Then replaces the BIM with CrossMap's output.
# This also includes the sort_bed rule.
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
        bed = temp(config["outputs"]["output_dir"] + "crossmapped/data.bed"),
        bim = temp(config["outputs"]["output_dir"] + "crossmapped/data.bim"),
        fam = temp(config["outputs"]["output_dir"] + "crossmapped/data.fam"),
        pgen = config["outputs"]["output_dir"] + "crossmapped/data.pgen",
        pvar = config["outputs"]["output_dir"] + "crossmapped/data.pvar",
        psam = config["outputs"]["output_dir"] + "crossmapped/data.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["crossmap"]["crossmap_time"]]
    threads: config["crossmap"]["crossmap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "crossmapped/data",
        chain_file = config["refs"]["ref_dir"] + get_chain_file()
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
            --make-bed \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --bed {output.bed} \
            --bim {output.bim} \
            --fam {output.fam} \
            --make-pgen \
            --max-alleles 2 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {input.psam} {output.psam}
        """

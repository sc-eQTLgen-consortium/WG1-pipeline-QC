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
# This also includes the rule sort_bed.
# Note: some scary code here: replacing pvar can desynchronize the binary genotype data and the .pvar/.psam indexes if used improperly.
rule crossmap:
    input:
        pgen = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.pgen",
        pvar = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.pvar",
        psam = config["outputs"]["output_dir"] + ("wgs_filtered/" if config["settings"]["is_wgs"] else "pre_processed/") + "data.psam",
    output:
        inbed = config["outputs"]["output_dir"] + "crossmapped/input.bed",
        outbed = config["outputs"]["output_dir"] + "crossmapped/output.bed",
        outbed_unmap = config["outputs"]["output_dir"] + "crossmapped/output.bed.unmap",
        crossmapped_pvar = config["outputs"]["output_dir"] + "crossmapped/crossmapped.pvar",
        excluded_ids = config["outputs"]["output_dir"] + "crossmapped/excluded_ids.txt",
        pgen = config["outputs"]["output_dir"] + "crossmapped/data.pgen",
        pvar = config["outputs"]["output_dir"] + "crossmapped/data.pvar",
        psam = config["outputs"]["output_dir"] + "crossmapped/data.psam",
        log = config["outputs"]["output_dir"] + "crossmapped/data.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["crossmap"]["crossmap_memory"] * config["crossmap"]["crossmap_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["crossmap"]["crossmap_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["crossmap"]["crossmap_time"]]
    threads: config["crossmap"]["crossmap_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        chain_file = config["refs"]["ref_dir"] + get_chain_file(),
        out_tmp = config["outputs"]["output_dir"] + "crossmapped/data",
        script = config["inputs"]["scripts_dir"] + "crossmap_pvar.py",
        out = config["outputs"]["output_dir"] + "crossmapped/data",
    log: config["outputs"]["output_dir"] + "log/crossmap.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        singularity exec --bind {params.bind} {params.sif} CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --pvar {input.pvar} \
            --bed {output.outbed} \
            --out_pvar {output.crossmapped_pvar} \
            --out_exclude {output.excluded_ids}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {output.crossmapped_pvar} \
            --psam {input.psam} \
            --exclude {output.excluded_ids} \
            --sort-vars \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

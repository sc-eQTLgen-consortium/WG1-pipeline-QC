#!/usr/bin/env python

# Input: PLINK binary on genome build 38
# Output: PLINK binary on genome build 38 with updated sex, ancestry and split per ancestry


def get_input_path():
    if config["inputs"]["genome_build"] in ["hg18", "b36", "hg19", "b37"]:
        return "crossmapped/"
    elif config["settings"]["is_wgs"]:
        return "wgs_filtered/"
    else:
        return "pre_processed/"


# Only using chromosome X to reduce space.
rule check_sex:
    input:
        pgen = config["outputs"]["output_dir"] + get_input_path() + "data.pgen",
        pvar = config["outputs"]["output_dir"] + get_input_path() + "data.pvar",
        psam = config["outputs"]["output_dir"] + get_input_path() + "data.psam"
    output:
        bed = config["outputs"]["output_dir"] + "check_sex/data.bed",
        bim = config["outputs"]["output_dir"] + "check_sex/data.bim",
        fam = config["outputs"]["output_dir"] + "check_sex/data.fam",
        log = config["outputs"]["output_dir"] + "check_sex/data.log",
        hh = config["outputs"]["output_dir"] + "check_sex/check_sex.hh",
        nosex = config["outputs"]["output_dir"] + "check_sex/check_sex.nosex",
        sexcheck = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck",
        sexcheck_log = config["outputs"]["output_dir"] + "check_sex/check_sex.log",
        sex_check_tab = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck.tsv",
        man_sex_select = config["outputs"]["output_dir"] + "manual_selection/sex_update_remove.tsv",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["ancestry_sex_qc"]["check_sex_memory"] * config["ancestry_sex_qc"]["check_sex_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["check_sex_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["check_sex_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["check_sex_time"]]
    threads: config["ancestry_sex_qc"]["check_sex_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out_conversion = config["outputs"]["output_dir"] + "check_sex/data",
        out_secheck = config["outputs"]["output_dir"] + "check_sex/check_sex",
    log: config["outputs"]["output_dir"] + "log/check_sex.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --merge-par \
            --chr X \
            --make-bed \
            --out {params.out_conversion}
        singularity exec --bind {params.bind} {params.sif} plink \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --bed {output.bed} \
            --bim {output.bim} \
            --fam {output.fam} \
            --check-sex \
            --out {params.out_secheck}
        singularity exec --bind {params.bind} {params.sif} touch {output.nosex}
        singularity exec --bind {params.bind} {params.sif} sed 's/^ \+//g; s/ \+/\t/g' {output.sexcheck} > {output.sex_check_tab}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NR==1{{print "#"$1,$2,$3,$4,$5,$6,"UPDATE/REMOVE/KEEP"}}' {output.sex_check_tab} > {output.man_sex_select}
        singularity exec --bind {params.bind} {params.sif} grep "PROBLEM" {output.sex_check_tab} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,""}}' >> {output.man_sex_select}
        """


rule prune_1000g:
    input:
        pvar = config["outputs"]["output_dir"] + get_input_path() + "data.pvar",
    output:
        snps_data = config["outputs"]["output_dir"] + "common_snps/snps_data.tsv",
        snps_1000g = config["outputs"]["output_dir"] + "common_snps/snps_1000g.tsv",
        common_pgen = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pgen",
        common_pvar = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pvar",
        common_psam = config["outputs"]["output_dir"] + "common_snps/subset_1000g.psam",
        common_log = config["outputs"]["output_dir"] + "common_snps/subset_1000g.log",
        prune_in_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_1000g_pruning.prune.in",
        prune_out_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_1000g_pruning.prune.out",
        prune_log = config["outputs"]["output_dir"] + "prune_1000g/subset_1000g_pruning.log",
        pruned_pgen = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pgen",
        pruned_pvar = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
        pruned_psam = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.psam",
        pruned_log = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["ancestry_sex_qc"]["prune_1000g_memory"] * config["ancestry_sex_qc"]["prune_1000g_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["prune_1000g_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["prune_1000g_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["prune_1000g_time"]]
    threads: config["ancestry_sex_qc"]["prune_1000g_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["scripts_dir"] + "overlap_pvars.py",
        pgen_1000g = config["refs"]["ref_dir"] + config["refs_extra"]["relative_1000g_files"] + ".pgen",
        pvar_1000g = config["refs"]["ref_dir"] + config["refs_extra"]["relative_1000g_files"] + ".pvar",
        psam_1000g = config["refs"]["ref_dir"] + config["refs_extra"]["relative_1000g_files"] + ".psam",
        indep_window_size = config["ancestry_sex_qc_extra"]["prune_1000g_indep_window_size"],
        indep_step_size = config["ancestry_sex_qc_extra"]["prune_1000g_indep_step_size"],
        pairwise_r2_threshold = config["ancestry_sex_qc_extra"]["prune_1000g_pairwise_r2_threshold"],
        out_common = config["outputs"]["output_dir"] + "common_snps/subset_1000g",
        out_pruning = config["outputs"]["output_dir"] + "prune_1000g/subset_1000g_pruning",
        out_pruned = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g"
    log: config["outputs"]["output_dir"] + "log/prune_1000g.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --pvar1 {input.pvar} \
            --pvar2 {params.pvar_1000g} \
            --variants1 {output.snps_data} \
            --variants2 {output.snps_1000g}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {params.pgen_1000g} \
            --pvar {params.pvar_1000g} \
            --psam {params.psam_1000g} \
            --extract {output.snps_1000g} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out_common}
        
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {output.common_pgen} \
            --pvar {output.common_pvar} \
            --psam {output.common_psam} \
            --indep-pairwise {params.indep_window_size} {params.indep_step_size} {params.pairwise_r2_threshold} \
            --out {params.out_pruning}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {output.common_pgen} \
            --pvar {output.common_pvar} \
            --psam {output.common_psam} \
            --extract {output.prune_in_1000g} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out_pruned}
        """


# use PCA from plink for PCA and projection
# Note: --ac-founders is required in newer versions of plink2.
rule pca_1000g:
    input:
        pgen_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.psam"
    output:
        frq = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.eigenvec",
        eig = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.eigenval",
        log = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.log"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["ancestry_sex_qc"]["pca_1000g_memory"] * config["ancestry_sex_qc"]["pca_1000g_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["pca_1000g_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["pca_1000g_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["pca_1000g_time"]]
    threads: config["ancestry_sex_qc"]["pca_1000g_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs"
    log: config["outputs"]["output_dir"] + "log/pca_1000g.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --freq counts \
            --ac-founders \
            --pca allele-wts \
            --out {params.out}
        """


rule prune_data:
    input:
        pgen = config["outputs"]["output_dir"] + get_input_path() + "data.pgen",
        pvar = config["outputs"]["output_dir"] + get_input_path() + "data.pvar",
        psam = config["outputs"]["output_dir"] + get_input_path() + "data.psam",
        pvar_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
    output:
        snps_data = config["outputs"]["output_dir"] + "prune_1000g/snps_data.tsv",
        snps_1000g = config["outputs"]["output_dir"] + "prune_1000g/snps_1000g.tsv",
        pgen = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.psam",
        log = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["ancestry_sex_qc"]["prune_1000g_memory"] * config["ancestry_sex_qc"]["prune_1000g_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["prune_1000g_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["prune_1000g_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["prune_1000g_time"]]
    threads: config["ancestry_sex_qc"]["prune_1000g_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["scripts_dir"] + "overlap_pvars.py",
        out = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data"
    log: config["outputs"]["output_dir"] + "log/prune_data.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --pvar1 {input.pvar} \
            --pvar2 {input.pvar_1000g} \
            --variants1 {output.snps_data} \
            --variants2 {output.snps_1000g}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --extract {output.snps_data} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """


# TODO: OMP_NUM_THREADS not used?
# use plink pca results to plot with R
rule pca_project:
    input:
        pgen = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.psam",
        frq = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        pgen_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.psam"
    output:
        projected_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_data_pcs.sscore",
        log = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_data_pcs.log",
        projected_1000g_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        log_1000g = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.log",
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["ancestry_sex_qc"]["pca_project_memory"] * config["ancestry_sex_qc"]["pca_project_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["pca_project_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["pca_project_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["pca_project_time"]]
    threads: config["ancestry_sex_qc"]["pca_project_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_data_pcs",
        out_1000g = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected"
    log: config["outputs"]["output_dir"] + "log/pca_project.log"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --read-freq {input.frq} \
            --score {input.eig_all} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --read-freq {input.frq} \
            --score {input.eig_all} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """


rule check_ancestry:
    input:
        projected_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_data_pcs.sscore",
        projected_1000g_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.sscore"
    output:
        anc_fig = report(config["outputs"]["output_dir"] + "check_ancestry/Ancestry_PCAs.png", category="Ancestry", caption="../report_captions/ancestry_pca.rst"),
        pca_anc_cluster_conversion = config["outputs"]["output_dir"] + "check_ancestry/ancestry_cluster_conversion.tsv",
        pca_anc_check = config["outputs"]["output_dir"] + "manual_selection/ancestry_update_remove.tsv",
        anc_maf_select = config["outputs"]["output_dir"] + "manual_selection/ancestry_mafs.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["check_ancestry_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["check_ancestry_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["check_ancestry_time"]]
    threads: config["ancestry_sex_qc"]["check_ancestry_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["scripts_dir"] + "check_ancestry.R",
        out_info = config["outputs"]["output_dir"] + "check_ancestry/",
        out_manual = config["outputs"]["output_dir"] + "manual_selection/"
    log: config["outputs"]["output_dir"] + "log/check_ancestry.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --scores {input.projected_scores} \
            --onekg_scores {input.projected_1000g_scores} \
            --out_info {params.out_info} \
            --out_manual {params.out_manual}
        """


rule summary_ancestry_sex:
    input:
        psam = config["outputs"]["output_dir"] + get_input_path() + "data.psam",
        sex_check_tab = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck.tsv",
        man_sex_select = config["outputs"]["output_dir"] + "manual_selection/sex_update_remove.tsv",
        man_anc_select = config["outputs"]["output_dir"] + "manual_selection/ancestry_update_remove.tsv"
    output:
        sex_summary = report(config["outputs"]["output_dir"] + "metrics/sex_summary.png", category="Ancestry and Sex Summary", caption="../report_captions/sex_summary.rst"),
        ancestry_summary = report(config["outputs"]["output_dir"] + "metrics/ancestry_summary.png", category="Ancestry and Sex Summary", caption="../report_captions/ancestry_summary.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["summary_ancestry_sex_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["ancestry_sex_qc"]["summary_ancestry_sex_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["ancestry_sex_qc"]["summary_ancestry_sex_time"]]
    threads: config["ancestry_sex_qc"]["summary_ancestry_sex_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["scripts_dir"] + "sex_ancestry_summaries.R",
        out = config["outputs"]["output_dir"] + "metrics/",
    log: config["outputs"]["output_dir"] + "log/summary_ancestry_sex.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --psam {input.psam} \
            --sex_check {input.sex_check_tab} \
            --sex_decisions {input.man_sex_select} \
            --ancestry_decisions {input.man_anc_select} \
            --out {params.out}
        """


# Note that we use the updated PSAM which was created in the Snakefile using
# config["outputs"]["output_dir"] + get_input_path() + "data.psam"
# as input.
# Still, some scary code here: replacing psam can desynchronize the binary genotype data and the .pgen/.pvar indexes if used improperly.
rule split_by_ancestry:
    input:
        pgen = config["outputs"]["output_dir"] + get_input_path() + "data.pgen",
        pvar = config["outputs"]["output_dir"] + get_input_path() + "data.pvar"
    output:
        keep = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_individuals.psam",
        bed = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bed",
        bim = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.bim",
        fam = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.fam",
        log = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset.log"
    resources:
        plink_mem_mb = lambda wildcards, attempt: (attempt * config["generic"]["process_plink_memory"] * config["generic"]["process_plink_threads"] - config["settings_extra"]["plink_memory_buffer"]) * 1000,
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["generic"]["process_plink_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["generic"]["process_plink_time"]]
    threads: config["generic"]["process_plink_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        psam = config["outputs"]["output_dir"] + "manual_selection/updated.psam",
        chr_list = ",".join(CHROMOSOMES),
        out = config["outputs"]["output_dir"] + "split_by_ancestry/{ancestry}_subset"
    log: config["outputs"]["output_dir"] + "log/split_by_ancestry.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} grep {wildcards.ancestry} {params.psam} > {output.keep}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --memory {resources.plink_mem_mb} \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {params.psam} \
            --keep {output.keep} \
            --chr {params.chr_list} \
            --make-bed \
            --out {params.out}
        """
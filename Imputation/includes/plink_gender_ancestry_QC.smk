#!/usr/bin/env python


rule indiv_missingness:
    input:
        pgen = plink_gender_ancestry_data + ".pgen",
        pvar = plink_gender_ancestry_data + ".pvar",
        psam = plink_gender_ancestry_data + ".psam",
    output:
        bed = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pgen",
        bim = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pvar",
        fam = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.psam",
        log = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.log",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["indiv_missingness_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["indiv_missingness_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["indiv_missingness_time"]]
    threads: config["plink_gender_ancestry_QC"]["indiv_missingness_threads"]
    params:
       bind = config["inputs"]["bind_path"],
       sif = config["inputs"]["singularity_image"],
       out = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness",
       mind = config["plink_gender_ancestry_QC"]["indiv_missingness_mind"]
    log: config["outputs"]["output_dir"] + "log/indiv_missingness.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --mind {params.mind} \
            --out {params.out}
        """


# TODO: first make-bed just to use plink v1 check-sex?
# TODO: tmp file not needed?
rule check_sex:
    input:
        pgen = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pgen",
        pvar = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pvar",
        psam = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.psam",
    output:
        bed = config["outputs"]["output_dir"] + "check_sex/check_sex.bed",
        bim = config["outputs"]["output_dir"] + "check_sex/check_sex.bim",
        fam = config["outputs"]["output_dir"] + "check_sex/check_sex.fam",
        hh = config["outputs"]["output_dir"] + "check_sex/check_sex.hh",
        nosex = config["outputs"]["output_dir"] + "check_sex/check_sex.nosex",
        sexcheck_tmp = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck",
        sexcheck = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["check_sex_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["check_sex_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["check_sex_time"]]
    threads: config["plink_gender_ancestry_QC"]["check_sex_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "check_sex/check_sex",
    log: config["outputs"]["output_dir"] + "log/check_sex.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --make-bed \
            --max-alleles 2 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink \
            --threads {threads} \
            --bed {output.bed} \
            --bim {output.bim} \
            --fam {output.fam} \
            --check-sex \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.nosex}
        singularity exec --bind {params.bind} {params.sif} sed 's/^ \+//g' {output.sexcheck_tmp} | \
            singularity exec --bind {params.bind} {params.sif} sed 's/ \+/\t/g' > {output.sexcheck}
        """


# TODO: might have duplicates
### Pull just common SNPs between two groups ###
rule common_snps:
    input:
        pgen = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pgen",
        pvar = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pvar",
        psam = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.psam"
    output:
        snps_data = config["outputs"]["output_dir"] + "common_snps/snps_data.tsv",
        snps_1000g = config["outputs"]["output_dir"] + "common_snps/snps_1000g.tsv",
        pgen = config["outputs"]["output_dir"] + "common_snps/subset_data.pgen",
        pvar = config["outputs"]["output_dir"] + "common_snps/subset_data.pvar",
        psam = config["outputs"]["output_dir"] + "common_snps/subset_data.psam",
        log = config["outputs"]["output_dir"] + "common_snps/subset_data.log",
        pgen_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.psam",
        log_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.log",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["common_snps_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["common_snps_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["common_snps_time"]]
    threads: config["plink_gender_ancestry_QC"]["common_snps_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        pgen_1000g = "/opt/1000G/all_phase3_filtered.pgen",
        pvar_1000g = "/opt/1000G/all_phase3_filtered.pvar",
        psam_1000g = "/opt/1000G/all_phase3_filtered.psam",
        out = config["outputs"]["output_dir"] + "common_snps/subset_data",
        out_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g"
    log: config["outputs"]["output_dir"] + "log/common_snps.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.pvar} {params.pvar_1000g} > {output.snps_1000g}
        singularity exec --bind {params.bind} {params.sif} awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {params.pvar_1000g} {input.pvar} > {output.snps_data}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --extract {output.snps_data} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {params.pgen_1000g} \
            --pvar {params.pvar_1000g} \
            --psam {params.psam_1000g} \
            --extract {output.snps_1000g} \
            --make-pgen \
            --out {params.out_1000g}
        """


### Prune with --indep,
rule prune_1000g:
    input:
        pgen_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.psam",
        pgen = config["outputs"]["output_dir"] + "common_snps/subset_data.pgen",
        pvar = config["outputs"]["output_dir"] + "common_snps/subset_data.pvar",
        psam = config["outputs"]["output_dir"] + "common_snps/subset_data.psam",
    output:
        prune_in_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.prune.in",
        prune_out_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.prune.out",
        pgen_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.psam",
        log_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.log",
        data_1000g_key = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data_1000g_key.txt",
        prune_out = config["outputs"]["output_dir"] + "prune_1000g/subset_data.prune.out",
        pgen = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.psam",
        log = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.log",
        pvar_old = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data_original.pvar",
        SNPs2keep = config["outputs"]["output_dir"] + "prune_1000g/SNPs2keep.txt",
        pvar_temp = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data_temp.pvar"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["prune_1000g_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["prune_1000g_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["prune_1000g_time"]]
    threads: config["plink_gender_ancestry_QC"]["prune_1000g_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g",
        out = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data"
    log: config["outputs"]["output_dir"] + "log/prune_1000g.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --indep-pairwise 50 5 0.5 \
            --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --extract {output.prune_out_1000g} \
            --make-pgen \
            --out {params.out_1000g}

        if [[ $(grep "##" {input.pvar} | wc -l) > 0 ]]
        then
            singularity exec --bind {params.bind} {params.sif} grep "##" {input.pvar} > {output.data_1000g_key}
        fi

        singularity exec --bind {params.bind} {params.sif} awk -F"\\t" 'BEGIN{{OFS=FS = "\\t"}} NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.pvar_1000g} {input.pvar} | \
            singularity exec --bind {params.bind} {params.sif} grep -v "##" >> {output.data_1000g_key}
        singularity exec --bind {params.bind} {params.sif} grep -v "##" {output.data_1000g_key} | \
            singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' > {output.prune_out}

        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --extract {output.prune_out} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        
        singularity exec --bind {params.bind} {params.sif} cp {output.pvar} {output.pvar_old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.pvar_old} | \
            singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        singularity exec --bind {params.bind} {params.sif} grep "#CHROM" {output.data_1000g_key} > {output.pvar}
        singularity exec --bind {params.bind} {params.sif} grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.pvar}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.pvar} > {output.pvar_temp}
        singularity exec --bind {params.bind} {params.sif} grep "##" {output.pvar_1000g} > {output.pvar}
        singularity exec --bind {params.bind} {params.sif} cat {output.pvar_temp} >> {output.pvar}
        """
        
        
### put in contingency for duplicated snps - remove from both 1000G and your dataset
rule final_pruning: 
    input:
        pgen = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_data.psam",
    output:
        pgen = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.psam",
        log = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.log",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["final_pruning_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["final_pruning_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["final_pruning_time"]]
    threads: config["plink_gender_ancestry_QC"]["final_pruning_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data"
    log: config["outputs"]["output_dir"] + "log/final_pruning.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --rm-dup 'force-first' \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """


### use PCA from plink for PCA and projection
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
        log = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.log",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_1000g_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_1000g_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["pca_1000g_time"]]
    threads: config["plink_gender_ancestry_QC"]["pca_1000g_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs"
    log: config["outputs"]["output_dir"] + "log/pca_1000g.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


# TODO: is it fine we are using the 1000g acount frq for the input file?
# TODO: OMP_NUM_THREADS not used?
### use plink pca results to plot with R ###
rule pca_project:
    input:
        pgen = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.pgen",
        pvar = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.pvar",
        psam = config["outputs"]["output_dir"] + "final_pruning/final_subset_pruned_data.psam",
        frq = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        pgen_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pgen",
        pvar_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.pvar",
        psam_1000g = config["outputs"]["output_dir"] + "prune_1000g/subset_pruned_1000g.psam"
    output:
        projected_scores = config["outputs"]["output_dir"] + "pca_projection/final_subset_pruned_data_pcs.sscore",
        log = config["outputs"]["output_dir"] + "pca_projection/final_subset_pruned_data_pcs.log",
        projected_1000g_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        log_1000g = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.log"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_project_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_project_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["pca_project_time"]]
    threads: config["plink_gender_ancestry_QC"]["pca_project_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "pca_projection/final_subset_pruned_data_pcs",
        out_1000g = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected"
    log: config["outputs"]["output_dir"] + "log/pca_project.log"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --read-freq {input.frq} \
            --score {input.eig_all} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.pgen_1000g} \
            --pvar {input.pvar_1000g} \
            --psam {input.psam_1000g} \
            --read-freq {input.frq} \
            --score {input.eig_all} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """


# TODO sure this is the right PSAM?
rule pca_projection_assign:
    input:
        projected_scores = config["outputs"]["output_dir"] + "pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = config["outputs"]["output_dir"] + "pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        psam = plink_gender_ancestry_data + ".psam",
        psam_1000g = config["outputs"]["output_dir"] + "common_snps/subset_1000g.psam",
        sexcheck = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck.tsv",
    output:
        anc_fig = report(config["outputs"]["output_dir"] + "pca_sex_checks/Ancestry_PCAs.png", category = "Ancestry", caption = "../report_captions/ancestry_pca.rst"),
        pca_sex_check = config["outputs"]["output_dir"] + "pca_sex_checks/sex_update_remove.tsv",
        pca_anc_check = config["outputs"]["output_dir"] + "pca_sex_checks/ancestry_update_remove.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_projection_assign_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["pca_projection_assign_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["pca_projection_assign_time"]]
    threads: config["plink_gender_ancestry_QC"]["pca_projection_assign_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/PCA_Projection_Plotting.R",
        out = config["outputs"]["output_dir"] + "pca_sex_checks/",
    log: config["outputs"]["output_dir"] + "log/pca_projection_assign.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --scores {input.projected_scores} \
            --onekg_scores {input.projected_1000g_scores} \
            --psam {input.psam} \
            --onekg_psam {input.psam_1000g} \
            --sexcheck {input.sexcheck} \
            --out {params.out}
        """


rule summary_ancestry_sex:
    input:
        sexcheck = config["outputs"]["output_dir"] + "check_sex/check_sex.sexcheck.tsv",
        pca_sex_check = config["outputs"]["output_dir"] + "pca_sex_checks/sex_update_remove.tsv",
        psam = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.psam",
        pca_anc_check = config["outputs"]["output_dir"] + "pca_sex_checks/ancestry_update_remove.tsv"
    output:
        sex_summary = report(config["outputs"]["output_dir"] + "metrics/sex_summary.png", category = "Ancestry and Sex Summary", caption = "../report_captions/sex_summary.rst"),
        ancestry_summary = report(config["outputs"]["output_dir"] + "metrics/ancestry_summary.png", category = "Ancestry and Sex Summary", caption = "../report_captions/ancestry_summary.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["summary_ancestry_sex_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["summary_ancestry_sex_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["summary_ancestry_sex_time"]]
    threads: config["plink_gender_ancestry_QC"]["summary_ancestry_sex_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG1-pipeline-QC/Imputation/scripts/sex_ancestry_summaries.R",
        out = config["outputs"]["output_dir"] + "metrics/",
    log: config["outputs"]["output_dir"] + "log/summary_ancestry_sex.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} \
            --sex_check {input.sexcheck} \
            --sex_decisions {input.pca_sex_check} \
            --psam {input.psam} \
            --ancestry_decisions {input.pca_anc_check} \
            --out {params.out}
        """


rule update_sex_ancestry:
    input:
        indiv_miss_pgen = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pgen",
        indiv_miss_pvar = config["outputs"]["output_dir"] + "indiv_missingness/indiv_missingness.pvar",
        ancestry_updated_psam = config["outputs"]["output_dir"] + "pca_sex_checks/ancestry_updated.psam",
        update_sex = config["outputs"]["output_dir"] + "separate_indivs/sex_update_indivs.tsv",
        remove_indiv = config["outputs"]["output_dir"] + "separate_indivs/remove_indivs.tsv",
    output:
        pgen = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.pgen",
        pvar = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.pvar",
        psam = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.psam",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["update_sex_ancestry_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["update_sex_ancestry_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["update_sex_ancestry_time"]]
    threads: config["plink_gender_ancestry_QC"]["update_sex_ancestry_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated"
    log: config["outputs"]["output_dir"] + "log/update_sex_ancestry.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pgen {input.indiv_miss_pgen} \
            --pvar {input.indiv_miss_pvar} \
            --psam {input.ancestry_updated_psam} \
            --update-sex {input.update_sex} \
            --remove {input.remove_indiv} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

rule subset_ancestry:
    input:
        pgen = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.pgen",
        pvar = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.pvar",
        psam = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated.psam"
    output:
        keep = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_individuals.psam",
        pgen = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.pgen",
        pvar = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.pvar",
        psam = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset.psam"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["subset_ancestry_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["subset_ancestry_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["plink_gender_ancestry_QC"]["subset_ancestry_time"]]
    threads: config["plink_gender_ancestry_QC"]["subset_ancestry_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        infile = config["outputs"]["output_dir"] + "update_sex_ancestry/sex_ancestry_updated",
        out = config["outputs"]["output_dir"] + "subset_ancestry/{ancestry}_subset"
    log: config["outputs"]["output_dir"] + "log/subset_ancestry.{ancestry}.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} grep {wildcards.ancestry} {input.psam} > {output.keep}
        singularity exec --bind {params.bind} {params.sif} plink2 \
            --threads {threads} \
            --pfile {params.infile} \
            --keep {output.keep} \
            --max-alleles 2 \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """



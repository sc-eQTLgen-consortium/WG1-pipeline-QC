#!/usr/bin/env python

if options_dict["ref"] == "hg38":
    rule hg38_liftover:
        input:
            vcf = input_dict["vcf"]
        output:
            output_dict["output_dir"] + "/liftover_hg38_to_hg19/hg19.vcf"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["hg38_liftover_memory"],
            disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["hg38_liftover_memory"]
        threads: plink_gender_ancestry_QC_dict["hg38_liftover_threads"]
        params:
            bind = input_dict["bind_paths"],
            out = output_dict["output_dir"] + "/liftover",
            sif = input_dict["singularity_image"],
            fasta = ref_dict["fasta19"]
        shell:
            """
            singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.mem_per_thread_gb}g -jar /opt/picard/build/libs/picard.jar LiftoverVcf \
                I={input.vcf} \
                O={output} \
                CHAIN=/opt/liftover_refs/hg38ToHg19.over.chain \
                REJECT={params.out}/LiftOver_rejected_variants.vcf \
                R={params.fasta}
            """


elif options_dict["ref"] == "hg19":
    print("Looks like your vcf is already on hg19, no need to liftover for QC and imputation.")
else:
    print("The parameter that you put in for the inputs:ref: in the yaml file is not recognized. It should be either hg19 or hg38.")


if options_dict["ref"] == "hg38":
    input_vcf = output_dict["output_dir"] + "/liftover/hg19.vcf"
elif options_dict["ref"] == "hg19":
    input_vcf = input_dict["vcf"]
else:
    print("There's a problem with the reference name used in the inputs of your yaml file (inputs:ref:). Accepted options are either hg19 or hg38")

rule vcf_to_plink:
    input:
        vcf = input_vcf,
        fam = input_dict["psam_file"]
    output:
        bed = output_dict["output_dir"] + "/plink_hg19/hg19_input.pgen",
        bim = output_dict["output_dir"] + "/plink_hg19/hg19_input.pvar",
        fam = output_dict["output_dir"] + "/plink_hg19/hg19_input.psam",
        indiv_file = output_dict["output_dir"] + "/plink_hg19/individual_file.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["vcf_to_plink_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["vcf_to_plink_memory"]
    threads: plink_gender_ancestry_QC_dict["vcf_to_plink_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/plink_hg19/hg19_input",
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' {input.fam} > {output.indiv_file}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --vcf {input.vcf} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out} --fam {input.fam} --id-delim _ --indiv-sort f {output.indiv_file} --max-alleles 2
        """
        # singularity exec --bind {params.bind} {params.sif} plink --vcf {input.vcf} --make-bed --out {params.out} --fam {input.fam} --indiv-sort f {output.indiv_file}


rule indiv_missingness:
    input:
        bed = output_dict["output_dir"] + "/plink_hg19/hg19_input.pgen",
        bim = output_dict["output_dir"] + "/plink_hg19/hg19_input.pvar",
        fam = output_dict["output_dir"] + "/plink_hg19/hg19_input.psam",
    output:
        bed = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pgen",
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pvar",
        fam = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.psam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["indiv_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["indiv_missingness_memory"]
    threads: plink_gender_ancestry_QC_dict["indiv_missingness_threads"]
    params:
       bind = input_dict["bind_paths"],
       infile = output_dict["output_dir"] + "/plink_hg19/hg19_input",
       out = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness",
       mind = plink_gender_ancestry_QC_dict["indiv_missingness_mind"],
       sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --make-pgen 'psam-cols='fid,parents,sex,phenos --mind {params.mind} --out {params.out}
        """

rule check_sex:
    input:
        bed = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pgen",
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pvar",
        fam = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.psam",
    output:
        bed = output_dict["output_dir"] + "/check_sex/check_sex.bed",
        bim = output_dict["output_dir"] + "/check_sex/check_sex.bim",
        fam = output_dict["output_dir"] + "/check_sex/check_sex.fam",
        hh = output_dict["output_dir"] + "/check_sex/check_sex.hh",
        log = output_dict["output_dir"] + "/check_sex/check_sex.log",
        nosex = output_dict["output_dir"] + "/check_sex/check_sex.nosex",
        sexcheck = output_dict["output_dir"] + "/check_sex/check_sex.sexcheck",
        sexcheck_tab = output_dict["output_dir"] + "/check_sex/check_sex.sexcheck.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["check_sex_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["check_sex_memory"]
    threads: plink_gender_ancestry_QC_dict["check_sex_threads"]
    params:
        bind = input_dict["bind_paths"],
        infile = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness",
        out = output_dict["output_dir"] + "/check_sex/check_sex",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --make-bed --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink --threads {threads} --bfile {params.out} --check-sex --out {params.out}
        singularity exec --bind {params.bind} {params.sif} touch {output.nosex}
        singularity exec --bind {params.bind} {params.sif} sed 's/^ \+//g' {output.sexcheck} | singularity exec --bind {params.bind} {params.sif} sed 's/ \+/\t/g' > {output.sexcheck_tab}
        """



### Pull just common SNPs between two groups ###
rule common_snps:
    input:
        bed = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pgen",
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pvar",
        fam = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.psam",
        log = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.log",
        bim_1000 = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Prepare1000G/plink2_resources_hg19/all_phase3_filtered.pvar"
    output:
        snps_data = output_dict["output_dir"] + "/common_snps/snps_data.tsv",
        snps_1000g = output_dict["output_dir"] + "/common_snps/snps_1000g.tsv",
        bed = output_dict["output_dir"] + "/common_snps/subset_data.pgen",
        bim = output_dict["output_dir"] + "/common_snps/subset_data.pvar",
        fam = output_dict["output_dir"] + "/common_snps/subset_data.psam",
        bed_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.psam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["common_snps_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["common_snps_memory"]
    threads: plink_gender_ancestry_QC_dict["common_snps_threads"]
    params:
        bind = input_dict["bind_paths"],
        infile = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness",
        infile_1000g = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Prepare1000G/plink2_resources_hg19/all_phase3_filtered",
        out = output_dict["output_dir"] + "/common_snps/subset_data",
        out_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g",
        sif = input_dict["singularity_image"]
    shell:
        """
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim} {input.bim_1000} > {output.snps_1000g}
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim_1000} {input.bim} > {output.snps_data}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {output.snps_data} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.snps_1000g} --make-pgen --out {params.out_1000g}
        """

### Prune with --indep,
rule prune_1000g:
    input:
        bed_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.psam",
        bim = output_dict["output_dir"] + "/common_snps/subset_data.pvar",
        bed = output_dict["output_dir"] + "/common_snps/subset_data.pgen",
        fam = output_dict["output_dir"] + "/common_snps/subset_data.psam",
    output:
        # prune_in_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.prune.in",
        prune_out_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.prune.out",
        prune_out = output_dict["output_dir"] + "/common_snps/subset_data.prune.out",
        bed_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.psam",
        bed = output_dict["output_dir"] + "/common_snps/subset_pruned_data.pgen",
        bim = output_dict["output_dir"] + "/common_snps/subset_pruned_data.pvar",
        bim_temp = output_dict["output_dir"] + "/common_snps/subset_pruned_data_temp.pvar",
        bim_old = output_dict["output_dir"] + "/common_snps/subset_pruned_data_original.pvar",
        fam = output_dict["output_dir"] + "/common_snps/subset_pruned_data.psam",
        data_1000g_key = output_dict["output_dir"] + "/common_snps/subset_pruned_data_1000g_key.txt",
        SNPs2keep = output_dict["output_dir"] + "/common_snps/SNPs2keep.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["prune_1000g_memory"]
    threads: plink_gender_ancestry_QC_dict["prune_1000g_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        out_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g",
        infile_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g",
        infile = output_dict["output_dir"] + "/common_snps/subset_data",
        out = output_dict["output_dir"] + "/common_snps/subset_pruned_data"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} \
            --indep-pairwise 50 5 0.5 \
            --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.prune_out_1000g} --make-pgen --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} grep "##" {input.bim} > {output.data_1000g_key}
        singularity exec --bind {params.bind} {params.sif} awk -F"\t" 'BEGIN{{OFS=FS = "\t"}} NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.bim_1000g} {input.bim} | singularity exec --bind {params.bind} {params.sif} grep -v "##" >> {output.data_1000g_key}
        singularity exec --bind {params.bind} {params.sif} grep -v "##" {output.data_1000g_key} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' > {output.prune_out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {output.prune_out} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {output.bim} {output.bim_old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.bim_old} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        singularity exec --bind {params.bind} {params.sif} grep "#CHROM" {output.data_1000g_key} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        singularity exec --bind {params.bind} {params.sif} grep "##" {output.bim_1000g} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} cat {output.bim_temp} >> {output.bim}
        """

### use PCA from plink for PCA and projection
rule pca_1000g:
    input:
        bed_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g.psam" 
    output:
        out = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.acount",
        eig_all = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.eigenvec",
        eig = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.eigenval",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_1000g_memory"]
    threads: plink_gender_ancestry_QC_dict["pca_1000g_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g",
        out = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


### use plink pca results to plot with R ###
rule pca_project:
    input:
        bed = output_dict["output_dir"] + "/common_snps/subset_pruned_data.pgen",
        bim = output_dict["output_dir"] + "/common_snps/subset_pruned_data.pvar",
        fam = output_dict["output_dir"] + "/common_snps/subset_pruned_data.psam",
        frq = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.acount",
        scores = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs.eigenvec.allele"
    output:
        projected_scores = output_dict["output_dir"] + "/common_snps/subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs_projected.sscore"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_project_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["pca_project_memory"]
    threads: plink_gender_ancestry_QC_dict["pca_project_threads"]
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/common_snps/subset_pruned_data",
        infile_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g",
        out = output_dict["output_dir"] + "/common_snps/subset_pruned_data_pcs",
        out_1000g = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs_projected"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """

rule pca_projection_assign:
    input:
        projected_scores = output_dict["output_dir"] + "/common_snps/subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["output_dir"] + "/common_snps/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = output_dict["output_dir"] + "/common_snps/subset_1000g.psam",
        psam = input_dict["psam_file"],
        sexcheck = output_dict["output_dir"] + "/check_sex/check_sex.sexcheck.tsv",
    output:
        sexcheck = output_dict["output_dir"] + "/pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = output_dict["output_dir"] + "/pca_sex_checks/ancestry_update_remove.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    threads: 2
    params:
        variables = output_dict["output_dir"] + "/pca_sex_checks/variables.tsv",
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        outdir = output_dict["output_dir"] + "/pca_sex_checks/",
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-pipeline-QC/Imputation/scripts/PCA_Projection_Plotting.R"
        # script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.outdir} > {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_1000g_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.fam_1000g} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.psam} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.sexcheck} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.variables}
        """

rule separate_indivs:
    input:
        sexcheck = output_dict["output_dir"] + "/pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = output_dict["output_dir"] + "/pca_sex_checks/ancestry_update_remove.tsv"
    output:
        update_sex = output_dict["output_dir"] + "/separate_indivs/sex_update_indivs.tsv",
        remove_indiv = output_dict["output_dir"] + "/separate_indivs/remove_indivs.tsv",
        remove_indiv_temp = output_dict["output_dir"] + "/separate_indivs/remove_indivs_temp.tsv",
        update_ancestry = output_dict["output_dir"] + "/separate_indivs/update_anc_indivs.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
    threads: 1
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} grep "UPDATE" {input.sexcheck} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$4)}}' | singularity exec --bind {params.bind} {params.sif} sed 's/SNPSEX/SEX/g'> {output.update_sex}
        singularity exec --bind {params.bind} {params.sif} grep "REMOVE" {input.sexcheck} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}'> {output.remove_indiv_temp}
        singularity exec --bind {params.bind} {params.sif} grep "REMOVE" {input.anc_check} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' >> {output.remove_indiv_temp}
        singularity exec --bind {params.bind} {params.sif} sort -u {output.remove_indiv_temp} > {output.remove_indiv}
        singularity exec --bind {params.bind} {params.sif} grep "UPDATE" {input.anc_check} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$3,$4,$5,$7)}}' > {output.update_ancestry}
        """


rule update_sex_ancestry:
    input:
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.pgen",
        psam = input_dict["psam_file"],
        update_sex = output_dict["output_dir"] + "/separate_indivs/sex_update_indivs.tsv",
        remove_indiv = output_dict["output_dir"] + "/separate_indivs/remove_indivs.tsv",
        update_ancestry = output_dict["output_dir"] + "/separate_indivs/update_anc_indivs.tsv"
    output:
        bed = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.pgen",
        bim = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.pvar",
        psam = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam",
        psam_temp = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam_temp",
        unique_ancestries = output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["update_sex_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["update_sex_memory"]
    threads: plink_gender_ancestry_QC_dict["update_sex_threads"]
    params:
        infile = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness",
        bind = input_dict["bind_paths"],
        out = output_dict["output_dir"] + "/update_sex_ancestry/update_sex",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --update-sex {input.update_sex} --remove {input.remove_indiv} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} mv {output.psam} {output.psam_temp}
        singularity exec --bind {params.bind} {params.sif} awk 'NR==FNR{{a[$1,$2,$3,$4,$5]=$0;next}} ($1,$2,$3,$4,$5) in a {{$0=a[$1,$2,$3,$4,$5]}}1' {input.update_ancestry} {output.psam_temp} | singularity exec --bind {params.bind} {params.sif} sed 's/PCA_Assignment/Ancestry/g' > {output.psam}
        singularity exec --bind {params.bind} {params.sif} tail -n+2 {output.psam} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($6)}}' | singularity exec --bind {params.bind} {params.sif} sort -u > {output.unique_ancestries}
        """

rule subset_ancestry:
    input:
        unique_ancestries = output_dict["output_dir"] + "/update_sex_ancestry/uniq_acestries.tsv",
        psam = output_dict["output_dir"] + "/update_sex_ancestry/update_sex.psam"
    output:
        keep = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_individuals.psam",
        pgen = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pgen",
        psam = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.psam",
        pvar = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pvar"
    resources: 
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        bind = input_dict["bind_paths"],
        sif = input_dict["singularity_image"],
        infile = output_dict["output_dir"] + "/update_sex_ancestry/update_sex",
        out = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} grep {wildcards.ancestry} {input.psam} > {output.keep}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --keep {output.keep} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """

rule snp_missingness:
    input:
        pgen = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pgen",
        psam = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.psam",
        pvar = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset.pvar"
    output:
        bed = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.pgen",
        bim = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.pvar",
        fam = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.psam",
        log = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.log",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["snp_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["snp_missingness_memory"]
    threads: plink_gender_ancestry_QC_dict["snp_missingness_threads"]
    params:
        bind = input_dict["bind_paths"],
        infile = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset",
        out = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness",
        snp_rate = plink_gender_ancestry_QC_dict["snp_missingness_snp_rate"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.infile} --make-pgen --geno {params.snp_rate} --out {params.out}
        """

rule maf:
    input:
        bed = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.pgen",
        bim = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.pvar",
        fam = output_dict["output_dir"] + "/snp_missingness/{ancestry}_snp_missingness.psam",
    output:
        bed = output_dict["output_dir"] + "/maf/{ancestry}_maf.pgen",
        bim = output_dict["output_dir"] + "/maf/{ancestry}_maf.pvar",
        fam = output_dict["output_dir"] + "/maf/{ancestry}_maf.psam",
        freq = output_dict["output_dir"] + "/maf/{ancestry}_freq.afreq"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["maf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["maf_memory"]
    threads: plink_gender_ancestry_QC_dict["maf_threads"]
    params:
        infile = output_dict["output_dir"] + "/subset_ancestry/{ancestry}_subset",
        bind = input_dict["bind_paths"],
        freq = output_dict["output_dir"] + "/maf/{ancestry}_freq",
        out = output_dict["output_dir"] + "/maf/{ancestry}_maf",
        chr_coding = output_dict["output_dir"] + "/maf/{ancestry}_chr_coding",
        maf = plink_gender_ancestry_QC_dict["maf"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --freq --out {params.freq}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --maf {params.maf} --allow-extra-chr --make-pgen --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {output.bim} {output.bim}.original
        singularity exec --bind {params.bind} {params.sif} sed -i "s/^X/23/g" {output.bim} 
        singularity exec --bind {params.bind} {params.sif} sed -i "s/ID=X/ID=23/g" {output.bim} 
        singularity exec --bind {params.bind} {params.sif} sed -i "s/^Y/24/g" {output.bim} 
        singularity exec --bind {params.bind} {params.sif} sed -i "s/ID=Y/ID=24/g" {output.bim} 
        singularity exec --bind {params.bind} {params.sif} sed -i "s/^MT/26/g" {output.bim} 
        singularity exec --bind {params.bind} {params.sif} sed -i "s/ID=MT/ID=26/g" {output.bim} 
        """

rule hwe:
    input:
        bed = output_dict["output_dir"] + "/maf/{ancestry}_maf.pgen",
        bim = output_dict["output_dir"] + "/maf/{ancestry}_maf.pvar",
        fam = output_dict["output_dir"] + "/maf/{ancestry}_maf.psam",
    output:
        bed = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.pgen",
        bim = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.pvar",
        fam = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.psam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["hwe_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["hwe_memory"]
    threads: plink_gender_ancestry_QC_dict["hwe_threads"]
    params: 
        bind = input_dict["bind_paths"],
        hwe_p = plink_gender_ancestry_QC_dict["hwe_pvalue"],
        maf = output_dict["output_dir"] + "/maf/{ancestry}_maf",
        out = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.maf} --hardy --hwe {params.hwe_p} --make-pgen --out {params.out}
        """

rule het:
    input:
        freq = output_dict["output_dir"] + "/maf/{ancestry}_freq.afreq",
        bed = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.pgen",
        bim = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.pvar",
        fam = output_dict["output_dir"] + "/hwe/{ancestry}_hwe.psam",
        script = input_dict["pipeline_dir"] + "/scripts/filter_het.R"
    output:
        inds = output_dict["output_dir"] + "/het/{ancestry}_het.inds",
        het = output_dict["output_dir"] + "/het/{ancestry}_het.het"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["het_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["het_memory"]
    threads: plink_gender_ancestry_QC_dict["het_threads"]
    params:
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het/{ancestry}_het",
        sif = input_dict["singularity_image"],
        # script = "/opt/WG1-pipeline-QC/Imputation/scripts/filter_het.R"
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-pipeline-QC/Imputation/scripts/filter_het.R"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.hwe} --het --read-freq {input.freq} --out {params.out}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.out}.het {output.inds}
        """


rule het_filter:
    input:
        inds = output_dict["output_dir"] + "/het/{ancestry}_het.inds",
        het = output_dict["output_dir"] + "/het/{ancestry}_het.het"
    output:
        bim = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pgen",
        bed = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pvar",
        fam = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.psam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_gender_ancestry_QC_dict["het_filter_memory"]
    threads: plink_gender_ancestry_QC_dict["het_filter_threads"]
    params:
        bind = input_dict["bind_paths"],
        hwe = output_dict["output_dir"] + "/hwe/{ancestry}_hwe",
        out = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.hwe} --remove {input.inds} --make-pgen --out {params.out}
        """





### Probably don't need this since using grm instead
    # rule ibd:
    #     input:
    #         bim = output_dict["output_dir"] + "/het_filter/het_filter.bim",
    #         bed = output_dict["output_dir"] + "/het_filter/het_filter.bed",
    #         fam = output_dict["output_dir"] + "/het_filter/het_filter.fam"
    #     output:
    #         genome = output_dict["output_dir"] + "/ibd/ibd.genome",
    #         hh = output_dict["output_dir"] + "/ibd/ibd.hh"
    #     resources:
    #         mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
    #         disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
    #     threads: 1
    #     params:
    #         het_filter = output_dict["output_dir"] + "/het_filter/het_filter",
    #         ibd = output_dict["output_dir"] + "/ibd/ibd",
    #         sif = input_dict["singularity_image"]
    #     shell:
    #         """
    #         singularity exec --bind {params.bind} {params.sif} plink --bfile {params.het_filter} --genome --out {params.ibd} --noweb
    #         """

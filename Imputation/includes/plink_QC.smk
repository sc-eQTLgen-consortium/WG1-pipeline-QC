#!/usr/local/envs/py36/bin python3

if options_dict["ref"] == "hg38":
    rule hg38_liftover:
        input:
            vcf = input_dict["vcf"]
        output:
            output_dict["output_dir"] + "/liftover_hg38_to_hg19/hg19.vcf"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["hg38_liftover_memory"],
            disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["hg38_liftover_memory"]
        threads: plink_QC_dict["hg38_liftover_threads"]
        params:
            out = output_dict["output_dir"] + "/liftover",
            sif = input_dict["singularity_image"],
            fasta = ref_dict["fasta19"]
        shell:
            """
            singularity exec {params.sif} java -Xmx{resources.mem_per_thread_gb}g -jar /opt/picard/build/libs/picard.jar LiftoverVcf \
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
        fam = plink_QC_dict["fam_file"]
    output:
        bed = output_dict["output_dir"] + "/plink_hg19/hg19_input.bed",
        bim = output_dict["output_dir"] + "/plink_hg19/hg19_input.bim",
        fam = output_dict["output_dir"] + "/plink_hg19/hg19_input.fam",
        indiv_file = output_dict["output_dir"] + "/plink_hg19/individual_file.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["vcf_to_plink_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["vcf_to_plink_memory"]
    threads: plink_QC_dict["vcf_to_plink_threads"]
    params:
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/plink_hg19/hg19_input",
    shell:
        """
        singularity exec {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' {input.fam} > {output.indiv_file}
        singularity exec {params.sif} plink --vcf {input.vcf} --make-bed --out {params.out} --fam {input.fam}
        """
        # singularity exec {params.sif} plink --vcf {input.vcf} --make-bed --out {params.out} --fam {input.fam} --indiv-sort f {output.indiv_file}

rule snp_missingness:
    input:
        bed = output_dict["output_dir"] + "/plink_hg19/hg19_input.bed",
        bim = output_dict["output_dir"] + "/plink_hg19/hg19_input.bim",
        fam = output_dict["output_dir"] + "/plink_hg19/hg19_input.fam"
    output:
        bed = output_dict["output_dir"] + "/snp_missingness/snp_missingness.bed",
        bim = output_dict["output_dir"] + "/snp_missingness/snp_missingness.bim",
        fam = output_dict["output_dir"] + "/snp_missingness/snp_missingness.fam",
        hh = output_dict["output_dir"] + "/snp_missingness/snp_missingness.hh",
        log = output_dict["output_dir"] + "/snp_missingness/snp_missingness.log",
        nosex = output_dict["output_dir"] + "/snp_missingness/snp_missingness.nosex"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["snp_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["snp_missingness_memory"]
    threads: plink_QC_dict["snp_missingness_threads"]
    params:
        out = output_dict["output_dir"] + "/snp_missingness/snp_missingness",
        snp_rate = plink_QC_dict["snp_missingness_snp_rate"],
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {input} --make-bed --geno {params.snp_rate} --out {params.out} --noweb
        """

rule indiv_missingness:
    input:
        bed = output_dict["output_dir"] + "/snp_missingness/snp_missingness.bed",
        bim = output_dict["output_dir"] + "/snp_missingness/snp_missingness.bim",
        fam = output_dict["output_dir"] + "/snp_missingness/snp_missingness.fam",
        hh = output_dict["output_dir"] + "/snp_missingness/snp_missingness.hh",
        log = output_dict["output_dir"] + "/snp_missingness/snp_missingness.log",
        nosex = output_dict["output_dir"] + "/snp_missingness/snp_missingness.nosex"
    output:
        bed = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.bed",
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.bim",
        fam = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.fam",
        hh = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.hh",
        irem = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.irem",
        log = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.log",
        nosex = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.nosex"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["indiv_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["indiv_missingness_memory"]
    threads: plink_QC_dict["indiv_missingness_threads"]
    params:
       out = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness",
       mind = plink_QC_dict["indiv_missingness_mind"],
       sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {input} --make-bed --mind {params.mind} --out {params.out} --noweb
        """

rule check_sex:
    input:
        bed = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.bed",
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.bim",
        fam = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.fam",
        hh = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.hh",
        irem = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.irem",
        log = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.log",
        nosex = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.nosex"
    output:
        hh = output_dict["output_dir"] + "/check_sex/check_sex.hh",
        log = output_dict["output_dir"] + "/check_sex/check_sex.log",
        nosex = output_dict["output_dir"] + "/check_sex/check_sex.nosex",
        sexcheck = output_dict["output_dir"] + "/check_sex/check_sex.sexcheck"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["check_sex_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["check_sex_memory"]
    threads: plink_QC_dict["check_sex_threads"]
    params:
        out = output_dict["output_dir"] + "/check_sex/check_sex",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {input} --check-sex --out {params.out} --noweb
        """

rule create_sex_check_inds_file:
    input:
    output:
        update_sex = output_dict["output_dir"] + "/check_sex/update_sex.inds",
        wrong_sex = output_dict["output_dir"] + "/check_sex/wrong_sex_remove_individuals.inds"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * 1,
        disk_per_thread_gb=lambda wildcards, attempt: attempt * 1
    threads: 1
    params:
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} printf "FID\tIID\tSEX\n" > {output.update_sex}
        singularity exec {params.sif} printf "FID\tIID\tSEX\n" > {output.wrong_sex}
        """

rule update_sex:
    input:
        bim = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.bim",
        nosex = output_dict["output_dir"] + "/indiv_missingness/indiv_missingness.nosex",
        nosex_sc = output_dict["output_dir"] + "/check_sex/check_sex.nosex",
        sexcheck_sc = output_dict["output_dir"] + "/check_sex/check_sex.sexcheck",
        sexcheck_inds = output_dict["output_dir"] + "/check_sex/check_sex.inds",
        sexcheck_inds_remove = output_dict["output_dir"] + "/check_sex/wrong_sex_remove_individuals.inds"
    output:
        bed = output_dict["output_dir"] + "/update_sex/update_sex.bed",
        bim = output_dict["output_dir"] + "/update_sex/update_sex.bim",
        fam = output_dict["output_dir"] + "/update_sex/update_sex.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["update_sex_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["update_sex_memory"]
    threads: plink_QC_dict["update_sex_threads"]
    params:
        out = output_dict["output_dir"] + "/update_sex/update_sex",
        chr_coding = output_dict["output_dir"] + "/update_sex/chr_coding",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {input} --update-sex {input.sexcheck_inds} --remove {input.sexcheck_inds_remove} --make-bed --out {params.out}
        singularity exec {params.sif} printf "s/X/23/\ns/Y/24/\ns/MT/26/\n" > {params.chr_coding}
        singularity exec {params.sif} paste <(singularity exec {params.sif} sed -f {params.chr_coding} <(singularity exec {params.sif} cut -f1 {params.out}.bim)) <(singularity exec {params.sif} cut -f2- ${params.out}.bim) > {params.out}_tmp.bim
        """

rule maf:
    input:
        bed = output_dict["output_dir"] + "/update_sex/update_sex.bed",
        bim = output_dict["output_dir"] + "/update_sex/update_sex.bim",
        fam = output_dict["output_dir"] + "/update_sex/update_sex.fam"
    output:
        bed = output_dict["output_dir"] + "/maf/maf.bed",
        bim = output_dict["output_dir"] + "/maf/maf.bim",
        fam = output_dict["output_dir"] + "/maf/maf.fam",
        freq = output_dict["output_dir"] + "/maf/freq.frq"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["maf_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["maf_memory"]
    threads: plink_QC_dict["maf_threads"]
    params:
        freq = output_dict["output_dir"] + "/maf/freq",
        out = output_dict["output_dir"] + "/maf/maf",
        chr_coding = output_dict["output_dir"] + "/maf/chr_coding",
        maf = 0.01,
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {input} --freq --out {params.freq} --noweb
        singularity exec {params.sif} plink --bfile {input} --maf {params.maf} --allow-extra-chr --make-bed  --out {params.out}
        singularity exec {params.sif} printf "s/X/23/g" > {params.chr_coding}
        singularity exec {params.sif} printf "s/Y/24/g" >> {params.chr_coding}
        singularity exec {params.sif} printf "s/MT/26/g" >> {params.chr_coding}
        singularity exec {params.sif} paste <(singularity exec {params.sif} sed -f {params.chr_coding} <(singularity exec {params.sif} cut -f1 {output.bim})) <(singularity exec {params.sif} cut -f2- {output.bim}) > {params.out}_tmp.bim
        singularity exec {params.sif} {params.out}_tmp.bim {output.bim}
        """

rule hwe:
    input:
        bed = output_dict["output_dir"] + "/maf/maf.bed",
        bim = output_dict["output_dir"] + "/maf/maf.bim",
        fam = output_dict["output_dir"] + "/maf/maf.fam"
    output:
        bed = output_dict["output_dir"] + "/hwe/hwe.bed",
        bim = output_dict["output_dir"] + "/hwe/hwe.bim",
        fam = output_dict["output_dir"] + "/hwe/hwe.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["hwe_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["hwe_memory"]
    threads: plink_QC_dict["hwe_threads"]
    params:
        maf = output_dict["output_dir"] + "/maf/maf",
        out = output_dict["output_dir"] + "/hwe/hwe",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {params.maf} --hardy --hwe --make-bed --out {params.out} --noweb
        """

rule het:
    input:
        bed = output_dict["output_dir"] + "/hwe/hwe.bed",
        bim = output_dict["output_dir"] + "/hwe/hwe.bim",
        fam = output_dict["output_dir"] + "/hwe/hwe.fam",
        script = input_dict["pipeline_dir"] + "/scripts/filter_het.R"
    output:
        inds = output_dict["output_dir"] + "/het/het.inds",
        het = output_dict["output_dir"] + "/het/het.het"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["het_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["het_memory"]
    threads: plink_QC_dict["het_threads"]
    params:
        hwe = output_dict["output_dir"] + "/hwe/hwe",
        out = output_dict["output_dir"] + "/het/het",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {params.hwe} --het --out {params.out} --noweb
        singularity exec {params.sif} Rscript filter_het.R {params.out}.het {output.inds}
        """


rule het_filter:
    input:
        inds = output_dict["output_dir"] + "/het/het.inds",
        het = output_dict["output_dir"] + "/het/het.het"
    output:
        bim = output_dict["output_dir"] + "/het_filter/het_filter.bim",
        bed = output_dict["output_dir"] + "/het_filter/het_filter.bed",
        fam = output_dict["output_dir"] + "/het_filter/het_filter.fam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * plink_QC_dict["het_filter_memory"]
    threads: plink_QC_dict["het_filter_threads"]
    params:
        hwe = output_dict["output_dir"] + "/hwe/hwe",
        out = output_dict["output_dir"] + "/het_filter/het_filter",
        sif = input_dict["singularity_image"]
    shell:
        """
        singularity exec {params.sif} plink --bfile {params.hwe} --remove {input.inds} --make-bed --out {params.out} --noweb
        """






#### Commands decided not to use - keeping tempsorarily in the case that we might want to use them
    # rule pca:
    #     input:
    #         bim = output_dict["output_dir"] + "/het_filter/het_filter.bim",
    #         bed = output_dict["output_dir"] + "/het_filter/het_filter.bed",
    #         fam = output_dict["output_dir"] + "/het_filter/het_filter.fam"
    #     output:
    #         map_file = output_dict["output_dir"] + "/pca/pca.map",
    #         pca = output_dict["output_dir"] + "/pca/pca.pca",
    #         ped = output_dict["output_dir"] + "/pca/pca.ped"
    #     resources:
    #         mem_per_thread_gb=lambda wildcards, attempt: attempt * 5,
    #         disk_per_thread_gb=lambda wildcards, attempt: attempt * 5
    #     threads: 1
    #     params:
    #         het_filter = output_dict["output_dir"] + "/het_filter/het_filter",
    #         pca = output_dict["output_dir"] + "/pca/pca",
    #         sif = input_dict["singularity_image"]
    #     shell:
    #         """
    #         singularity exec {params.sif} plink --bfile {params.het_filter} --recode --out {params.pca} --noweb
    #         singularity exec {params.sif} printf "genotypename:\t{params.pca}.ped\nsnpname:\t{params.pca}.map\nindivname:\t{params.pca}.ped\nevecoutname:\t{params.pca}.evec\nevaloutname:\t{params.pca}.eval\naltnormstyle:\tNO\nnumoutevec:\t5\nfamilynames:\tNO\nnumoutlieriter:\t0\n" > {params.pca}.smartpca
    #         singularity exec {params.sif} smartpca -p {params.pca}.smartpca | tee {params.pca}.pca
    #         """

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
    #         singularity exec {params.sif} plink --bfile {params.het_filter} --genome --out {params.ibd} --noweb
    #         """

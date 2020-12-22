#!/usr/local/envs/py36/bin python3


rule forward_strand:
    input:
        fasta = ref_dict["fasta19"],
        bim = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pvar",
        bed = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.pgen",
        fam = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter.psam"
    output:
        amb = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_snpflip.ambiguous",
        snpflip = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_snpflip.reverse",
        forward_bed = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand.bed",
        forward_bim = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand.bim",
        forward_fam = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand.fam",
        non_amb_bed = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_noambiguous.bed",
        non_amb_bim = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_noambiguous.bim",
        non_amb_fam = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_noambiguous.fam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["forward_strand_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["forward_strand_memory"]
    threads: HRC_check_dict["forward_strand_threads"]
    params:
        het_filt = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        bind = input_dict["bind_paths"],
        non_amb = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_noambiguous",
        forward_strand = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand",
        snpflip = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand_snpflip",
        sif = input_dict["singularity_image"]
    shell:
        """
        ## Make bed files
        singularity exec --bind {params.bind} {params.sif} plink2 --pfile {params.het_filt} --make-bed --out {params.het_filt}

        # Get strand diagnostics
        singularity exec --bind {params.bind} {params.sif} snpflip -b {params.het_filt}.bim -f {input.fasta} -o {params.snpflip}

        # Remove ambiguous SNPs
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.het_filt} --exclude {output.amb} --make-bed --out {params.non_amb} --noweb

        # Flip SNPs with reverse strand genotypes
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.non_amb} --flip {output.snpflip} --make-bed --out {params.forward_strand} --noweb
        """

rule hrc_check:
    input:
        output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand.bim"
    output:
        freq = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_freq.frq",
        script = output_dict["output_dir"] + "/hrc_check/{ancestry}/Run-plink.sh"
    resources: 
        mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_check_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_check_memory"]
    threads: HRC_check_dict["hrc_check_threads"]
    params:
        prefix = output_dict["output_dir"] + "/het_filter/{ancestry}_het_filter",
        sif = input_dict["singularity_image"],
        freq = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_freq",
        bind = input_dict["bind_paths"],
        out = output_dict["output_dir"] + "/hrc_check/{ancestry}/"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.prefix} --freq --out {params.freq} --noweb
        singularity exec --bind {params.bind} {params.sif} perl /opt/HRC-1000G-check/HRC-1000G-check-bim.pl -b {input} -f {output.freq} -r /opt/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o {params.out}
        """

rule hrc_fix:
    input:
        output_dict["output_dir"] + "/hrc_check/{ancestry}/Run-plink.sh"
    output:
        freq = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_forward_strand-updated-chr23.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_fix_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_fix_memory"]
    threads: HRC_check_dict["hrc_fix_threads"]
    params:
        out = output_dict["output_dir"] + "/hrc_check/{ancestry}/",
        script = output_dict["output_dir"] + "/hrc_check/{ancestry}/Run-plink.sh",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_paths"],
        freq = output_dict["output_dir"] + "/hrc_check/{ancestry}/{ancestry}_freq"
    shell:
        """
        cd {params.out}
        singularity exec --bind {params.bind} {params.sif} sh {params.script}
        """



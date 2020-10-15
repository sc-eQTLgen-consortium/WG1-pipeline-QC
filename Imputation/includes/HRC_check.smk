#!/usr/local/envs/py36/bin python3


rule forward_strand:
    input:
        fasta = ref_dict["fasta19"],
        bim = output_dict["output_dir"] + "/grm_subset/grm_subset.bim"
    output:
        amb = output_dict["output_dir"] + "/forward_strand/forward_strand_snpflip.ambiguous",
        snpflip = output_dict["output_dir"] + "/forward_strand/forward_strand_snpflip.reverse",
        forward_bed = output_dict["output_dir"] + "/forward_strand/forward_strand.bed",
        forward_bim = output_dict["output_dir"] + "/forward_strand/forward_strand.bim",
        forward_fam = output_dict["output_dir"] + "/forward_strand/forward_strand.fam",
        non_amb_bed = output_dict["output_dir"] + "/forward_strand/forward_strand_noambiguous.bed",
        non_amb_bim = output_dict["output_dir"] + "/forward_strand/forward_strand_noambiguous.bim",
        non_amb_fam = output_dict["output_dir"] + "/forward_strand/forward_strand_noambiguous.fam",
        chr_coding = output_dict["output_dir"] + "/forward_strand/forward_strand_chr_coding.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["forward_strand_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["forward_strand_memory"]
    threads: HRC_check_dict["forward_strand_threads"]
    params:
        grm_subset = output_dict["output_dir"] + "/grm_subset/grm_subset",
        bind = input_dict["bind_paths"],
        non_amb = output_dict["output_dir"] + "/forward_strand/forward_strand_noambiguous",
        forward_strand = output_dict["output_dir"] + "/forward_strand/forward_strand",
        snpflip = output_dict["output_dir"] + "/forward_strand/forward_strand_snpflip",
        sif = input_dict["singularity_image"]
    shell:
        """
        # Get strand diagnostics
        singularity exec --bind {params.bind} {params.sif} snpflip -b {input.bim} -f {input.fasta} -o {params.snpflip}

        # Remove ambiguous SNPs
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.grm_subset} --exclude {output.amb} --make-bed --out {params.non_amb} --noweb

        # Flip SNPs with reverse strand genotypes
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.non_amb} --flip {output.snpflip} --make-bed --out {params.forward_strand} --noweb

        ### Update X and MT chromosome names
        singularity exec --bind {params.bind} {params.sif} echo "s/^23/X/" > {output.chr_coding}
        singularity exec --bind {params.bind} {params.sif} echo 's/ID=23/ID=X/g' >> {output.chr_coding}
        singularity exec --bind {params.bind} {params.sif} echo "s/^26/MT/" >> {output.chr_coding}
        singularity exec --bind {params.bind} {params.sif} echo 's/ID=26/ID=MT/g' >> {output.chr_coding}
        singularity exec --bind {params.bind} {params.sif} paste <(singularity exec --bind {params.bind} {params.sif} sed -f {output.chr_coding} <(singularity exec --bind {params.bind} {params.sif} cut -f1 {output.forward_bim})) <(singularity exec --bind {params.bind} {params.sif} cut -f2- {output.forward_bim}) > {params.forward_strand}_tmp.bim
        """



# ### *** need to get the grm_loadings for projection - add to pipeline somehow... *** ###
# rule projection:
#     input:
#         ../data/hapmap3/grm/grm_loadings
#     output:

#     resources:
#         mem_per_thread_gb=lambda wildcards, attempt: attempt * 10,
#         disk_per_thread_gb=lambda wildcards, attempt: attempt * 10
#     threads: 1
#     params:
#         forward_strand = output_dict["output_dir"] + "/forward_strand/forward_strand",
#         projection = output_dict["output_dir"] + "/projection/projection",
#         sif = input_dict["singularity_image"]
#     shell:
#         """
#         singularity exec --bind {params.bind} {params.sif} gcta64 --bfile {params.forward_strand} --project-loading ../data/hapmap3/grm/grm_loadings 6 --out {params.projection}
#         """



rule hrc_check:
    input:
        output_dict["output_dir"] + "/grm_subset/grm_subset.bim"
    output:
        freq = output_dict["output_dir"] + "/hrc_check/freq.frq",
        script = output_dict["output_dir"] + "/hrc_check/Run-plink.sh"
    resources: 
        mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_check_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_check_memory"]
    threads: HRC_check_dict["hrc_check_threads"]
    params:
        prefix = output_dict["output_dir"] + "/grm_subset/grm_subset",
        sif = input_dict["singularity_image"],
        freq = output_dict["output_dir"] + "/hrc_check/freq",
        bind = input_dict["bind_paths"],
        out = output_dict["output_dir"] + "/hrc_check"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink --bfile {params.prefix} --freq --out {params.freq} --noweb
        singularity exec --bind {params.bind} {params.sif} perl /opt/HRC-1000G-check/HRC-1000G-check-bim.pl -b {input} -f {output.freq} -r /opt/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -o {params.out}
        """

# rule hrc_fix:
#     input:
#         output_dict["output_dir"] + "/hrc_check/Run-plink.sh"
#     output:
#         freq = output_dict["output_dir"] + "/hrc_check/freq.frq"
#     resources:
#         mem_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_fix_memory"],
#         disk_per_thread_gb=lambda wildcards, attempt: attempt * HRC_check_dict["hrc_fix_memory"]
#     threads: HRC_check_dict["hrc_fix_threads"]
#     params:
#         out = output_dict["output_dir"] + "/hrc_check/",
#         sif = input_dict["singularity_image"],
#         bind = input_dict["bind_paths"],
#         freq = output_dict["output_dir"] + "/hrc_check/freq"
#     shell:
#         """
#         cd {params.out}
#         singularity exec --bind {params.bind} {params.sif} sh {params.out}/Run-plink.sh
#         """


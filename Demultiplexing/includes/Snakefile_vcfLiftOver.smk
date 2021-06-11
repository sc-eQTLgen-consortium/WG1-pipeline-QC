#!/usr/bin/env python
shell.executable('bash')

rule hg19_liftover:
        input:
            vcf = input_dict["snp_genotypes_filepath"]
        output:
            output_dict["output_dir"] + "/liftover_hg19_to_hg38/hg38.vcf"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * liftover_dict["hg19_liftover_memory"],
            disk_per_thread_gb=lambda wildcards, attempt: attempt * liftover_dict["hg19_liftover_memory"]
        threads: liftover_dict["hg19_liftover_threads"]
        params:
            bind = bind_path,
            out = output_dict["output_dir"] + "/liftover",
            sif = input_dict["singularity_image"],
            fasta = ref_dict["fasta_filepath"]
        shell:
            """
            singularity exec --bind {params.bind} {params.sif} java -Xmx{resources.mem_per_thread_gb}g -jar /opt/picard/build/libs/picard.jar LiftoverVcf \
                I={input.vcf} \
                O={output} \
                CHAIN=/opt/liftover_refs/hg19ToHg38.over.chain \
                REJECT={params.out}/LiftOver_rejected_variants.vcf \
                R={params.fasta}
            """
#!/bin/bash

export SINGULARITY_BINDPATH="/directflow/SCCGGroupShare/projects/,/data:/mnt"

nohup snakemake --snakefile /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/scripts/ONEK1K/hg38/refSNVs/Snakefile \
    --rerun-incomplete --jobs 200 --use-singularity --restart-times 4 \
    --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/logs -o /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/logs -j y -V" \
    > /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/logs/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &





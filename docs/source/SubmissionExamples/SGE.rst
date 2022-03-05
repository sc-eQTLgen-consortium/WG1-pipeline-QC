
SGE Example
-----------

This is an additional example for an SGE cluster.

.. code-block:: bash

    nohup \
      snakemake \
        --snakefile $IMPUTATION_SNAKEFILE \
        --configfile $IMPUTATION_CONFIG \
        --rerun-incomplete \
        --jobs 20 \
        --use-singularity \
        --restart-times 2 \
        --keep-going \
        --cluster \
            "qsub -S /bin/bash \
            -q short.q \
            -r yes \
            -pe smp {threads} \
            -l tmp_requested={resources.disk_per_thread_gb}G \
            -l mem_requested={resources.mem_per_thread_gb}G \
            -e $LOG \
            -o $LOG \
            -j y \
            -V" \
      > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

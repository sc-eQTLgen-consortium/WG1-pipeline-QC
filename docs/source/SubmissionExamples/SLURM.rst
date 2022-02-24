

SLURM Examples
--------------

.. code-block:: bash

  nohup \
    snakemake \
      --snakefile $IMPUTATION_SNAKEFILE \
      --configfile $IMPUTATION_CONFIG \
      --rerun-incomplete \
      --jobs 48 \
      --use-singularity \
      --restart-times 2 \
      --keep-going \
      --cluster \
         "sbatch \
         --qos debug \
         -N 1 \
         --ntasks 1 \
         --cpus-per-task 48 \
         -o $LOG/%{rule}.out \
         --export ALL" \
       > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


Another SLURM example where file latency causes issues with snakemakes ability to detect if a job is completed (note the ``--latency-wait`` parameter):

.. code-block:: bash

  nohup 
    snakemake \
      --snakefile $IMPUTATION_SNAKEFILE \
      --configfile $IMPUTATION_CONFIG \
      --rerun-incomplete \
      --jobs 1 \
      --use-singularity \
      --restart-times 2 \
      --keep-going \
      --latency-wait 30 \
      --cluster \
          "sbatch \
	  --qos regular \
	  -N {threads} \
	  --mem={resources.mem_per_thread_gb}G \
	  --tmp={resources.disk_per_thread_gb}G \
	  -o $LOG/{rule}.out \
	  --export ALL \
	  --time=05:59:59" \
	> $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

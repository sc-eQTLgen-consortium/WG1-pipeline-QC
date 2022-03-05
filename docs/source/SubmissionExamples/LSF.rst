

LSF Example
------------

Here is an example of a submission for a LSF cluster. Of course, you may have to update it based on how your cluster has been set up.

.. code-block:: bash

  nohup \
  snakemake \
  --snakefile $SCEQTL_PIPELINE_SNAKEFILE \
  --configfile $SCEQTL_PIPELINE_CONFIG \
  --rerun-incomplete \
  --jobs 20 \
  --use-singularity \
  --restart-times 2 \
  --keep-going \
  --cluster \
    "bsub \
    -W 24:00 \
    -x \
    -M 10000 \
    -e $LOG \
    -o $LOG" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

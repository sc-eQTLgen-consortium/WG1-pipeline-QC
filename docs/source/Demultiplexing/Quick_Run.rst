.. _Demultiplexing_Quickstart-docs:

Quick Start
=============

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

If you are quite comfortable with Snakemake_, these instructions should be sufficient to run the pipeline.
If you need further details, feel free to go to the more detailed :ref:`Running the Pipeline section <Demultiplexing_Pipeline-docs>`.

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)


.. note::

  With the implementation of newer versions of the pipeline, it is important to make sure your singularity image alignes with the version of the pipeline documentation that you are currently using.
  To check the version of your singluation image plase run:

  .. code-block:: bash

      singularity inspect WG1-pipeline-QC_wgpipeline.sif 

  which will tell you the image version you are currenlty using and, therefore, the relevant documentation for that image.



Quick Run
----------

You will need to run the snakemake pipeline three times (which enables some manual user input that is required). We don't recommend running this pipeline locally because because it will require many parallel jobs, some of which could use upwards of 250GB of memory. 

#. Run the snakemake pipeline (change the qsub command based on your system and preferences)

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


#. Once those jobs have finished, run the snakemake pipeline again (change the qsub command based on your system and preferences)

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


1. Now comes the manual part of the pipeline. Let's start with the ``DoubletDetection`` output. Take a look at the "convergence_test.pdf" in the ``DoubletDetection`` folder in each Pool directory. Ideally, the doublet number predicted should converge toward a common number with each iteration similar to the figure below:

   .. figure:: https://user-images.githubusercontent.com/44268007/104434976-ccf8fa80-55db-11eb-9f30-00f71e4592d4.png
     :width: 400
    
   - In order to indicate whether the pool passed or failed your manual inspection, go to the "DoubletDetection_manual_selection.tsv" located in ``outdir/manual_selections``. This is a tab separated file that has the pools in the first column and a second column header to indicate whether or not the sample passed or failed the manual inspection. For this example, this is what our tsv looks like:
    
   Before user input:

   +------------+-----------------------------+
   |Pool        |  DoubletDetection_PASS_FAIL |
   +============+=============================+
   |test_dataset|                             |
   +------------+-----------------------------+

   After user input:

   +------------+-----------------------------+
   |Pool        |  DoubletDetection_PASS_FAIL |
   +============+=============================+
   |test_dataset|                        PASS |
   +------------+-----------------------------+

.. admonition:: Note
  :class: hint
  
  If the number of doublets do not converge, you can go to the :ref:`Manual Inspection of DoubletDetection and Scrublet Results<manual_selection-docs>` sections to see how to rerun to obtain convergence

1. Next let's check the scrublet results to see if the thresholding was automatically well chosen. Remember that we ran scrublet for each pool with 4 different percentile variable genes: 80, 85, 90 and 95. Take a look at the "doublet_score_histogram.png" in each of the ``scrublet`` directories in each of the pool directories. You want to see that the threshold that was automatically selected nicely separates a bimodal distribution of simulated doublets like below:

   .. figure:: https://user-images.githubusercontent.com/44268007/104436850-016db600-55de-11eb-8f75-229338f7bac7.png

   - In order to identify which scrublet results should be used for downstream analyses, you need to decide which percentile variable gene threshold resulted in the  best simulated doublet bimodal distribution with an effectively set threshold and provide that information in the ``outdir/manual_selections/scrublet_percentile_manual_selection.tsv`` file. For this example, the contents of our ``scrublet_percentile_manual_selection.tsv`` look like this:
      
   +------------+----------------------+
   |Pool        |  scrublet_Percentile |
   +============+======================+
   |test_dataset|                      |
   +------------+----------------------+


   - Enter the percentile variable gene threshold number that resulted in the best bimodal distribution and effectively selected a threshold for the doublet score into the second column of ``scrublet_percentile_manual_selection.tsv``. In our case, the best distribution and threshold selection was for 95th percentile variable genes so we enter the number 95 next to our pool:

   +------------+----------------------+
   |Pool        |  scrublet_Percentile |
   +============+======================+
   |test_dataset| 95                   |
   +------------+----------------------+

   .. admonition:: Note
     :class: hint
      
     If the distribution of the doublet scores do not have two clear peaks, you can go to the :ref:`Manual Inspection of DoubletDetection and Scrublet Results<manual_selection-docs>` sections to see how to rerun to try and get better doublet calling

1. Once you have completed those manual steps, you can run the snakemake pipeline for the final time (change the qsub command based on your system and preferences)

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


A number of QC figures of the singlet droplets have also been produced. 
These can be used to discuss possible QC thresholds with the WG1 and before final QC filtering. 
Let's move to the :ref:`QC Filtering Section <QC_Figures-docs>` to discuss the figures produced and next next steps for additional QC filtering.

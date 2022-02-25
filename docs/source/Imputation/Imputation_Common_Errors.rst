.. _Imputation_Errors-docs:

Common Errors and How to Fix Them
=====================================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues



#. My ``final_pruning`` or my ``sex4imputation`` and ``subset_ancestry`` rules have returned core dumps and segmentation faults. How do I fix this?

   These rules depend on OpenBLAS and need some variables to be set when running from a conda environment. Try running these options before resubmitting the command:

   .. code-block:: bash

    export OMP_NUM_THREADS=1
    export USE_SIMPLE_THREADED_LEVEL3=1


2. The ``harmonize_hg38`` step is failing due to lack of memory
    
   This happens due to system presets for java which is over-riding the memory that is needed for the job. Check what your `$JAVA_OPTS` are set them to a larger Xmx before restarting the pipeline. For example:

   .. code-block:: bash

    JAVA_OPTS="-Xms256m -Xmx25g"




If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

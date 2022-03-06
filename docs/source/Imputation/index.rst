.. _Imputation_Background-docs:

Imputation
=====================================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _minimac4: https://genome.sph.umich.edu/wiki/Minimac4
.. _eagle: https://alkesgroup.broadinstitute.org/Eagle/#x1-320005.3.2
.. _Singularity: https://singularity.lbl.gov/archive/docs/v2-2/index.html
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

 
In preparation for sceQTL detection as part of the sceQTL-Gen consortium, we want to ensure that all the genotypes have been imputed using the same reference, imputation method and processing steps. 
The imputed SNP genotypes will also be needed for demultiplexing and doublet removal. To that end, we have put together instructions for processing and imputing SNP genotype data using eagle_ for phasing and minimac4_ for imputation with the `1000G hg38 reference <https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>`__.

This section of the documentation will provide instructions to run the Imputation pipeline with Snakemake_ from a provided Singularity_ image.
This helps provide consistency across different groups and HPCs.
The first step will be to organize the :ref:`required inputs <Imputation_Input-docs>` for the pipeline.


If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)


.. toctree::
   :maxdepth: 2
   :caption: Imputation
   :hidden:
   
   Imputation_Required_Input
   Imputation_Required_Software
   SNP_Genotype_Imputation
   Imputation_Common_Errors

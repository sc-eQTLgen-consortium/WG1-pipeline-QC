.. _Imputation_Background-docs:

Background and Introduction
=====================================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _minimac4: https://genome.sph.umich.edu/wiki/Minimac4
.. _eagle: https://alkesgroup.broadinstitute.org/Eagle/#x1-320005.3.2

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)


In preparation for sceQTL detection as part of the sceQTL-Gen consortium, we want to ensure that all the genotypes have been imputed using the same reference, imputation method and processing steps. 
The imputed SNP genotypes will also be needed for demultiplexing and doublet removal. To that end, we have put together instructions for processing and imputing SNP genotype data using eagle_ for phasing and minimac4_ for imputation with the `1000G hg38 reference <https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>`__.
.. _Demultiplexing_Introduction-docs:

Background and Introduction
=====================================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues



In order to enable consistency for QC of the datasets used for the sceQTL-Gen Consortium, demultiplexing and doublet detection has been built into a snakemake pipeline leveraging Singularity environments to maintain consistency across different compute structures. We identified five softwares to use for the sceQTL-Gen Consortium by testing the combination of 10 demultiplexing and doublet detection softwares with multiple intersectional methods (manuscript in process). The five softwares that will be used for demultiplexing and doublet detection in this consortium include two SNP-based demultiplexing and doublet detection softwares:

  - popscle-demuxlet_

  - souporcell_

and three transcriptome-based doublet detection softwares:

  - DoubletDetection

  - scds

  - scrublet

The complete pipeline was built in Snakemake in order to provide reproducibility across labs and users with a Singularity image to enable consistency in softwares across systems. Most of the softwares (popscle-demuxlet, souporcell and scds) run without need for any user interaction besides providing input files. However, both DoubletDetection and scrublet require users to check that the thresholding used is effective for the data. Therefore, the pipeline has been built to stop and only run certain jobs after users have provided input regarding each pool they are demultiplexing. Below are the instructions for running the demultiplexing and doublet detection software.


If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

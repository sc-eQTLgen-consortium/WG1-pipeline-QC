.. _Demultiplexing_Introduction-docs:

Demultiplexing
=====================================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _popscle-demuxlet: https://github.com/statgen/popscle
.. _souporcell: https://github.com/wheaton5/souporcell
.. _DoubletDetection: https://github.com/JonathanShor/DoubletDetection
.. _scds: https://www.bioconductor.org/packages/release/bioc/html/scds.html#:~:text=In%20single%20cell%20RNA%20sequencing,in%20scRNA%2Dseq%20data%20computationally.&text=%E2%80%9Cscds%3A%20Computational%20Annotation%20of%20Doublets,RNA%20Sequending%20Data.%E2%80%9D%20bioRxiv.
.. _scrublet: https://github.com/swolock/scrublet



In order to enable consistency for QC of the datasets used for the sceQTL-Gen Consortium, demultiplexing and doublet detection has been built into a snakemake pipeline leveraging Singularity environments to maintain consistency across different compute structures. 
We identified five softwares to use for the sceQTL-Gen Consortium by testing the combination of 10 demultiplexing and doublet detection softwares with multiple intersectional methods (manuscript in process). 
The five softwares that will be used for demultiplexing and doublet detection in this consortium include two SNP-based demultiplexing and doublet detection softwares:

- popscle-demuxlet_

- souporcell_

and three transcriptome-based doublet detection softwares:

- DoubletDetection_

- scds_

- scrublet_

The complete pipeline was built in Snakemake in order to provide reproducibility across labs and users with a Singularity image to enable consistency in softwares across systems. 
Most of the softwares (popscle-demuxlet, souporcell and scds) run without need for any user interaction besides providing input files. 
However, both DoubletDetection and scrublet require users to check that the thresholding used is effective for the data. 
Therefore, the pipeline has been built to stop and only run certain jobs after users have provided input regarding each pool they are demultiplexing. 
Let's first :ref:`prepare the data <Demultiplexing_Input-docs>` and :ref:`software <Demultiplexing_Software-docs>` that we will need to run this pipeline.


If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

.. toctree::
   :maxdepth: 1
   :caption: Demultiplexing and Doublet Removal
   :hidden:
   
   Demultiplexing_Required_Input
   Demultiplexing_Required_Software
   Running_Pipeline
   Quick_Run
   Demultiplexing_Common_Errors

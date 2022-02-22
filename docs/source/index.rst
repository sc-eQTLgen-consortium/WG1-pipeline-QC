.. WG1-pipeline-QC documentation master file, created by
   sphinx-quickstart on Sun Feb 20 15:55:22 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation for WG1 of the sceQTL-Gen Consortium!
===================================================================



.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues

.. figure:: https://user-images.githubusercontent.com/44268007/89252548-35b96f80-d659-11ea-97e9-4b4176df5f08.png


The purpose of this repository is to provide references and instructions for preparation of data for the sceQTL-Gen Consortium. Please note that you can run this pipeline in parallel to the [Working Group 2 (Cell Classification)](https://github.com/sc-eQTLgen-consortium/WG2-pipeline-classification) pipeline both of which will be used for [Working Group 3 (eQTL Detection)](https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL).

There are four major steps that this group is addressing with data preprocessing:
1. The first step is to :ref:`impute the SNP genotypes <Imputation_Background-docs>`. This will be used for demultiplexing and for eQTL detection.
2. The second step is to :ref:`demultiplex and identify doublets<>`. This will allow droplets containing single cells to be assigned to an individual and droplets that contain two cells to be removed.
3. The third step is to :ref:`analyze the quality metrics<>` of the data. These data and results should be fully discussed with members of WG1 too choose effective thresholds for each dataset
4. The last step is the :ref:`final data preparation<>` for WG2 to conduct cell classification

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)



.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. toctree::
   :maxdepth: 2
   :caption: Imputation
   :hidden:
   
   Imputation_Background_Introduction
   Imputation_Required_Input
   Imputation_Required_Software
   SNP_Genotype_Imputation
   Imputation_Common_Errors


.. toctree::
   :maxdepth: 2
   :caption: Demultiplexing and Doublet Removal
   :hidden:
   
   Demultiplexing_Background_Introduction
   Demultiplexing_Required_Input
   Demultiplexing_Required_Software
   Running_Pipeline
   Quick_Run
   Demultiplexing_Common_Errors


.. toctree::
   :maxdepth: 2
   :caption: QC Filtering
   :hidden:
   
   QC_Figures
   Filtering
   Final_Object




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

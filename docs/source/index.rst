.. WG1-pipeline-QC documentation master file, created by
   sphinx-quickstart on Sun Feb 20 15:55:22 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation for WG1 of the sceQTL-Gen Consortium!
===================================================================



.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues

.. figure:: https://user-images.githubusercontent.com/44268007/89252548-35b96f80-d659-11ea-97e9-4b4176df5f08.png
  :width: 300

The purpose of this repository is to provide references and instructions for preparation of data for the sceQTL-Gen Consortium. Please note that you can run this pipeline in parallel to the `Working Group 2 (Cell Classification) <https://github.com/sc-eQTLgen-consortium/WG2-pipeline-classification>`__ pipeline both of which will be used for `Working Group 3 (eQTL Detection) <https://github.com/sc-eQTLgen-consortium/WG3-pipeline-QTL>`__.

There are four major steps that this group is addressing with data preprocessing:

#. The first step is to :ref:`impute the SNP genotypes <Imputation_Background-docs>`. This will be used for demultiplexing and for eQTL detection.

#. The second step is to :ref:`demultiplex and identify doublets<Demultiplexing_Introduction-docs>`. This will allow droplets containing single cells to be assigned to an individual and droplets that contain two cells to be removed.

#. The third step is to :ref:`analyze the quality metrics<QC_Filtering-docs>` of the data. These data and results should be fully discussed with members of WG1 too choose effective thresholds for each dataset

#. The last step is the :ref:`final data preparation<QC_Filtering-docs>` to fileter out doublets.


If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)



.. toctree::
   :titlesonly:
   :hidden:
   
   Imputation/index

.. toctree::
   :maxdepth: 2
   :hidden:
   
   Demultiplexing/index


.. toctree::
   :maxdepth: 2
   :hidden:
   
   QC/index



.. toctree::
   :maxdepth: 2
   :hidden:
   
   SubmissionExamples/index


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

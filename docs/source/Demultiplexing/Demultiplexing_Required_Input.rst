.. _Demultiplexing_Input-docs:

========================
Required Input
========================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _popscle-demuxlet: https://github.com/statgen/popscle
.. _souporcell: https://github.com/wheaton5/souporcell

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

This section explains the data and it's structure which will be required for the Demultiplexing pipeline.

*Please note that all data must be aligned to hg38*


Here is a table of the input you will need:

+--------------------------------------------+----------------------------------+
| Input                                      | Example Included with Test Data? |
+============================================+==================================+
| Sample Table tsv                           | Yes                              |
+--------------------------------------------+----------------------------------+
| Cellranger output directory                | Yes                              |
+--------------------------------------------+----------------------------------+
| Reference genotypes vcf                    | Yes                              |
+--------------------------------------------+----------------------------------+
| Per pool files that contain individual IDs | Yes                              |
+--------------------------------------------+----------------------------------+

For more detailed instructions for each required input, please see below

Test Data
=========

*We have provided a test dataset that contains one pool of a 10x run. This dataset will be used for all example steps below. 
The gzipped directory that contains all the required files can be downloaded. This gzipped directory is ~40Gb. 
If you don't want to download such a large amount of data, we have included a significantly down-sized and sub-sampled version of this dataset provided in the Singularity image (in the `Required Software <Demultiplexing_Software-docs>`) section). 
However, please be aware that the smaller dataset will not provide the expected results due to the downsampling but will enable the pipeline to be tested faster. 
If you want the complete test dataset, use these download instructions:*


.. code-block:: bash

  wget https://www.dropbox.com/s/3oujqq98y400rzz/TestData4PipelineFull.tar.gz
  wget https://www.dropbox.com/s/5n7u723okkf5m3l/TestData4PipelineFull.tar.gz.md5



After downloading the tar.gz directory, it is best to make sure the md5sum of the ``TestData4PipelineFull.tar.gz`` file matches the md5sum in the ``TestData4PipelineFull.tar.gz.md5``:

.. code-block:: bash

  md5sum TestData4PipelineFull.tar.gz > downloaded_TestData4PipelineFull.tar.gz.md5
  diff -s TestData4PipelineFull.tar.gz.md5 downloaded_TestData4PipelineFull.tar.gz.md5


which should return:

.. code-block:: bash

  Files TestData4PipelineFull.tar.gz.md5 and downloaded_TestData4PipelineFull.tar.gz.md5 are identical



Here is the structure of the unzipped TestData4PipelineFull directory (downloaded from dropbox):

.. code-block:: bash
    
  TestData4PipelineFull
  ├── donor_list.txt
  ├── individuals_list_dir
  │   └── test_dataset.txt
  ├── samplesheet.txt
  ├── test_dataset
  │   ├── outs
  │   │   └── filtered_gene_bc_matrices
  │   │       └── Homo_sapiens_GRCh38p10
  │   │           ├── barcodes.tsv
  │   │           ├── genes.tsv
  │   │           └── matrix.mtx
  │   ├── possorted_genome_bam.bam
  │   └── possorted_genome_bam.bam.bai
  └── test_dataset.vcf


Required Data
=============

Sample Table
------------

A tsv file that has Pool names in the first column and the number of individuals per pool in the second column

- Tab separated

- It is assumed that the pool names used here will be somewhere in the directory names for each pool

- This file must have a header

- The Sample Table provided in the test dataset is the `TestData4PipelineFull/samplesheet.txt` file:
  
  +---------------+---------------+
  | Pool          | N Individuals |
  +===============+===============+
  | test_dataset  | 14            |
  +---------------+---------------+

Reference Genotypes Vcf
------------------------

- The vcf should be imputed, filtered for minor allele frequency >= 0.05 and filtered for SNPs that overlap exons. Instructions on preparation of this file are on the `SNP Genotype Imputation <Imputation_Background-docs>`__.

- The vcf provided in the test dataset is the `TestData4PipelineFull/test_dataset.vcf` file

.. admonition:: Important
  :class: caution
  
  This file must NOT be gzipped as souporcell cannot handle vcf.gz files

.. admonition:: Important
  :class: caution
  
  ``popscle`` will error if the order of the chromosomes in this vcf do not match those in your bam file or if your bam uses "chr" encoding ("chr1" instead of "1"). Please check for these possible discrepances and fix the order in the vcf if they do not match. Example code for this is available in the third entry of :ref:`Common Errors and How to Fix them <Demultiplexing_Errors-docs>`.
    

Cellranger output directory
---------------------------

The pipeline assumes a cellranger output or a similar directory structure to below that contains these files:
  
- Bam of aligned reads from single cells

- ``matrix.mtx`` (or ``matrix.mtx.gz``)

- ``genes.tsv`` (or ``features.tsv.gz``)

- ``barcodes.tsv`` (or ``barcodes.tsv.gz``)

Assumed structure for finding the bam and counts file directories:
 
.. code-block:: bash

  parent_data_directory
  ├──Pool1
  │  ├──bam_file.bam
  │  ├──filtered_counts_matrix_dir
  │      ├──barcodes.tsv                     # or barcodes.tsv.gz
  │      ├──genes.tsv                        # or features.tsv.gz
  │      └──matrix.mtx                       # or matrix.mtx.gz
  │  └──...
  ├──Pool2
  │  ├──bam_file.bam
  │  ├──filtered_counts_matrix_dir
  │      ├──barcodes.tsv                     # or barcodes.tsv.gz
  │      ├──genes.tsv                        # or features.tsv.gz
  │      └──matrix.mtx                       # or matrix.mtx.gz
  │  └──...
  └──...


We make the following assumptions when finding files:

- The names of the pool directories are the same as those input into the Sample Table **or** the names of the pools in the Sample Table are contained somewhere within the name of the pool directories that contain the bam and matrix files

- There is only one bam file within the Pool directory

- The matrix, barcode and feature files to be used are downstream of a directory that contains the string "filtered" in the name

The test dataset cellranger output directory is ``TestData4PipelineFull/test_dataset``


Individuals Per Pool
--------------------

Directory that contains one file per pool that has individual IDs for that pool

- Directory should contain a file for each pool that has the ID of each individual that matches the ID used in the reference genotypes vcf

- Each individual ID should be separated by a new line

- No header

- Assumed that the file name contains the pool name somewhere within it

In the test dataset, this file is `TestData4PipelineFull/individuals_list_dir/donor_list.txt`:

.. code-block:: bash

  113_113
  349_350
  352_353
  39_39
  40_40
  41_41
  42_42
  43_43
  465_466
  596_597
  597_598
  632_633
  633_634
  660_661


Next Steps
------------

Now that you have the data prepared, we can move on to getting the :ref:`required software <Demultiplexing_Software-docs>` for the demultiplexing pipeline.
.. _Imputation_Input-docs:

Required Input
========================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _plink2: https://www.cog-genomics.org/plink/2.0/formats

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

This table illustrates the input that the pipeline requires to run and whether it is provided or needs to be prepared and provided by the user. 

+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+
| Input                                      | .. centered:: Category    | .. centered:: User-provided               | Developer-provided                       |
+============================================+===========================+===========================================+==========================================+
| Genotyped SNPs plink2 pfiles (on hg19)     | .. centered:: Study Data  | .. centered:: |:heavy_check_mark:|        | .. centered:: |:heavy_multiplication_x:| |
|                                            |                           |                                           | .. centered:: (Example dataset provided) |
+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+
| 1000G hg38 imputation reference            | .. centered:: Reference   | .. centered:: |:heavy_multiplication_x:|  | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+
| Singularity image                          | .. centered:: Software    | .. centered:: |:heavy_multiplication_x:|  | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+
| Snakemake                                  | .. centered:: Software    | .. centered:: |:heavy_multiplication_x:|  | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+
| scipy                                      | .. centered:: Software    | .. centered:: |:heavy_multiplication_x:|  | .. centered:: |:heavy_check_mark:|       |
+--------------------------------------------+---------------------------+-------------------------------------------+------------------------------------------+


.. _Imputation_Reference-docs:

Reference
---------
We will be using the `1000G hg38 reference <https://www.internationalgenome.org/data-portal/data-collection/30x-grch38>`__ to impute the data and have prepared the reference. You can access it by running 

.. code-block:: bash

  wget https://www.dropbox.com/s/l60a2r3e4vo78mn/eQTLGenImpRef.tar.gz
  wget https://www.dropbox.com/s/eci808v0uepqgcz/eQTLGenImpRef.tar.gz.md5

After downloading the reference, it is best to make sure the md5sum of the ``eQTLGenImpRef.tar.gz`` file matches the md5sum in the ``eQTLGenImpRef.tar.gz.md5``:
  
.. code-block:: bash

  md5sum eQTLGenImpRef.tar.gz > downloaded_eQTLGenImpRef.tar.gz.md5
  diff -s eQTLGenImpRef.tar.gz.md5 downloaded_eQTLGenImpRef.tar.gz.md5


which should return:

.. code-block:: bash

  Files eQTLGenImpRef.tar.gz.md5 and downloaded_eQTLGenImpRef.tar.gz.md5 are identical


If you get anything else, the download was probably incomplete and you should try to download the file again. Then, unpack the contents of the file:

.. code-block:: bash

  tar xvzf eQTLGenImpRef.tar.gz


.. admonition:: Note
  :class: hint

  Some HPCs limit the amount of time that a command can run on a head node, causing it to stop/fail part way through so it is best to untar by using a submission script.


Now you should have the references that are needed to impute the SNP genotype data. You will have the following directory structure:

.. code-block:: bash

  hg38
  ├── imputation
  ├── phasing
  │   ├── genetic_map
  │   └── phasing_reference
  ├── ref_genome_QC
  └── ref_panel_QC


.. _Imputation_Data-docs:

Data
-------
*We have provided a test dataset that can be used to test the pipeline and we have built it in to the singularity image (below).
It will be used for the example below and can be used to test the pipeline. You can also download it directly from https://www.dropbox.com/s/uy9828g1r1jt5xy/ImputationTestDataset_plink.tar.gz and check complete download with https://www.dropbox.com/s/q49gppt7uu75wxr/ImputationTestDataset_plink.tar.gz.md5*

For your own dataset, you will need to make sure you have all the following files in the correct formats. 
You can check the test dataset for an example.

.. _plink2_ref-docs:

Plink2 reference SNP genotype pfiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your reference SNP genotype data will need to be supplied in the plink2_ format which includes 3 files: ``data.pgen``, ``data.psam``, ``data.pvar``

.. admonition:: Important
  :class: caution

  Your chromosome encoding in the ``data.pvar`` file must **not** use 'chr'.
  For example, chromosome 1 would be encoded as '1', not 'chr1'.
  The pipeline will check for this before running and will not run if it finds 'chr' chromsome encoding.


.. admonition:: Important
  :class: caution

  The ``data.psam`` file needs to be in a specific format since it will be important for: 

  - Comparing reported sexes with SNP-genotype predicted sexes

  - Comparing reported ancestries with 1000 Genomes-projected ancestry predictions

  - Creating a per-individual meta-data file for use in WG3 (eQTL detection)


The psam must be tab separated with the following headers and contents should look like this (and requires these headings):

+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
| #FID |  IID   |  PAT   |  MAT  |  SEX | Provided_Ancestry | genotyping_platform                  | array_available | wgs_available | wes_available | age | age_range | Study  | smoking_status | hormonal_contraception_use_currently | menopause | pregnancy_status |
+======+========+========+=======+======+===================+======================================+=================+===============+===============+=====+===========+========+================+======================================+===========+==================+
| 113  |   113  |   0    |     0 |   1  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 78  | 70        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|349   |  350   |   0    |    0  |   1  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 81  | 80        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|352   |   353  |   0    |    0  |   2  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 89  | 80        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|39    |   39   |    0   |     0 |   2  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 56  | 50        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|40    |    40  |    0   |    0  |   2  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 53  | 50        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|41    |    41  |    0   |    0  |   1  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 63  | 60        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|42    |   42   |   0    |    0  |   2  |   EUR             | IlluminaInfiniumGlobalScreeningArray | Y               |  N            | N             | 76  | 70        | OneK1K | NA             | NA                                   | NA        | NA               |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+
|...   | ...    | ...    | ...   | ...  | ...               | ...                                  | ...             | ...           | ...           | ... | ...       | ...    | ...            | ...                                  | ...       | ...              |
+------+--------+--------+-------+------+-------------------+--------------------------------------+-----------------+---------------+---------------+-----+-----------+--------+----------------+--------------------------------------+-----------+------------------+


Key for column contents:

- **#FID**: Family ID

- **IID**: Within-family ID

- **PAT**: Within-family ID of father ('0' if father isn't in dataset)

- **MAT**: Within-family ID of mother ('0' if mother isn't in dataset)

- **SEX**: Sex code ('1' = male, '2' = female, '0' = unknown)

- **Provided_Ancestry**: reported ancestry ('AFR' = African, 'AMR' = Ad Mixed American, 'EAS' = East Asian, 'EUR' = European, 'SAS' = South Asian). If you don't know, use 'NA'.

- **genotyping_platform**: array genotyping was done on

- **array_available**: 'Y' or 'N'; whether SNP genotype array is available for this sample

- **wgs_available**: 'Y' or 'N'; whether whole genome sequencing is available

- **wes_available**: 'Y' or 'N'; whether whole exome sequencing is available

- **age**: age in years of integer, NA if unknown.

- **age_range**: age in decades - lower bound, NA if unknown.

- **Study**: name of the study this donor was included in.

- **smoking_status**: Whether the donor smokes or smoked in the past. Options are:'yes': smokes at time of sample collection, 'past': smoked in the past but not at time of sample collection, 'no': never smoked, 'NA': unknown smoking status.

- **hormonal_contraception_use_currently**: whether the donor is currently using hormonal contraception. Options are: 'yes' (currently using hormonal contraception), 'no' (not  currently using hormonal contraception) or 'NA' (unknown status of contraception use). Note that male donors must be coded as 'NA'.

- **menopause**: Donor menopause status at the time of sample collection. Options are 'pre' (have not yet gone through menopause), 'menopause' (currently going through menopause), 'post' (completed menopause) or 'NA' (unknown menopause status or male). *Note:* that male donors must be coded as 'NA'.

- **pregnancy_status**: Donor pregnancy status at the time of sample collection. Options are 'yes' (pregnant at time of sample collection), 'no' (not pregnant at time of sample collection) or 'NA' (unknown pregnancy status or male). *Note:* that male donors must be coded as 'NA'.

- Any additional metadata can be added as additional columns


.. admonition:: Important
  :class: caution

  The ``data.psam`` file will be used to generate a per-individual meta-data file for use in WG3 (eQTL detection) and will be uploaded to a shared own cloud.
  As such, it is important that you carefully consider whether any individual IDs need to be anonymized.



Next Steps
------------

Now that you have the required inputs organized, you can move on to the :ref:`Required Software <Imputation_Software-docs>` for the imputation pipeline.
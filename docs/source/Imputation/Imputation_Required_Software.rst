.. _Imputation_Software-docs:

Required Software
=================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _Singularity: https://singularity.lbl.gov/archive/docs/v2-2/index.html
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)

This table illustrates the input that the pipeline requires to run and whether it is provided or needs to be prepared and provided by the user. 

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



Singularity Image
-----------------

You will also need the Singularity_ image (has all softwares needed for the pipeline installed for consistency across different systems).
Most HPCs have  Singularity_ installed so you should be able to run  Singularity_ images without any issue. 
If you are unsure, try running `singularity --help`. If  Singularity_ is not installed on your HPC, reach out to your system administrators and see if they can install it. 
Otherwise, feel reach out to us and we can find another solution. 


1. Download the singularity image (~6Gb) with the following command:

   .. code-block:: bash

    wget https://www.dropbox.com/s/ib78j88g5my15fs/WG1-pipeline-QC_imputation.sif
    wget https://www.dropbox.com/s/he89x6xg915b4c3/WG1-pipeline-QC_imputation.sif.md5

   Again, please make sure that the image is uncorrupted by comparing the md5sum


   .. code-block:: bash

    md5sum WG1-pipeline-QC_imputation.sif > downloaded_WG1-pipeline-QC_imputation.sif.md5
    diff -s WG1-pipeline-QC_imputation.sif.md5 downloaded_WG1-pipeline-QC_imputation.sif.md5


   .. note::

    With the implementation of newer versions of the pipeline, it is important to make sure your singularity image alignes with the version of the pipeline documentation that you are currently using.
    To check the version of your singluation image plase run:

    .. code-block:: bash

        singularity inspect WG1-pipeline-QC_imputation.sif

    which will tell you the image version you are currenlty using and, therefore, the relevant documentation for that image.


2. Next, please set up the pipeline with the following command; you will need to provide an absolute directory path to allow singularity to find the correct directories - this must be a directory somewhere above where you current working directory or the current working directory itself:


   .. code-block:: bash

    singularity run --bind <absolute_directory_path> --app setup WG1-pipeline-QC_imputation.sif .


   .. admonition:: Note
    :class: hint

    The pipeline expects certain files pulled from the image to be in the same directory as the singularity image so you will have to rerun the setup steps if you move the image


.. _Imputation_Software_test_data-docs:

Optional
^^^^^^^^
  
If you want to test the pipeline on your system you can use the test dataset we have provided in the image, run the following to get the test dataset from the image to your local directory:

.. code-block:: bash

  singularity run --bind <absolute_directory_path> --app test_dataset WG1-pipeline-QC_imputation.sif .
  tar -xzf ImputationTestDataset_plink.tar.gz

.. _install_snakemake-docs:

Install Snakemake
-----------------
You will also need Snakemake_ and scipy to run the pipeline. We highly recommend you use th conda environment that we have prepared with all the requirements. The conda environment yaml (`snakemake.yaml`) should have been copied into your local directory when you set up the singularity image above.

To make a conda environment from this yaml, you can run:

.. code-block:: bash

  conda env create -f snakemake.yaml -n wg1_snakemake

Then, to activate this environment, run:

.. code-block:: bash

  conda activate wg1_snakemake


If you would prefer to install snakemake and scipy yourself, you can follow the instructions for installing `Snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`__ and then install scipy with ``pip install scipy``.



Next Steps
------------

Now that you have all the required files organized and the software installed, you can move on to running the :ref:`SNP genotype imputation pipeline <Imputation-docs>`.
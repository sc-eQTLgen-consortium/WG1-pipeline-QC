.. _Demultiplexing_Software-docs:

Required Software
===========================
.. _Singularity: https://singularity.lbl.gov/archive/docs/v2-2/index.html
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

This section will provide instructions for installing the required software to run the demultiplexing and doublet detecting pipeline.

Here is a table listing all the required software to run this pipeline:

+--------------------+-----------------------------------+
| Software           | Requires Additional Set up Steps? |
+====================+===================================+
| Singularity Image  | Yes                               |
+--------------------+-----------------------------------+
| Snakemake          | No                                |
+--------------------+-----------------------------------+
| scipy              | No                                |
+--------------------+-----------------------------------+



Singularity Image
-----------------

You will need the Singularity_ image (has all softwares needed for the pipeline installed for consistency across different systems). 
Most HPCs have Singularity_ installed so you should be able to run Singularity_ images without any issue. If you are unsure, try running ``singularity --help``. 
If Singularity_ is not installed on your HPC, reach out to your system administrators and see if they can install it. Otherwise, feel free to open an issue and we can find another solution. 
To install the singularity image with the pipeline:

#. Make a directory where you want to house the pipeline 

#. Move to the directory that you just made to house the pipeline

#. Download the singularity bucket from singularity hub.

   .. code-block:: bash

    wget https://www.dropbox.com/s/t8txw3yttj14we2/WG1-pipeline-QC_wgpipeline.sif
    wget https://www.dropbox.com/s/n2vugd6w7glyuft/WG1-pipeline-QC_wgpipeline.sif.md5


   After downloading the image, it is best to make sure the md5sum of the `WG1-pipeline-QC_wgpipeline.sif` file matches the md5sum in the `WG1-pipeline-QC_wgpipeline.sif.md5`:

   .. code-block:: bash

    md5sum WG1-pipeline-QC_wgpipeline.sif > downloaded_WG1-pipeline-QC_wgpipeline.sif.md5
    diff -s WG1-pipeline-QC_wgpipeline.sif.md5 downloaded_WG1-pipeline-QC_wgpipeline.sif.md5


   which should return:

   .. code-block:: bash

      Files WG1-pipeline-QC_wgpipeline.sif.md55 and downloaded_WG1-pipeline-QC_wgpipeline.sif.md5 are identical

   .. note::

      With the implementation of newer versions of the pipeline, it is important to make sure your singularity image alignes with the version of the pipeline documentation that you are currently using.
      To check the version of your singluation image plase run:

      .. code-block:: bash

         singularity inspect WG1-pipeline-QC_wgpipeline.sif 

      which will tell you the image version you are currenlty using and, therefore, the relevant documentation for that image.


#. Set up the pipeline with the following command; you will need to provide an absolute directory path to allow singularity to find the correct directories - this must be a directory somewhere above where you current working directory or the current working directory itself:

   .. code-block:: bash

    singularity run --bind <absolute_directory_path> --app setup WG1-pipeline-QC_wgpipeline.sif .

   - If you want to use the smaller test dataset to test the pipeline and installation, please also run:

   .. code-block:: bash

    singularity run --bind <absolute_directory_path> --app test_dataset WG1-pipeline-QC_wgpipeline.sif .
    tar xzvf TestData4PipelineSmall.tar.gz

   - The directory structure is exactly the same as the full test dataset but the parent directory is called ``TestData4PipelineSmall``

   - This should copy all the files from the singularity bucket that are needed to run the pipeline to the current directory
   
   .. admonition:: Note
      :class: hint
     
      The pipeline assumes that the files are in the same directory as the singularity image so it is important that you do not copy them to a different location - follow the instructions above and you shouldn't have any issues


Snakemake
---------

You will also need Snakemake_ and scipy to run the pipeline. You can either use a conda environment that we have prepared with all the requirements (recommended) or install these yourself. 
You likely already created an environment when preparing for the :ref:`SNP Genotype Imputation<Imputation_Background-docs>` steps. If you didn't, you can find the steps for the Snakemake_ conda environment :ref:`here <install_snakemake-docs>`.



Next Steps
------------

Now you have prepared the required data and installed the required software, we can move on to :ref:`running the pipeline <Demultiplexing_Pipeline-docs>`.
If you are familiar and comfortable with Snakemake_ pipelines and structures, you can use the :ref:`Quick Run instructions <Demultiplexing_Quickstart-docs>`.
However, they contain far less detail than the main sections so we only recommend using the :ref:`Quick Run instructions <Demultiplexing_Quickstart-docs>` if you have quite a bit of experience with Snakemake_ pipelines.
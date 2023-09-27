.. _Demultiplexing_Errors-docs:

Common Errors and How to Fix Them
=====================================

The following is a list of common errors experienced by users of this pipeline and ways to solve those errors.

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _Snakemake: 

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)



#. **I am getting an error that says that the file doesn't exist but when I type ``ll /path/to/file``, it is there.**

   This error will occur if the directory or parent directory to the file has not been exported to the ``$SINGULARITY_BINDPATH``. Please check that you have indeed exported a parent directory to the file to the ``$SINGULARITY_BINDPATH``.

#. **Snakemake is telling me that my directory cannot be locked:**

   .. code-block:: bash

      Building DAG of jobs...
      Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
      /path/to/directory
      If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the ``--unlock`` argument.

   This error will occur for one of two reasons:

   - If you still have jobs running. Snakemake won't allow you to start more than one instance at a time. This is to ensure that it doesn't start the same job more than once since it depends on the file outputs to identify what has been run. 

   - One of the jobs failed without providing an exit code for snakemake. If you check your jobs that are running on the cluster and see that there are no snakemake jobs, you can unlock the directory by running the same command with ``--unlock``. This will unlock the directory and then you can run the snakemake command.

   - On our HPC, this happens often for popscle-pileup and results in a core dump file. If you find the same on your system, you may need to increase the memory allotment in the configuration yaml


#. **I am getting an error from the `popscle_pileup` step telling me that my chromosome ordering is different in my vcf and my bam:**

   .. code-block:: bash

      FATAL ERROR - 
      [E:int32_t cmdCramDscPileup(int32_t, char**)] Your VCF/BCF files and SAM/BAM/CRAM files have different ordering of chromosomes. 
      SAM/BAM/CRAM file has 19 before 2, but VCF/BCF file has 19 after 2


   - This error happens when the order of the chromosomes in the SNP genotype vcf are not in the same order as those in your bam file.

   - To fix this, you will have to first reorder the chromosomes in the vcf header and then reorder the chromosomes in the SNP genotypes by the chromosome order in the header. If your bam uses "chr" encoding for chromosomes ("chr1" instead of "1"), then you will also have to update the vcf to have a "chr" before each chromosome number. Here are some instructions on how to do that:

     #. Update the header with ``bcftools reheader`` using the fai from the reference that was used to align reads for the single cell data. We recommend using bcftools from the singularity image to do this since this functionality is only available in more recentt ``bcftools``:

        .. code-block:: bash

            singularity exec --bind $BIND_DIR $SIF bcftools reheader -f $FAI imputed_hg38_qc_filtered_exons_sorted.vcf > imputed_hg38_qc_filtered_exons_sorted_reheader.vcf


        .. admonition:: Note
           :class: hint
               
           ``$BIND_DIR`` is the path to a parent directory of all the files being used in the command, ``$SIF`` is the path to the singularity image.if not all your files are below the directory you are running from and ``$FAI`` is the path to your fast fai file.

     #. *Optional*: if your bam file uses chr encoding for the chromosomes (*ie* "chr1" instead of "1"), you will also need to add chr onto the beginning of each chromosome in the vcf which you can do with the following code:

        .. code-block:: bash

            awk '{
            if($0 !~ /^#/)
                  print "chr"$0;
            else print $0
            }' imputed_hg38_qc_filtered_exons_sorted_reheader.vcf > imputed_hg38_qc_filtered_exons_sorted_reheader_chr.vcf

     #. Finally, you can update the order of the SNPs in the file to match the order of the chromosomes in the header with picard. Again, you can run this directly from the singularity image: 

        .. code-block:: bash

            singularity exec --bind $BIND_DIR $SIF java -jar /opt/picard/build/libs/picard.jar SortVcf \
            I=imputed_hg38_qc_filtered_exons_sorted_reheader.vcf \
            O=imputed_hg38_qc_filtered_exons_sorted_reheader_reorder.vcf


        .. admonition:: Note
         :class: hint
         Where ``$BIND_DIR`` is the path to a parent directory of all the files being used in the command and ``$SIF`` is the path to the singularity image.if not all your files are below the directory you are running from.


#. **Some of my rules are running fine but DoubletDetection and scrublet are returning core dumps and segmentation faults, I don't know what to do?**
    
   - We have found that there are some environment variables that need to be set for OpenBLAS when running from snakemake built in conda. You can set these variables by running the following code in your terminal before executing your snakemake command:
      
     .. code-block:: bash

         export OMP_NUM_THREADS=1
         export USE_SIMPLE_THREADED_LEVEL3=1


#. **I am receiving an error from snakemake about ``PuLP`` and ``coincbc``:**

   .. code-block:: bash

      WorkflowError:
      You need to install at least one LP solver compatible with PuLP (e.g. coincbc). See https://coin-or.github.io/pulp for details. Alternatively, run Snakemake_ with ``--scheduler greedy``.

   but I have both PuLP and coincbc installed and when I run with ``--scheduler greedy``, it doesn't work.

   - We have noticed this issue with snakmake version 5.25.0. PuLP was introduced in snakemake version 5.23.0 so you can try to install an earlier version of snakemake in order to omit this error. 

   - You can also try using ``drmaa`` by installing it and using ``--drmaa`` in place of your ``--cluster`` option and remove the "qsub" or equivalent portion of the command


If you have any questions or issues while running the pipeline, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au).


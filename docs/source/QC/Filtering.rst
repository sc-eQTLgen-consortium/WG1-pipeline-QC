.. _Filtering:

Filtering Based on QC Metrics
===============================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues
.. _Seurat: https://satijalab.org/seurat/

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)


Once appropriate QC thresholds for your data have been identified after discussions with other members from WG1, the cells that fall outside of those thresholds need to be removed. 
We will use the Seurat_ package functions to perform QC filtering and and produce a final seurat object that can be directly used by WG2 for cell type prediction.


Since there are a few possible methods that thresholds could be selected depending on your dataset, we will provide pseudo-code for multiple possibilities. Please choose the scenario that best fits your data:

Filter All Pools Together
-------------------------

If the characteristics of all the pools in your dataset are pretty similar, it may be deemed appropriate (through discussions with WG1) to filter all pools on a common threshold. 
If this is the case, either an absolute threshold may be chosen or a certain MAD from the median may be chosen for each QC metric. 
We will provide some pseudo-code below for each case. If your pools are significantly different on these QC metrics, it may be desirable to filter each pool with a different threshold and you can jump to the :ref:`Filter Each Pool Separately section <filter_separate-docs>`.


Absolute Threshold
^^^^^^^^^^^^^^^^^^^

Here is some example code that you can use if you would like to filter all your pools with absolute thresholds (ie a specific numbers instead of MADs from the median) for each of the QC metrics. 
You will have to change the thresholding values so that they are appropriate for your data. If you would like to instead filter by MADs from the median for each QC metric across all pools, you can jump to the :ref:`MAD Threshold section <mad_threshold_all-docs>`.

.. code-block:: R

  ### Load in required packages
  .libPaths("/usr/local/lib/R/site-library")  ## The location of R packages in the Singularity image - so R doesn't try to lead packages from outside the image
  library(Seurat)

  ### Read in the seurat object
  seurat <- readRDS("seurat_object_all_pools_singlet_barcodes_final_assignments.rds")

  ### Check the size of the seurat object
  print(seurat)

  ### Filter based on your absolute thresholds
  ## Substitute the values for nFeature_RNA (number of genes), nCount_RNA (number of UMI) and percent.mt as you need
  seurat_final <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 12 & nCount_RNA > 250)

  ### Check the size of the seurat object after filtering to ensure that cells have been removed
  print(seurat_final)

  ### Save the Seurat Object
  saveRDS(seurat_final, "/path/to/out/dir/seurat_singlets_QCfiltered.rds")


.. _mad_threshold_all-docs:

MAD Threshold
^^^^^^^^^^^^^^^^^

Here is some example code in order to filter by the same MAD from the median across all the pools in your dataset - you can choose a different MAD number for each QC metric. 
If you instead want to filter by different MAD thresholds in each pool, you can jump down to the :ref:`Filter Each Pool Separately: MAD Threshold section <mad_threshold_pools-docs>`.

.. code-block:: R

  ### Load in required packages
  .libPaths("/usr/local/lib/R/site-library")  ## The location of R packages in the Singularity image - so R doesn't try to lead packages from outside the image
  library(Seurat)

  ### Define a function to identify cells that are outliers based on certain MAD from the median
  mad_function <- function(seurat, column, number_mad){
      mad <- mad(seurat@meta.data[,column])
      low <- median(seurat@meta.data[,column]) - number_mad*mad
      high <- median(seurat@meta.data[,column]) + number_mad*mad
      print("The lower bound is:")
      print(low)
      print("The upper bound is:")
      print(high)
      seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
      return(seurat)
  }

  ### Read in the seurat object
  seurat <- readRDS("seurat_object_all_pools_singlet_barcodes_final_assignments.rds")

  ### Identify the cells that are outliers using MAD function
  seurat <- mad_function(seurat = seurat, column = "percent.mt", number_mad = 3)
  seurat <- mad_function(seurat = seurat, column = "nCount_RNA", number_mad = 3)
  seurat <- mad_function(seurat = seurat, column = "nFeature_RNA", number_mad = 3)

  ##### Remove the outliers #####
  seurat_final <- subset(seurat, subset = percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "NotOutlier" & nFeature_RNA_mad == "NotOutlier") 

  ### Check the size of the seurat object after filtering to ensure that cells have been removed
  print(seurat_final)

  ### Save the Seurat Object
  saveRDS(seurat_final, "/path/to/out/dir/seurat_singlets_QCfiltered.rds")


.. _filter_separate-docs:

Filter Each Pool Separately
---------------------------
If the QC figures reveal that some of your pools have quite different QC metric characteristics, it may be better to filter each of the pools with different thresholds - this can be decided through discussion with WG1.
If this is the case for your dataset, you can follow the sample code below.
As mentioned in the code below, a for loop would be beneficial for each of the QC metrics in order to reduce the risk of potential errors. 

Absolute Threshold
^^^^^^^^^^^^^^^^^^^^^^^
If you need to filter each of your pools on different QC metric thresholds and would like to use absolute thresholds (ie exact numbers as opposed to MAD from the median), you can use some of the example code below. 
If you have decided to use MAD from the median per pool instead, you can jump down to the :ref:`next section <mad_threshold_pools-docs>`.

.. code-block:: R

  .libPaths("/usr/local/lib/R/site-library")  ## The location of R packages in the Singularity image - so R doesn't try to lead packages from outside the image
  ### Load in required packages
  library(Seurat)

  ### Read in the seurat object
  seurat <- readRDS("seurat_object_all_pools_singlet_barcodes_final_assignments.rds")

  ### Switch the identity of the seurat object to the pools
  Idents(seurat) <- "Pool"

  ### Create a list to store seurat objects of each individual pool
  seurat_list <- list()

  ### Separate the one Seurat object into separate objects for each pool that are stored in a list
  for (pool in unique(seurat@meta.data$Pool)){
      seurat_list[[pool]] <- subset(seurat, idents = pool)
  }

  ### Then you can filter each pool based on the filtering thresholds selected
  seurat_list[["Pool1"]] <- subset(seurat_list[["Pool1"]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 12 & nCount_RNA > 250)
  seurat_list[["Pool2"]] <- subset(seurat_list[["Pool2"]], subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 15 & nCount_RNA > 300)
  ## ... continue this until all pools have been filtered. You can also put all the thresholding values into a dataframe, lists, or vectors and write a for loop to do this for all the pools

  ### Merge each of the seurat objects together again
  seurat_final <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

  ### Check the size of the seurat object after filtering to ensure that cells have been removed
  print(seurat_final)

  ### Save the Seurat Object
  saveRDS(seurat_final, "/path/to/out/dir/seurat_singlets_QCfiltered.rds")


.. _mad_threshold_pools-docs:

MAD Threshold
^^^^^^^^^^^^^^^
Here is some example code to filter by different MAD from the median per pool. Of course, you will have the change the number of MADs depending on each of the pools in your dataset. 

.. code-block:: R

  ### Load in required packages
  .libPaths("/usr/local/lib/R/site-library")  ## The location of R packages in the Singularity image - so R doesn't try to lead packages from outside the image
  library(Seurat)

  ### Define a function to identify cells that are outliers based on certain MAD from the median
  ## This will return a seurat object that contains new columns that indicate whether each cell is an outlier for that QC metric
  mad_function <- function(seurat, column, number_mad){
      mad <- mad(seurat@meta.data[,column])
      low <- median(seurat@meta.data[,column]) - number_mad*mad
      high <- median(seurat@meta.data[,column]) + number_mad*mad
      print("The lower bound is:")
      print(low)
      print("The upper bound is:")
      print(high)
      seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
      return(seurat)
  }

  ### Read in the seurat object
  seurat <- readRDS("seurat_object_all_pools_singlet_barcodes_final_assignments.rds")


  ### Switch the identity of the seurat object to the pools
  Idents(seurat) <- "Pool"

  ### Create a list to store seurat objects of each individual pool
  seurat_list <- list()

  ### Separate the one Seurat object into separate objects for each pool that are stored in a list
  for (pool in unique(seurat@meta.data$Pool)){
      seurat_list[[pool]] <- subset(seurat, idents = pool)
  }

  ### Calculate the MAD and cells to filter for each of the QC metrics of interest
  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "percent.mt", number_mad = 3)
  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "nCount_RNA", number_mad = 3)
  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "nFeature_RNA", number_mad = 3)

  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "percent.mt", number_mad = 2)
  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "nCount_RNA", number_mad = 2)
  seurat_list[["Pool1"]] <- mad_function(seurat = seurat_list[["Pool1"]], column = "nFeature_RNA", number_mad = 2)
  ## ... continue this until all pools have calculated all the MAD outliers for each QC metric . You can also put all the MAD thresholds into a dataframe, lists, or vectors and write a for loop to do this for all the pools

  ### Then you can filter each pool based on the filtering thresholds selected
  seurat_list[["Pool1"]] <- subset(seurat_list[["Pool1"]], subset = seurat_final <- subset(seurat, subset = percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "NotOutlier" & nFeature_RNA_mad == "NotOutlier")
  seurat_list[["Pool2"]] <- subset(seurat_list[["Pool2"]], subset = seurat_final <- subset(seurat, subset = percent.mt_mad == "NotOutlier" & nCount_RNA_mad == "NotOutlier" & nFeature_RNA_mad == "NotOutlier") 
  ## ... continue this until all pools have been filtered. You can also put all the thresholding values into a dataframe, lists, or vectors and write a for loop to do this for all the pools

  ### Check the size of the seurat object after filtering to ensure that cells have been removed
  print(seurat_final)

  ### Save the Seurat Object
  saveRDS(seurat_final, "/path/to/out/dir/seurat_singlets_QCfiltered.rds")


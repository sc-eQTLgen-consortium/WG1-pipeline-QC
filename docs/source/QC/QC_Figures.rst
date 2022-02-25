.. _QC_Figures-docs:

QC Figures
============

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)



In the last step, :ref:`Demultiplexing and Doublet Removal <Demultiplexing_Introduction-docs>`, you identified droplets that contained a single cell and assigned them to the correct individual. 
Next, we will use typical scRNA-seq QC metrics to assess the quality of those singlets and remove cells that don't pass those QC thresholds. 
For clarity, in the :ref:`Demultiplexing and Doublet Removal <Demultiplexing_Introduction-docs>` example we used just one pool. 
To look at the QC metrics figures, we will be using the results from two pools so we can easily see a comparison between the two pools. 
The QC metrics that the pipeline provides figures for by default are:

#. Mitochondrial percent

#. Ribosomal percent

#. Number UMIs

#. Number features

You will see three violin plots for each of these QC metrics. These have the pools on the x-axis and the QC metric on the y-axis. 

- One of the plots are the violin plots without any additional information:

  .. list-table::

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89234171-85347700-d62a-11ea-9d56-f1335eb4d126.png
          :width: 100%

          Number UMIs

      - .. figure:: https://user-images.githubusercontent.com/44268007/89234276-cd539980-d62a-11ea-9a36-192b64131a7f.png
          :width: 100%

          Number Genes Detected

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89235309-199fd900-d62d-11ea-9b45-bb43692d7227.png
          :width: 100%

          Ribosomal Percent

      - .. figure:: https://user-images.githubusercontent.com/44268007/89235400-3dfbb580-d62d-11ea-9207-f283fa0415b3.png
          :width: 100%

          Mitochondrial Percent




- One of the plots will have the lines indicating the median and different median absolute deviations (MADs) from the median for **all the pools together**. This will be useful to discuss the possibility of using a single filtering threshold across all pools. The median is in grey, 1 MAD is in blue, 2 MADs is in purple, 3 MADs is in red, 4 MADs is in orange and 5 MADs is in yellow:

  .. list-table::

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89236470-ca0edc80-d62f-11ea-9e1a-684e14139783.png
          :width: 100%

          Number UMIs

      - .. figure:: https://user-images.githubusercontent.com/44268007/89236475-ce3afa00-d62f-11ea-8a9b-ff5ed7b327e3.png
          :width: 100%

          Number Genes Detected

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89236397-9a5fd480-d62f-11ea-9a35-2a37c779dfaf.png
          :width: 100%

          Ribosomal Percent

      - .. figure:: https://user-images.githubusercontent.com/44268007/89236362-81572380-d62f-11ea-82cc-a24959fddbf2.png
          :width: 100%

          Mitochondrial Percent




- One of the plots will have the lines indicating the median and different median absolute deviations (MADs) from the median for **each pool separately**. This will be helpful to demonstrate possible pool-specific QC metrics and will be helpful in trying to decide if the pools shoudld each be filtered with a different threshold or if a common filtering threshold can be used for all the pools. The median is in grey, 1 MAD is in blue, 2 MADs is in purple, 3 MADs is in red, 4 MADs is in orange and 5 MADs is in yellow:
  
  .. list-table::

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89235743-f32e6d80-d62d-11ea-9c27-b95a8667080f.png
          :width: 100%

          Number UMIs

      - .. figure:: https://user-images.githubusercontent.com/44268007/89235781-03464d00-d62e-11ea-806b-3fe43fa673ab.png
          :width: 100%

          Number Genes Detected

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89243036-182cdb80-d642-11ea-9099-59c514085b8b.png
          :width: 100%

          Ribosomal Percent

      - .. figure:: https://user-images.githubusercontent.com/44268007/89243057-25e26100-d642-11ea-9f7e-9b27020ea184.png
          :width: 100%

          Mitochondrial Percent



Then there are a number of scatter plots that compare the QC metrics to one another:

- There will be two plots that compare number of features per cell vs the percent of mitochondrial content per cell where the cell is colored by pool (left) or includes lines that indicate the median and MADs from the median across all the cells for both the number of features per cell (x-axis) vs the percent of mitochondrial (y-axis). The median is in grey, 1 MAD is in blue, 2 MADs is in purple, 3 MADs is in red, 4 MADs is in orange and 5 MADs is in yellow.
  
  .. list-table::

    * - Colored by Pool:

      - With Median and MADs:

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89249983-c80a4500-d652-11ea-85ce-d66df6e1cb94.png
          :width: 100%

      - .. figure:: https://user-images.githubusercontent.com/44268007/89237825-5ec70980-d633-11ea-8691-5ba546c98ce4.png
          :width: 100%



      
- Two plots that compare the UMIs per cell (x-axis) vs the percent of mitochondrial genes per cell (y-axis) where the cell is colored by pool (left) or includes lines that indicate the median and MADs from the median across all the cells for both the UMIs per cell (x-axis) vs the percent of mitochondrial genes per cell (y-axis) (right). The median is in grey, 1 MAD is in blue, 2 MADs is in purple, 3 MADs is in red, 4 MADs is in orange and 5 MADs is in yellow.

  .. list-table::

    * - Colored by Pool:

      - With Median and MADs:

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89241690-bd45b500-d63e-11ea-993e-0784922f2273.png
          :width: 100%

      - .. figure:: https://user-images.githubusercontent.com/44268007/89241708-c8004a00-d63e-11ea-8de4-c3e1ef583a5e.png
          :width: 100%




- Two plots that compare the UMIs per cell (x-axis) vs the number of features per cell (y-axis) where the cell is colored by pool (left) or includes lines that indicate the median and MADs from the median across all the cells for both the UMIs per cell (x-axis) vs the percent of mitochondrial genes per cell (y-axis) (right). The median is in grey, 1 MAD is in blue, 2 MADs is in purple, 3 MADs is in red, 4 MADs is in orange and 5 MADs is in yellow.

  .. list-table::

    * - Colored by Pool:

      - With Median and MADs:

    * - .. figure:: https://user-images.githubusercontent.com/44268007/89242204-1f52ea00-d640-11ea-88b6-41e2b50f1b8c.png
          :width: 100%

      - .. figure:: https://user-images.githubusercontent.com/44268007/89242181-119d6480-d640-11ea-8d71-9585b25b0485.png
          :width: 100%

  



These figures can be used for discussion with the other members of WG1 of the sceQTL-Gen Consortium and should help identify relevant QC thresholds for this dataset. If you require additional QC figures, you can generate them using one of the seurat objects in the output directory. All the meta.data provided in each of these objects are provided toward the end of the :ref:`Demultiplexing and Doublet Removal <Demultiplexing_Introduction-docs>` section:

- ``seurat_object_all_pools_all_barcodes_all_metadata.rds`` to look at all the droplets from all pools and the assignments from each software

- ``seurat_object_all_pools_all_barcodes_final_assignments.rds`` to look at the droplets from all pools and the final assignment using the intersectional method

- ``seurat_object_all_pools_singlet_barcodes_final_assignments.rds`` to look at the singlets as defined by the intersectional method



Next Steps
------------

The next steps are to discuss the results with the consortium in order to identify filtering thresholds for your dataset.
Please reach out to Drew Neavin (d.neavin @ garvan.org.au) to set up a time for a discussion once you have completed these steps.
This process helps ensure consistency between different datasets used in the consortium.
After discussion with the consortium, you can move to the :ref:`QC Filtering section <Filtering>` which has pseudocode to help direct how to filter depending on the method selected in discussions with WG1.
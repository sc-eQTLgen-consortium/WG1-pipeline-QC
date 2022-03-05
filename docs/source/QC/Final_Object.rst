.. _QC_Final-docs:

Final Object
=============

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)



Final Seurat Object for WG2 Cell Prediction
---------------------------------------------

Congrats! The object you saved after :ref:`QC metric filtering <Filtering>` is ready to be used by working group two for cell classification. 
This final object is a seurat object that has all the genes (none of the genes have been filtered out yet) and the cells that were identified as singlets and passed QC metric filtering. 
There is also gene-level metadata and cell-level metadata saved in this object:

- Gene-level metadata is saved in the "RNA" assay of the seurat object and can be accessed with ``seurat[["RNA"]][[]]``. It contains the following columns:

  - Gene_ID: The HGNC symbol for this gene (Note that these were used for row names as well but since the Gene_IDs are not unique, they would have been made unique for the rownames)

  - ENSG_ID: The ENSG symbol for this gene

- Cell-level metadata can be accessed with ``seurat@meta.data``. It contains the following columns:

  - orig.ident: Name for the project. This will be "SeuratProject" by default

  - nCount_RNA: nUMIs per cell

  - nFeature_RNA: Number of genes per cell

  - Barcode: Original barcode ID; Note that this will be different from the rownames as seurat automatically changes the barcode names when multiple pools are merged together so that there are no issues with replicate barcode names across pools

  - Assignment: The droplet individual assignment from the intersectional method (individual ID, doublet or unassigned)

  - DropletType: The droplet type assignment from the intersectional method (singlet, doublet or unassigned)

  - Pool: Pool name that droplet was collected in   - percent.mt: Percent of genes that are mitochondrial per droplet

  - percent.rb: Percent of genes that are ribosomal per droplet

  - If you calculated MADs for filtering, columns indicating whether they are an outlier ("Outlier") or not ("NotOutlier") will also be present in the meta.data. Note that all "Outliers" should have been removed so only cells labled as "NotOutlier" should be present. 

These final QC-filterd seurat objects can now be integrated with the cell type classifications from `Working Group 2 <https://github.com/sc-eQTLgen-consortium/WG2-pipeline-classification>`__.

If you have any questions or issues while running the pipeline, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au).
# WG1-pipeline-QC
Part of the sc-eQTLgen consortium pipeline. Step 1, where the QC is done.

## Demultiplexing and Doublet Detection
### Background and Introduction
Demultiplexing and doublet detection has been built into a snakemake pipeline leveraging Singularity environments to maintain consistency across different compute structures. We identified five softwares to use for the sceQTL-Gen Consortium by testing the combination of 10 demultiplexing and doublet detection softwares with multiple intersectional methods (manuscript in process). The five softwares include two SNP-based demultiplexing and doublet detection softwares:

- [popscle-demuxlet](https://github.com/statgen/popscle/wiki)
- [souporcell](https://github.com/wheaton5/souporcell)

and three transcriptome-based doublet detection softwares:
- [DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
- [scds](https://github.com/kostkalab/scds)
- [scrublet](https://github.com/AllonKleinLab/scrublet)

The complete pipeline was built in Snakemake in order to provide reproducibility across labs and users with Singularity buckets to enable consistency in softwares across systems. Most of the softwares (popscle-demuxlet, souporcell and scds) run without need for any user interaction besides providing input files. However, both DoubletDetection and scrublet require users to check that the thresholding used is effective for the data. Therefore, the pipeline has been built to stop and only run certain jobs after users have provided input in the provided text files for each pool they are demultiplexing. Below are the intructions for running the demultiplexing and doublet detection software.

### Required Input

## QC Filtering


## Preparation for Working Group 2 - Cell Classification

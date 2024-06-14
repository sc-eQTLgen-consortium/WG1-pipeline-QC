# Changelog

All notable changes implemented in this branch compared to the main branch of [WG1-pipeline-QC](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC) (2022-03-12, SHA: `a9a81bd9b5aa12083dee47d2886a184d4a223196`) are documented in this file. 

Note that this branch is in beta and version 2.0.0 is not yet ready for release.


## [2.0.0] - Imputation - 2023-09-05

#### Additions
- Additional preflight checks to reduce failed runs
- Dynamic time usage per rule
- Support for VCF (per chromosome), as well as `PLINK` binary (per chromosome), as input
- Support for WGS data (`pre_processing.smk` steps + split per chromosome for `GenotypeHarmonizer`)
- Support for build 36 and 38 genotypes as input (as a result, also using build 38 1000G as reference for sex / ancestry check now)
- Functionality to split the output in datasets per ancestry (rule `split_per_dataset`)
- Added rule to combine phasing confidence stats
- Added option to ignore 1000G ancestry assignment and force an ancestry on all individuals

#### Fixes
- Fixed issue where some rules did not have dynamic time / memory usage
- Fixed issue where folders containing only numbers are read as integers instead of strings
- Fixed issue where missing values were not allowed in the PSAM for categorical variables
- Fixed issue where the number of samples printed was from the old PSAM instead of the updated PSAM
- Fixed issue where 1000G reference contains a duplicate variant (`X_68204280_rs1361839_C_T`). Fixed by switching to the build 38 reference
- Fixed issue where the MAF per ancestry selection was limited to two ancestries and any subsequent ancestries overwrote the setting of previous ancestries
- Fixed issue where `F_MISSING` was missing for rule `filter_preimpute_vcf` to filter on. Fixed by upgrading `bcftools` version to 1.18
- Fixed issue in `prune_1000g` where the variants in LD were extracted instead of the variants NOT in LD
- Fixed issue in `filter_het.R` where the heterozygosity passed samples included sampels below -N standard deviations


#### Changes
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / ouput of rules more explicit
- Moved all filepaths to settings file
- Moved all settings to settings file
- Reduced amount of temporary files
- Mark temporary files as `temp()` so they will be removed after use
- Moved demultiplexing rules to `Demultiplexing` pipeline
- Switched from singularity to Docker
- Updated all software versions
- Merged singularity image with `Demultiplexing`
- Moved (reference) data download from image to separate download scripts `download_data.sh` / `download_test_dataset.sh`
- Moved scripts outside of image for flexibility
- Reduced image size
- Changed `check_sex` rule to only use `X` chromosomes to reduce file size 
- Moved `update_sex_ancestry` rule to `Snakefile`
- Moved `genotype_donor_annotation` rule to `Snakefile`
- Moved `crossmap` rule before the sex and ancestry check
- Merged `indivindiv_missingness` rule into the `pre_processing.smk` / `wgs_filter.smk` rules
- Split `common_snps`, `prune_1000g`  rules into `prune_1000g`, doing all the pruning of the 1000G reference, and then only ones filtering the input data in `prune_data`
- Translated bash commands in `prune_1000g` into Python script `overlap_pvars.py` for clarity and robustness
- Removed `sort_bed` rule; functionality is now part of `crossmap`
- Removed `final_pruning` rule; duplicates were already removed earlier
- Changed `calculate_missingness` and `het` to use gzipped VCF files
- Moved creation of `sex_update_remove.tsv` to `check_sex` rule such that rule `pca_projection_assign` only relates to the ancestry check, and `check_sex` does all the sex check
- Split `qc_plots` rule into `combine_pools` (rds output) and `qc_plots` (figures) rules

## [2.0.0] - Demultiplexing - 2023-09-05

#### Additions
- Added additional preflight checks to reduce failed runs
- Added `DoubletFinder` doublet detection method
- Added `scDblFinder` doublet detection method
- Added Flag `is_multiplexed` to skip demultiplexing
- Added Flag `sc_data_type` denoting the single-cell data type to automatically select doublet detection method
- Added rule `rename_chrs` to enable the input of BAM files which have chromosomes encoded as `chr1` instead of `1`
- Added rule `reheader_vcf` to prevent different ordering of chromosomes errors
- Added fixes by Roy Oelen and Dan Kaptijn improving performance of `souporcell` including reduced memory usage and allowing gzipped barcodes as input
- Added additional fix build upon work from Roy Oelen and Dan Kaptijn to make `souporcell` also use the filtered BAM file
- Added parallel (re-)processing of all method with fully dynamic settings per run without overwriting previous results
- Added `plot_DoubletDetection`, `plot_DoubletFinder` and `plot_Scrublet` rule to create combined plots for multiple runs 
- Added demultiplexing rules from the `Imputation` pipeline
- Added adaptive `combine_results` rule from [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/v0.0.4/) v0.0.4 using the results of whichever doublet detection methods was used
- Split `souporcell_pipeline.py` into separate rules and refactored to use gzipped files for computational efficiency
- Added `skip_remap` for `souporcell`
- Made parameter `expected_doublet_rate` for `DoubletFinder` and `Scrublet` variable dependend on expected doublet rate
- Added rule `verifyBamID`, optionally using the `remap` functionality of `souporcell`, to check for sample swaps
- Added multiple checks verifying that output files are not empty to prevent running subsequent rules 

#### Fixes
- Fixed issue where some rules did not have dynamic time / memory usage
- Fixed issue where naive BAM file search can return temporary bam file in `CellRanger` v7.0.1
- Fixed issue where `Scrublet` uses `min_cells` value as parameter for `min_counts`
- Fixed issue where some R scripts had hard-coded `future.globals.maxSize`
- Fixed issue in `scds.R` where `bcds()` would fail if some barcodes had zero counts for the selected genes
- Fixed issue in `Assign_Indiv_by_Geno.R` where `Heatmap` would fail if there was only a single cluster / individual in the pool
- Fixed edge case where the Demuxafy `combine_results` rule fails if only one doublet detection method was run
- Fixed issue in `Singlet_QC_Figures.R` where gene symbols were stored as ENSG
- Fixed issue in `Singlet_QC_Figures.R` where combining pools from CellRanger v>=6.0.0 (i.e. containing variable numer of genes) results in a `number of rows of matrices must match` error
- Fixed issue in `expected_observed_individuals_doublets.R` where the figure width exceeds the maximum if there are too many pools

#### Changes
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / ouput of rules more explicit
- Moved all filepaths to settings file
- Moved all settings to settings file
- Reduced amount of temporary files
- Mark temporary files as `temp()` so they will be removed after use
- Switched from singularity to Docker
- Updated all software versions
- Merged singularity image with `Imputation`
- Moved (reference) data download from image to separate download scripts `download_data.sh` / `download_test_dataset.sh`
- Moved scripts outside of image for flexibility
- Reduced image size
- Using official souporcell conda environment
- Doublet detection methods now use the `h5` file as input (required to align with new `WG0-Propocessing` pipeline)
- Changed scripts to align with [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/v0.0.4/) v0.0.4 code
- Update rules `filter4demultiplexing`, `sort4demultiplexing`, `popscle_pileup`, `popscle_demuxlet`, `souporcell` to use / output gzipped VCF files
- Made all files gzipped
- Moved merging of the pools and related rule `barcode_qc_plots` to WG3, only keep merging of the droplet assignments metadata in this WG

#### Known issues
- option `assignment` in rule `combine_results`  (i.e. checking for sample swaps) only works if there is no demultiplexing method applied.
- due to the way the input_files for rule `all` depend on the status of the pipeline (multiple runs), options like `--delete-temp-output` to remove old temp files will not work since old rules are not rechecked
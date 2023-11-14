# Changelog

All notable changes implemented in this branch compared to the main branch of [sc-eQTLgen-consortium](https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC) (2022-03-12, SHA: `a9a81bd9b5aa12083dee47d2886a184d4a223196`) are documented in this file. 

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

#### Fixes
- Fixed issue where some rules did not have dynamic time / memory usage
- Fixed issue where folders containing only numbers are read as integers instead of strings
- Fixed issue where missing values were not allowed in the PSAM for categorical variables
- Fixed issue where the number of samples printed was from the old PSAM instead of the updated PSAM
- Fixed issue where 1000G reference contains a duplicate variant (`X_68204280_rs1361839_C_T`). Fixed by switching to the build 38 reference
- Fixed issue where the MAF per ancestry selection was limited to two ancestries and any subsequent ancestries overwrote the setting of previous ancestries
- Fixed issue where `F_MISSING` was missing for rule `filter_preimpute_vcf` to filter on. Fixed by upgrading `bcftools` version to 1.18
- Fixed issue in `prune_1000g` where the variants in LD were extracted instead of the variants NOT in lD

#### Changes
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / ouput of rules more explicit
- Moved all filepaths to settings file
- Moved all settings to settings file
- Reduced amount of temporary files
- Moved demultiplexing rules to `Demultiplexing` pipeline
- Switched from singularity to Docker
- Updated all software versions
- Merged singularity image with `Demultiplexing`
- Moved (reference) data download from image to separate download scripts `download_data.sh` / `download_test_dataset.sh`
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


## [2.0.0] - Demultiplexing - 2023-09-05

#### Additions
- Added additional preflight checks to reduce failed runs
- Added `DoubletFinder` doublet detection method
- Added `scDblFinder` doublet detection method
- Added Flag `is_multiplexed` to skip demultiplexing
- Added parallel (re-)processing of `DoubletDetection` and `Scrublet` runs with fully dynamic settings per run without overwriting previous results
- Added Flag `sc_data_type` denoting the single-cell data type to automatically select doublet detection method
- Added fixes by Roy Oelen and Dan Kaptijn improving performance of `souporcell` including reduced memory usage and allowing gzipped barcodes as input
- Added demultiplexing rules from the `Imputation` pipeline

#### Fixes
- Fixed issue where some rules did not have dynamic time / memory usage
- Fixed issue where naive BAM file search can return temporary bam file in `CellRanger` v7.0.1
- Fixed issue where `Scrublet` uses `min_cells` value as parameter for `min_counts`

#### Changes
- Refactor code to (mostly) PEP8
- Improved clarity of rules by making input / ouput of rules more explicit
- Moved all filepaths to settings file
- Moved all settings to settings file
- Reduced amount of temporary files
- Switched from singularity to Docker
- Updated all software versions
- Merged singularity image with `Imputation`
- Moved (reference) data download from image to separate download scripts `download_data.sh` / `download_test_dataset.sh`
- Reduced image size
- Doublet detection methods now use the `h5` file as input (required to align with new `WG0-Propocessing` pipeline)
- Changed scripts to align with `[Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/v0.0.4/) v0.0.4 code
- Added adaptive `combine_results` rule using the results of whichever doublet detection methods was used
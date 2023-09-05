#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/scds.R
.libPaths("/usr/local/lib/R/site-library")

suppressMessages(suppressWarnings(library(argparse)))
# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-o", "--out", required=TRUE, type="character", help="The output directory where results will be saved.")
parser$add_argument("-c", "--counts", required=TRUE, type="character", help="Path to the 10x filtered h5 file.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(scds)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

## Add max future globals size for large pools
options(future.globals.maxSize=(850 * 1024 ^ 2))

## Read in data
counts <- Read10X_h5(args$counts)

if (is.list(counts)){
	sce <- SingleCellExperiment(list(counts=counts[[grep("Gene", names(counts))]]))
} else {
	sce <- SingleCellExperiment(list(counts=counts))
}

## Annotate doublet using binary classification based doublet scoring:
sce <- bcds(sce, retRes=TRUE, estNdbl=TRUE)

## Annotate doublet using co-expression based doublet scoring:
try({
    sce <- cxds(sce, retRes=TRUE, estNdbl=TRUE)
})

### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
} else {
    print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
}

## Doublet scores are now available via colData:
colnames(Doublets) <- c("Barcode", "scds_score", "scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE", "singlet", Doublets$scds_DropletType)
Doublets$scds_DropletType <- gsub("TRUE", "doublet", Doublets$scds_DropletType)

message(paste0("Writing results to ", args$out, "/scds_doublets_singlets.tsv."))
write_delim(Doublets, paste0(args$out, "/scds_doublets_singlets.tsv"), "\t")


summary <- as.data.frame(table(Doublets$scds_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
message(paste0("Writing summary to ", args$out, "/scds_doublet_summary.tsv."))
write_delim(summary, paste0(args$out, "/scds_doublet_summary.tsv"), "\t")


#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/scDblFinder.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-o", "--out", required=TRUE, help="The output directory where results will be saved")
parser$add_argument("-c", "--counts", required=TRUE, type="character", help="Path to the 10x filtered h5 file.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(scDblFinder)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(tidyverse)))

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

## Calculate doublet ratio ###
doublet_ratio <- ncol(sce) / 1000 * 0.008

### Calculate Singlets and Doublets ###
sce <- scDblFinder(sce, dbr=doublet_ratio)

### Make a dataframe of the results ###
results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)

message(paste0("Writing results to ", args$out, "/scDblFinder_doublets_singlets.tsv."))
write_delim(results, file=paste0(args$out,"/scDblFinder_doublets_singlets.tsv"), delim="\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(results$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
message(paste0("Writing summary to ", args$out, "/scDblFinder_doublet_summary.tsv."))
write_delim(summary, paste0(args$out, "/scDblFinder_doublet_summary.tsv"), "\t")


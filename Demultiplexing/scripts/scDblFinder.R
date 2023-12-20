#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/scDblFinder.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(jsonlite)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--counts", required=TRUE, type="character", help="Path to the 10x filtered h5 file.")
parser$add_argument("--expected_doublet_scaling_factor", required=FALSE, default=8e-06, type="double", help="The fraction of droublets expected based on the number of nuclei recovered.")
parser$add_argument("--stdev_doublet_rate", required=FALSE, default=NULL, type="double", help="The uncertainty range in the doublet rate, interpreted as a +/- around `dbr`. During thresholding, deviation from the expected doublet rate will be calculated from these boundaries, and will be considered null within these boundaries. If NULL, will be 40% of `dbr`. Set to `dbr.sd=0` to disable the uncertainty around the doublet rate, or to `dbr.sd=1` to disable any expectation of the number of doublets (thus letting the thresholding be entirely driven by the misclassification of artificial doublets).")
parser$add_argument("--nfeatures", required=FALSE, default=1352, type="integer", help="The number of top features to use.")
parser$add_argument("--dims", required=FALSE, default=20, type="integer", help="The number of dimensions used.")
parser$add_argument("--keepUnidentifiable", dest="removeUnidentifiable", action='store_false', default=TRUE, help="Whether to remove artificial doublets of a combination that is generally found to be unidentifiable.")
parser$add_argument("--include_pcs", required=FALSE, default=19, type="integer", help="The number of top components to use (e.g. `includePCs=10`, equivalent to 1:10).")
parser$add_argument("--prop_markers", required=FALSE, default=0, type="double", help="The proportion of features to select based on marker identification.")
parser$add_argument("--score", required=FALSE, choices=c("xgb", "weighted", "ratio"), default="xgb", type="character", help="Score to use for final classification.")
parser$add_argument("--processing", required=FALSE, choices=c("default", "rawPCA", "rawFeatures", "normFeatures"), default="default", type="character", help="processing Counts (real and artificial) processing before KNN. Either 'default' (normalization and PCA), 'rawPCA' (PCA without normalization), 'rawFeatures' (no normalization/dimensional reduction), 'normFeatures' (uses normalized features, without PCA), returning a named matrix with cells as rows and components as columns.")
parser$add_argument("--metric", required=FALSE, choices=c("merror", "logloss", "auc", "aucpr"), default="logloss", type="character", help="Error metric to optimize during training.")
parser$add_argument("--nrounds", required=FALSE, default=0.25, type="double", help="Maximum rounds of boosting.")
parser$add_argument("--max_depth", required=FALSE, default=4, type="integer", help="Maximum depths of each tree.")
parser$add_argument("--iter", required=FALSE, default=3, type="integer", help="A positive integer indicating the number of scoring iterations (ignored if `score` isn't based on classifiers). At each iteration, real cells that would be called as doublets are excluding from the training, and new scores are calculated. Recommended values are 1 or 2.")
parser$add_argument("--multi_sample_mode", required=FALSE, choices=c("split", "singleModel", "singleModelSplitThres", "asOne"), default="split", type="character", help="Either 'split' (recommended if there is heterogeneity across samples), 'singleModel', 'singleModelSplitThres', or 'asOne'")
parser$add_argument("--mem", required=FALSE, default=500, type="integer", help="The maximum allowed size in GB")
parser$add_argument("--out", required=TRUE, help="The output directory where results will be saved.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

print("Options in effect:")
for (name in names(args)) {
	print(paste0("  --", name, " ", args[[name]]))
}
print("")

write(toJSON(args, pretty=TRUE, auto_unbox=TRUE, null="null"), paste0(args$out,"/scDblFinder_settings.json"))

suppressMessages(suppressWarnings(library(scCustomize)))
suppressMessages(suppressWarnings(library(scDblFinder)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(tidyverse)))

## Add max future globals size for large pools
options(future.globals.maxSize=(args$mem * 1000 * 1024^2))

## Read in data
counts <- tryCatch({
	print("Loading count matrix using Seurat - Read10X_h5()")
	counts <- Read10X_h5(args$counts)
},error = function(e){
	print("Failed, trying to load count matrix using scCustomize - Read_CellBender_h5_Mat()")
	counts <- Read_CellBender_h5_Mat(args$counts)
})

if (is.list(counts)){
	counts <- counts[[grep("Gene", names(counts))]]
}

paste0('Counts matrix shape: ', nrow(counts) ,' rows, ', ncol(counts), ' columns')

sce <- SingleCellExperiment(list(counts=counts))

### Calculate Singlets and Doublets ###
# Default: https://github.com/plger/scDblFinder/blob/devel/R/scDblFinder.R
# clusters=NULL
# samples=NULL
# clustCor=NULL
# artificialDoublets=NULL,
# knownDoublets=NULL
# knownUse="discard"
# dbr=NULL
# dbr.sd=NULL
# nfeatures=1000
# dims=20
# k=NULL
# removeUnidentifiable=TRUE
# includePCs=10
# propRandom=0
# propMarkers=0
# aggregateFeatures=FALSE,
# returnType="sce"
# score="xgb"
# processing="default"
# metric="logloss"
# nrounds=0.25
# max_depth=4
# iter=3
# trainingFeatures=NULL
# multiSampleMode=c("split","singleModel","singleModelSplitThres","asOne")
# threshold=TRUE
# verbose=is.null(samples)
# NOTE: defaults of nfeatures changed from 1000 to 1352 and includePCs changede from 10 to 19 on Nov 14, 2022.
# However we use an older version so I use those defaults.
sce <- scDblFinder(
	sce,
	clusters=NULL,
	samples=NULL,
	clustCor=NULL,
	artificialDoublets=NULL,
	knownDoublets=NULL,
	knownUse="discard",
	dbr=ncol(sce) * args$expected_doublet_scaling_factor,
	dbr.sd=args$stdev_doublet_rate,
	nfeatures=args$nfeatures,
	dims=args$dims,
	k=NULL,
	removeUnidentifiable=args$removeUnidentifiable,
	includePCs=args$include_pcs,
	propRandom=0,
	propMarkers=args$prop_markers,
	aggregateFeatures=FALSE,
	returnType="sce",
	score=args$score,
	processing=args$processing,
	metric=args$metric,
	nrounds=args$nrounds,
	max_depth=args$max_depth,
	iter=args$iter,
	trainingFeatures=NULL,
	multiSampleMode=args$multi_sample_mode,
	threshold=TRUE,
	verbose=TRUE
)

### Make a dataframe of the results ###
results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score)

message(paste0("Writing results to ", args$out, "scDblFinder_doublets_singlets.tsv.gz."))
write_delim(results, file=gzfile(paste0(args$out, "scDblFinder_doublets_singlets.tsv.gz")), delim="\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(results$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
message(paste0("Writing summary to ", args$out, "scDblFinder_doublet_summary.tsv.gz."))
write_delim(summary, gzfile(paste0(args$out, "scDblFinder_doublet_summary.tsv.gz")), "\t")


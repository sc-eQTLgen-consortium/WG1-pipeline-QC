#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/scds.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(jsonlite)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--counts", required=TRUE, type="character", help="Path to the 10x filtered h5 file.")
parser$add_argument("--bcds_ntop", required=FALSE, default=500, type="integer", help="Indicating number of top variance genes to consider.")
parser$add_argument("--bcds_srat", required=FALSE, default=1, type="integer", help="indicating ratio between orginal number of 'cells' and simulated doublets.")
parser$add_argument("--bcds_nmax", required=FALSE, default="tune", type="character", help="maximum number of training rounds; integer or 'tune'")
parser$add_argument("--cxds_ntop", required=FALSE, default=500, type="integer", help="Indimessageing number of top variance genes to consider.")
parser$add_argument("--cxds_binthresh", required=FALSE, default=0, type="integer", help="Minimum counts to consider a gene 'present' in a cell.")
parser$add_argument("--out", required=TRUE, type="character", help="The output directory where results will be saved.")

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

write(toJSON(args, pretty=TRUE, auto_unbox=TRUE, null="null"), paste0(args$out,"/scds_settings.json"))

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(scds)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

## Add max future globals size for large pools
options(future.globals.maxSize=(850 * 1024 ^ 2))

## Read in data
counts <- Read10X_h5(args$counts)

paste0('Counts matrix shape: ', nrow(counts) ,' rows, ', ncol(counts), ' columns')

if (is.list(counts)){
	sce <- SingleCellExperiment(list(counts=counts[[grep("Gene", names(counts))]]))
} else {
	sce <- SingleCellExperiment(list(counts=counts))
}

## Annotate doublet using binary classification based doublet scoring:
# Defaults: https://github.com/kostkalab/scds/blob/master/R/bcds.R
# ntop=500
# srat=1
# verb=FALSE
# retRes=FALSE
# nmax="tune"
# varImp=FALSE
# estNdbl=FALSE
sce <- bcds(
  sce,
  ntop=args$bcds_ntop,
  srat=args$bcds_srat,
  verb=TRUE,
  retRes=TRUE,
  nmax=args$bcds_nmax,
  varImp=FALSE,
  estNdbl=TRUE
)

## Annotate doublet using co-expression based doublet scoring:
# Defaults: https://github.com/kostkalab/scds/blob/master/R/cxds.R
# ntop=500
# binThresh=0
# verb=FALSE
# retRes=FALSE
# estNdbl=FALSE
try({
    sce <- cxds(
      sce,
      ntop=args$cxds_ntop,
      binThresh=args$cxds_binthresh,
      verb=TRUE,
      retRes=TRUE,
      estNdbl=TRUE)
})

### If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    # Defaults: https://github.com/kostkalab/scds/blob/master/R/cxds_bcds_hybrid.R
    # cxdsArgs=NULL
    # bcdsArgs=NULL
    # verb=FALSE
    # estNdbl=FALSE
    # force=FALSE
    sce <- cxds_bcds_hybrid(
      sce,
      cxdsArgs=NULL,
      bcdsArgs=NULL,
      verb=TRUE,
      estNdbl=TRUE,
      force=FALSE
    )
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
} else {
    print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
}

## Doublet scores are now available via colData:
colnames(Doublets) <- c("Barcode", "scds_score", "scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE", "singlet", Doublets$scds_DropletType)
Doublets$scds_DropletType <- gsub("TRUE", "doublet", Doublets$scds_DropletType)

message(paste0("Writing results to ", args$out, "scds_doublets_singlets.tsv.gz."))
write_delim(Doublets, gzfile(paste0(args$out, "scds_doublets_singlets.tsv.gz")), "\t")

summary <- as.data.frame(table(Doublets$scds_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
message(paste0("Writing summary to ", args$out, "scds_doublet_summary.tsv."))
write_delim(summary, paste0(args$out, "scds_doublet_summary.tsv"), "\t")


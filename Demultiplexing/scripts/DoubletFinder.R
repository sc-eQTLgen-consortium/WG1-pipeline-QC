#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/DoubletFinder.R and
# https://github.com/chris-mcginnis-ucsf/DoubletFinder and https://rpubs.com/kenneditodd/doublet_finder_example

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-o", "--out", required=TRUE, help="The output directory where results will be saved.")
parser$add_argument("-c", "--counts", required=TRUE, type="character", help="Path to the 10x filtered h5 file.")
parser$add_argument("-r", "--resolution", required=FALSE, default=0.8, type="double", help="Resultion to use for \'FindClusters()' function.")
parser$add_argument("-n", "--pN", required=FALSE, default=0.25, type="double", help="Number of doublets to simulate as a proportion of the pool size.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(scCustomize)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(DoubletFinder)))
suppressMessages(suppressWarnings(library(readr)))

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

## Add max future globals size for large pools
options(future.globals.maxSize=(850 * 1024 ^ 2))

## Read in data
counts <- Read10X_h5(args$counts)

## Construct Seurat object.
seu <- CreateSeuratObject(counts)

## Pre-process seurat object with standard seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, nfeatures.print=10)

# Find significant PCs
stdv <- seu[["pca"]]@stdev
sum.stdv <- sum(seu[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing=TRUE)[1] + 1
min.pc <- min(co1, co2)
paste0("Number of statistically-significant principal components: ", min.pc)

# finish pre-processing
seu <- RunUMAP(seu, dims=1:min.pc)
seu <- FindNeighbors(object=seu, dims=1:min.pc)
seu <- FindClusters(object=seu, resolution=find.clusters.resolution) # resolution = 0.8 is default, adds "RNA_snn_res.0.1" and "seurat_clusters" to meta.data

# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(seu, PCs=1:min.pc) # generates Rplots.pdf
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
plot <- ggplot(bcmvn, aes(pK, BCmetric)) +  geom_point()
ggsave(plot, filename=paste0(args$out,"/pKvBCmetric.png"))

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
paste0("PC neighborhood size used to compute pANN: ", optimal.pk)

## Homotypic doublet proportion estimate
annotations <- seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp.poi <- round(optimal.pk * nrow(seu@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
paste0("pANN therhold used to make final doublet/singla prediction: ", nExp.poi.adj)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
# seu ~ This is a fully-processed Seurat object (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).
# PCs ~ The number of statistically-significant principal components, specified as a range (e.g., PCs = 1:10)
# pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).
# pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values should be estimated using the strategy described below.
# nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions. This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.
seu <- doubletFinder_v3(seu=seu, PCs=1:min.pc, pN=pN, pK=optimal.pk, nExp=nExp.poi.adj) # adds "pANN_0.25_0.05_277" and "DF.classifications_0.25_0.05_277" to meta.data
doublets <- as.data.frame(seu@meta.data[grepl("pANN*|DF.classifications*", colnames(seu@meta.data))])
colnames(doublets) <-  c("DoubletFinder_score", "DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet", "singlet", doublets$DoubletFinder_DropletType) %>% gsub("Doublet", "doublet",.)

message(paste0("Writing results to ", args$out, "/DoubletFinder_doublets_singlets.tsv."))
write_delim(doublets, file = paste0(args$out, "/DoubletFinder_doublets_singlets.tsv"), delim="\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
message(paste0("Writing summary to ", args$out, "/DoubletFinder_doublet_summary.tsv."))
write_delim(summary, paste0(args$out, "/DoubletFinder_doublet_summary.tsv"), "\t")


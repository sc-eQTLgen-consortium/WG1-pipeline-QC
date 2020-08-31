.libPaths("/usr/local/lib/R/site-library")
library(dplyr)
library(tidyr)
library(tidyverse)
library(scds)
library(Seurat)
library(SingleCellExperiment)

args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
pool <- arguments[1,]
out <- arguments[2,]
tenX <- arguments[3,]

## Read in data
counts <- Read10X(as.character(tenX[1]), gene.column = 1)
sce <- SingleCellExperiment(list(counts=counts))

## Annotate doublet using co-expression based doublet scoring:
sce = cxds(sce, retRes = TRUE)

## Annotate doublet using binary classification based doublet scoring:
sce = bcds(sce, retRes = TRUE)

## Combine both annotations into a hybrid annotation
sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)

## Doublet scores are now available via colData:
Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType) 
Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)

message("writing output")
write_delim(Doublets, paste0(out,"/scds_doublets.txt"), "\t")
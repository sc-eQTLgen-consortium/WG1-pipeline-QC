#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/scripts/Singlet_QC_Figures.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-p", "--poolsheet", required=TRUE, help="")
parser$add_argument("-m", "--main_dir", required=TRUE, help="")
parser$add_argument("-rb", "--rb_genes", required=TRUE, help="")
parser$add_argument("-mt", "--mt_genes", required=TRUE, help="")
parser$add_argument("-o", "--out", required=TRUE, help="The output directory where results will be saved.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

## make sure the directory exists ###
dir.create(args$out, recursive = TRUE)

print("Options in effect:")
for (name in names(args)) {
	print(paste0("  --", name, " ", args[[name]]))
}
print("")

suppressMessages(suppressWarnings(library(scCustomize)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggforce)))
suppressMessages(suppressWarnings(library(ggnewscale)))

pools <- read_delim(args$poolsheet, delim = "\t")

##### Read in the mt and rb gene lists #####
RB_genes <- read_delim(args$rb_genes, delim = "\t")
MT_genes <- read_delim(args$mt_genes, delim = "\t")

##### Read in and format the software assignments #####
load_seurat_object <- function(row, assignment_path){
    print(paste0("  Creating Seurat object for pool ", row[["Pool"]]))
    counts <- tryCatch({
        Read10X_h5(row[["Counts"]])
    },error = function(e){
        Read_CellBender_h5_Mat(row[["Counts"]])
    })
    print(paste0("    Loaded matrix with shape (", length(rownames(counts)), ", ", length(colnames(counts)), ")"))

    assignment <- as.data.frame(read_delim(paste0(as.character(args$main_dir), row[["Pool"]], assignment_path), delim = "\t"))
    print(paste0("    Loaded assignments with shape (", length(rownames(assignment)), ", ", length(colnames(assignment)), ")"))

    assignment$Pool <- row[["Pool"]]
    assignment$Barcode <- gsub("-1", "", assignment$Barcode)
    assignment$Barcode <- paste0(assignment$Barcode, "_", row[["Pool"]])
    rownames(assignment) <- assignment$Barcode

    # change to a barcode unique across lanes
    colnames(counts) <- assignment$Barcode

    seurat <- CreateSeuratObject(counts, min.cells = 0, min.features = 0, meta.data = assignment)
    print("  Done")
    return(seurat)
}
paste0("Creating Seurat object with all metadata")
seurat_objects <- apply(pools,1, load_seurat_object, assignment_path="/CombinedResults/combined_results_w_combined_assignments.tsv")
seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
saveRDS(seurat, paste0(args$out, "seurat_object_all_pools_all_barcodes_all_metadata.rds"))
paste0("Saved ", args$out, "seurat_object_all_pools_all_barcodes_all_metadata.rds")

paste0("Creating Seurat object with final assignments")
seurat_objects <- apply(pools,1, load_seurat_object, assignment_path="/CombinedResults/Final_Assignments_demultiplexing_doublets.tsv")
seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
saveRDS(seurat, paste0(args$out, "seurat_object_all_pools_all_barcodes_final_assignments.rds"))
paste0("Saved ", args$out, "seurat_object_all_pools_all_barcodes_final_assignments.rds")

##### subset for just singlets #####
paste0("Creating Seurat object with final assignments of the singlets")
seurat <- subset(seurat, DropletType == "singlet")

##### Get the mitochondiral and ribosomal percentage QC metrics for each cell #####
if ((sum(which(MT_genes$GeneID %in% rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat)))) & (sum(which(MT_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    mt_features <- MT_genes$GeneID[MT_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
} else if ((sum(which(MT_genes$ENSG %in% rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (sum(which(MT_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    mt_features <- MT_genes$ENSG[MT_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mt_features)
} else if ((length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have mitochondrial genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

if ((sum(which(RB_genes$GeneID %in% rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat)))) & (sum(which(RB_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    rb_features <- RB_genes$GeneID[RB_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
} else if ((sum(which(RB_genes$ENSG %in% rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (sum(which(RB_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    rb_features <- RB_genes$ENSG[RB_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rb_features)
} else if ((length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have ribosomal genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

##### Save the object with singlets and QC metrics #####
saveRDS(seurat, paste0(args$out, "seurat_object_all_pools_singlet_barcodes_final_assignments.rds"))
paste0("Saved ", args$out, "seurat_object_all_pools_singlet_barcodes_final_assignments.rds")

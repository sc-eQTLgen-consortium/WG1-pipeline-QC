#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

suppressMessages(suppressWarnings(library(argparse)))
# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-b", "--basedir", required=TRUE, help="")
parser$add_argument("-p", "--pools", required=TRUE, help="")
parser$add_argument("-i", "--pool_indirs", required=TRUE, help="")
parser$add_argument("-o", "--out", required=TRUE, help="")
parser$add_argument("-r", "--rb_genes_file", required=TRUE, help="")
parser$add_argument("-m", "--mt_genes_file", required=TRUE, help="")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggforce)))
suppressMessages(suppressWarnings(library(ggnewscale)))

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

## Read in data
pools <- read_delim(as.character(args$pools[1]), delim="\t", col_names=c("Pool"))
pools_list <- sapply(pools$Pool, toString)

##### Read in the mt and rb gene lists #####
RB_genes <- read_delim(args$rb_genes_file, delim="\t")
MT_genes <- read_delim(args$mt_genes_file, delim="\t")

##### Readin in the file directory locations #####
dir_locations <- read_delim(args$pool_indirs, delim="\t")

##### Read in the genes file
if (file.exists(paste0(dir_locations$Matrix_Directories[1], "/genes.tsv"))){
    genes <- read_delim(paste0(dir_locations$Matrix_Directories[1], "/genes.tsv"), delim="\t", col_names=c("ENSG_ID", "Gene_ID"))
} else if (file.exists(paste0(dir_locations$Matrix_Directories[1], "/features.tsv.gz"))){
    genes <- read_delim(paste0(dir_locations$Matrix_Directories[1], "/features.tsv.gz"), delim="\t", col_names=c("ENSG_ID", "Gene_ID", "FeatureType"))
} else {
    print(print("We're having issues finding your gene file (gene.tsv or features.tsv.gz). Please make sure one of these exist in ", dir_locations$Matrix_Directories[1]))
}


##### Read in and format the software assignments #####
assignments_list <- lapply(pools_list, function(x){
    read_delim(paste0(as.character(args$basedir), "/", x, "/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"), delim="\t")
})
names(assignments_list) <- pools_list

assignments_list <- lapply(pools_list, function(x){
    assignments_list[[x]]$Pool <- x
    assignments_list[[x]]$Barcode <- gsub(" - 1", "", assignments_list[[x]]$Barcode) 
    assignments_list[[x]]$Barcode  <- paste0(assignments_list[[x]]$Barcode, "_", x)
    return(assignments_list[[x]])
})
assignments <- do.call(rbind,assignments_list)
assignments <- as.data.frame(assignments)
rownames(assignments) <- assignments$Barcode

##### Read in and format the final assignments #####
assignments_final_list <- lapply(pools_list, function(x){
    read_delim(paste0(as.character(args$basedir), "/", x, "/CombinedResults/Final_Assignments_demultiplexing_doublets.txt"), delim="\t")
})
names(assignments_final_list) <- pools_list

assignments_final_list <- lapply(pools_list, function(x){
    assignments_final_list[[x]]$Pool <- x
    assignments_final_list[[x]]$Barcode <- gsub("-1", "", assignments_final_list[[x]]$Barcode) 
    assignments_final_list[[x]]$Barcode  <- paste0(assignments_final_list[[x]]$Barcode, "_", x)
    return(assignments_final_list[[x]])
})
assignments_final <- do.call(rbind,assignments_final_list)
assignments_final <- as.data.frame(assignments_final)
rownames(assignments_final) <- assignments_final$Barcode

## Read in data
counts_list <- lapply(dir_locations$Matrix_Directories, function(x){
    Read10X(x, gene.column=2)
})
names(counts_list) <- pools_list

counts_list <- lapply(pools_list, function(x){
    colnames(counts_list[[x]]) <- gsub(" - 1", "", colnames(counts_list[[x]]))
    colnames(counts_list[[x]]) <- paste0(colnames(counts_list[[x]]), "_", x)
    return(counts_list[[x]])
})
counts <- do.call(cbind, counts_list)

##### Read in files  -  seurat with all metadata #####
seurat <- CreateSeuratObject(counts, meta.data=assignments)
seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], genes$Gene_ID, col.name="Gene_ID")
seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], genes$Gene_ID, col.name="ENSG_ID")
saveRDS(seurat,paste0(args$out, "/seurat_object_all_pools_all_barcodes_all_metadata.rds"))


##### create seurat object  -  seurat with metadata #####
seurat <- CreateSeuratObject(counts, meta.data=assignments_final)
seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], genes$Gene_ID, col.name="Gene_ID")
seurat[["RNA"]] <- AddMetaData(seurat[["RNA"]], genes$Gene_ID, col.name="ENSG_ID")
saveRDS(seurat,paste0(args$out, "/seurat_object_all_pools_all_barcodes_final_assignments.rds"))


##### subset for just singlets #####
seurat <- subset(seurat, DropletType == "singlet")

##### Get the mitochondiral and ribosomal percentage QC metrics for each cell #####
if ((sum(which(MT_genes$GeneID %in% rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat)))) & (sum(which(MT_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse="|"), rownames(seurat))))){
    mt_features <- MT_genes$GeneID[MT_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features=mt_features)
} else if ((sum(which(MT_genes$ENSG %in% rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (sum(which(MT_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse="|"), rownames(seurat))))){
    mt_features <- MT_genes$ENSG[MT_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features=mt_features)
} else if ((length(grep(paste0(MT_genes$ENSG, "\\.", collapse="|"), rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(MT_genes$ENSG, "\\.", collapse="|"), rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features=rownames(seurat)[grep(paste0(MT_genes$ENSG, "\\.", collapse="|"), rownames(seurat))])
} else {
    message("Either you do not have mitochondrial genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

if ((sum(which(RB_genes$GeneID %in% rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat)))) & (sum(which(RB_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse="|"), rownames(seurat))))){
    rb_features <- RB_genes$GeneID[RB_genes$GeneID %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features=rb_features)
} else if ((sum(which(RB_genes$ENSG %in% rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (sum(which(RB_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse="|"), rownames(seurat))))){
    rb_features <- RB_genes$ENSG[RB_genes$ENSG %in% rownames(seurat)]
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features=rb_features)
} else if ((length(grep(paste0(RB_genes$ENSG, "\\.", collapse="|"), rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(RB_genes$ENSG, "\\.", collapse="|"), rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features=rownames(seurat)[grep(paste0(RB_genes$ENSG, "\\.", collapse="|"), rownames(seurat))])
} else {
    message("Either you do not have ribosomal genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

##### Save the object with singlets and QC metrics #####
saveRDS(seurat,paste0(args$out, "/seurat_object_all_pools_singlet_barcodes_final_assignments.rds"))

##### Create a dataframe of MADs #####
MAD_df_list <- lapply(pools_list, function(x){
    as.data.frame(matrix(ncol=2, nrow=length((which(c("percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA") %in% colnames(seurat@meta.data))))))
})
names(MAD_df_list) <- pools_list

MAD_df_list <- lapply(names(MAD_df_list), function(x){
    rownames(MAD_df_list[[x]]) <- colnames(seurat@meta.data)[which(colnames(seurat@meta.data) %in% c("percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))]
    colnames(MAD_df_list[[x]]) <- c("Median", "MAD")
    for (QC in rownames(MAD_df_list[[x]])){
        MAD_df_list[[x]][QC, "Median"] <- median(seurat@meta.data[which(seurat@meta.data$Pool == x),QC], na.rm=TRUE)
        MAD_df_list[[x]][QC, "MAD"] <- mad(seurat@meta.data[which(seurat@meta.data$Pool == x),QC], center=MAD_df_list[[x]][QC, "Median"],  constant=1.4826, na.rm=TRUE,low=FALSE, high=FALSE)
    }
    MAD_df_list[[x]]$Pool <- x
    MAD_df_list[[x]]$QC_Metric <- rownames(MAD_df_list[[x]])
    return(MAD_df_list[[x]])
})
MAD_df <- do.call(rbind,MAD_df_list)
MAD_df$Pool <- as.factor(MAD_df$Pool)

MAD_df_All <- as.data.frame(matrix(ncol=2, nrow=length((which(c("percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA") %in% colnames(seurat@meta.data))))))
rownames(MAD_df_All) <- colnames(seurat@meta.data)[which(colnames(seurat@meta.data) %in% c("percent.rb", "percent.mt", "nCount_RNA", "nFeature_RNA"))]
colnames(MAD_df_All) <- c("Median", "MAD")
for (QC in rownames(MAD_df_All)){
    MAD_df_All[QC, "Median"] <- median(seurat@meta.data[,QC], na.rm=TRUE)
    MAD_df_All[QC, "MAD"] <- mad(seurat@meta.data[,QC], center=MAD_df_All[QC, "Median"],  constant=1.4826, na.rm=TRUE,low=FALSE, high=FALSE)
}
MAD_df_All$QC_Metric <- rownames(MAD_df_All)



# MAD_df_long <- pivot_longer(MAD_df, names_to="QC_metric", cols=c("mt.percent", "nUMI", "nFeature"), values_to="value")

##### Make figures #####
### Mt % ###
## Add multiple MAD lines ##
violins_MADperPOOL <- list()
violins_MAD_ALL <- list()
violins <- list()

for (QC in unique(MAD_df$QC_Metric)){
    if (QC == "percent.mt" | QC == "percent.rb"){
        violins_MADperPOOL[[QC]] <- ggplot(seurat@meta.data, aes( x=Pool, y=seurat@meta.data[,QC])) + 
                    geom_violin() + 
                    geom_sina(size=1, alpha=0.6) + 
                    theme_classic() + 
                    ylim(0, min(max(seurat@meta.data[,QC]), 100)) + 
                    labs(title=paste0(QC, " Violin Plot"), y=QC, x="Pool") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median, yend=Median, col="Median"), size=1) + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + MAD, yend=Median + MAD,col="1 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - MAD, yend=Median - MAD, col="1 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 2 * MAD, yend=Median + 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 2 * MAD, yend=Median - 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 3 * MAD, yend=Median + 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 3 * MAD, yend=Median - 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 4 * MAD, yend=Median + 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 4 * MAD, yend=Median - 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 5 * MAD, yend=Median + 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 5 * MAD, yend=Median - 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                    scale_color_manual("MAD", values=c("Median"="grey", "1 MAD"="blue3", "2 MAD"="darkviolet", "3 MAD"="firebrick1", "4 MAD"="darkorange1", "5 MAD"="gold1")) + 
                    theme(text=element_text(size=14),
                        axis.text.x=element_text(angle=45, hjust=1),
                        plot.title=element_text(hjust=0.5),
                        plot.subtitle=element_text(hjust=0.5))

        violins_MAD_ALL[[QC]] <- ggplot(seurat@meta.data, aes( x=Pool, y=seurat@meta.data[,QC])) + 
                    geom_violin() + 
                    geom_sina(size=1, alpha=0.6) + 
                    theme_classic() + 
                    ylim(0, min(max(seurat@meta.data[,QC]), 100)) + 
                    labs(title=paste0(QC, " Violin Plot"), y=QC, x="Pool") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median, yend=Median, col="Median"), size=1) + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + MAD, yend=Median + MAD,col="1 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - MAD, yend=Median - MAD, col="1 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 2 * MAD, yend=Median + 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 2 * MAD, yend=Median - 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 3 * MAD, yend=Median + 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 3 * MAD, yend=Median - 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 4 * MAD, yend=Median + 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 4 * MAD, yend=Median - 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 5 * MAD, yend=Median + 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                    geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 5 * MAD, yend=Median - 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                    scale_color_manual("MAD", values=c("Median"="grey", "1 MAD"="blue3", "2 MAD"="darkviolet", "3 MAD"="firebrick1", "4 MAD"="darkorange1", "5 MAD"="gold1")) + 
                    theme(text=element_text(size=14),
                        axis.text.x=element_text(angle=45, hjust=1),
                        plot.title=element_text(hjust=0.5),
                        plot.subtitle=element_text(hjust=0.5))
    } else {
        violins_MADperPOOL[[QC]] <- ggplot(seurat@meta.data, aes( x=Pool, y=seurat@meta.data[,QC])) + 
                            geom_violin() + 
                            geom_sina(size=1, alpha=0.6) + 
                            theme_classic() + 
                            ylim(0, NA) + 
                            labs(title=paste0(QC, " Violin Plot"), y=QC, x="Pool") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median, yend=Median, col="Median"), size=1) + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + MAD, yend=Median + MAD,col="1 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - MAD, yend=Median - MAD, col="1 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 2 * MAD, yend=Median + 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 2 * MAD, yend=Median - 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 3 * MAD, yend=Median + 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 3 * MAD, yend=Median - 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 4 * MAD, yend=Median + 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 4 * MAD, yend=Median - 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median + 5 * MAD, yend=Median + 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC), ], aes(x=as.numeric(as.factor(Pool)) - 0.5, xend=as.numeric(as.factor(Pool)) + 0.5, y=Median - 5 * MAD, yend=Median - 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                            scale_color_manual("MAD", values=c("Median"="grey", "1 MAD"="blue3", "2 MAD"="darkviolet", "3 MAD"="firebrick1", "4 MAD"="darkorange1", "5 MAD"="gold1")) + 
                            theme(text=element_text(size=14),
                                axis.text.x=element_text(angle=45, hjust=1),
                                plot.title=element_text(hjust=0.5),
                                plot.subtitle=element_text(hjust=0.5))

        violins_MAD_ALL[[QC]] <- ggplot(seurat@meta.data, aes( x=Pool, y=seurat@meta.data[,QC])) + 
                            geom_violin() + 
                            geom_sina(size=1, alpha=0.6) + 
                            theme_classic() + 
                            ylim(0, NA) + 
                            labs(title=paste0(QC, " Violin Plot"), y=QC, x="Pool") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median, yend=Median, col="Median"), size=1) + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + MAD, yend=Median + MAD,col="1 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - MAD, yend=Median - MAD, col="1 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 2 * MAD, yend=Median + 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 2 * MAD, yend=Median - 2 * MAD, col="2 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 3 * MAD, yend=Median + 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 3 * MAD, yend=Median - 3 * MAD, col="3 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 4 * MAD, yend=Median + 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 4 * MAD, yend=Median - 4 * MAD, col="4 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median + 5 * MAD, yend=Median + 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                            geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC), ], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool))) - 0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool))) + 0.5, y=Median - 5 * MAD, yend=Median - 5 * MAD, col="5 MAD"), size=1, linetype="longdash") + 
                            scale_color_manual("MAD", values=c("Median"="grey", "1 MAD"="blue3", "2 MAD"="darkviolet", "3 MAD"="firebrick1", "4 MAD"="darkorange1", "5 MAD"="gold1")) + 
                            theme(text=element_text(size=14),
                                axis.text.x=element_text(angle=45, hjust=1),
                                plot.title=element_text(hjust=0.5),
                                plot.subtitle=element_text(hjust=0.5))
    }
    violins[[QC]] <- ggplot(seurat@meta.data, aes( x=Pool, y=seurat@meta.data[,QC])) + 
                    geom_violin() + 
                    geom_sina(size=1, alpha=0.6) + 
                    theme_classic() + 
                    ylim(0, NA) + 
                    labs(title=paste0(QC, " Violin Plot"), y=QC, x="Pool") + 
                    theme(text=element_text(size=14),
                        axis.text.x=element_text(angle=45, hjust=1),
                        plot.title=element_text(hjust=0.5),
                        plot.subtitle=element_text(hjust=0.5))
}
 
for (QC in names(violins_MAD_ALL)){
    ggsave(violins_MADperPOOL[[QC]], filename=paste0(args$out, "/", QC, "_violin_MADper_Pool.png"), width=29.7, height=21 ,units=c("cm"))
    ggsave(violins_MAD_ALL[[QC]], filename=paste0(args$out, "/", QC, "_violin_MAD_All.png"), width=29.7, height=21 ,units=c("cm"))
    ggsave(violins[[QC]], filename=paste0(args$out, "/", QC, "_violin_noMADlines.png"), width=29.7, height=21 ,units=c("cm"))
}

MAD_df_All$MAD1_up <- MAD_df_All$Median +  MAD_df_All$MAD
MAD_df_All$MAD1_down <- MAD_df_All$Median  -  MAD_df_All$MAD
MAD_df_All$MAD2_up <- MAD_df_All$Median +  2 * MAD_df_All$MAD
MAD_df_All$MAD2_down <- MAD_df_All$Median  -  2 * MAD_df_All$MAD
MAD_df_All$MAD3_up <- MAD_df_All$Median +  3 * MAD_df_All$MAD
MAD_df_All$MAD3_down <- MAD_df_All$Median  -  3 * MAD_df_All$MAD
MAD_df_All$MAD4_up <- MAD_df_All$Median +  4 * MAD_df_All$MAD
MAD_df_All$MAD4_down <- MAD_df_All$Median  -  4 * MAD_df_All$MAD
MAD_df_All$MAD5_up <- MAD_df_All$Median +  5 * MAD_df_All$MAD
MAD_df_All$MAD5_down <- MAD_df_All$Median  - 5 * MAD_df_All$MAD

MAD_df_All_long <- pivot_longer(MAD_df_All, cols=c("Median", "MAD1_up", "MAD1_down", "MAD2_up", "MAD2_down", "MAD3_up", "MAD3_down", "MAD4_up", "MAD4_down", "MAD5_up", "MAD5_down"), names_to="Metric")
MAD_df_All_long <- separate(MAD_df_All_long, col=Metric, sep="_", into=c("MAD_metric", "direction"))

colors <- c("grey", "blue3", "darkviolet", "firebrick1", "darkorange1", "gold1")
names(colors)=unique(MAD_df_All_long$MAD_metric)

### Number Genes vs Number UMI per Cell ###
## Color by Individual ##
## Add multiple MAD lines ##
if ("percent.mt" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric) & ("nFeature_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_MTscatter_MAD <- ggplot(seurat@meta.data, mapping=aes(x=nCount_RNA, y=percent.mt)) + 
        geom_point(size=0.5, alpha=0.5) + 
        theme_classic() + 
        labs(x="Number UMI", y="Percent Mitochondial Genes") +    
        xlim(0, NA) + 
        ylim(0, NA) + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nCount_RNA"), ], aes(x=value, xend=value, color=MAD_metric), y=min(seurat@meta.data$percent.mt), yend=max(seurat@meta.data$percent.mt), linetype="longdash") + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "percent.mt"), ], aes(y=value, yend=value, color=MAD_metric), x=min(seurat@meta.data$nCount_RNA), xend=max(seurat@meta.data$nCount_RNA), linetype="longdash") + 
        scale_color_manual(values=colors)
    ggsave(pUMI_MTscatter_MAD, filename=paste0(args$out, "/UMI_vs_percentMT_QC_scatter_w_MADlines.png"))

    pNfeatures_MTscatter_MAD <- ggplot(seurat@meta.data, mapping=aes(x=nFeature_RNA, y=percent.mt)) + 
        geom_point(size=0.5, alpha=0.5) + 
        theme_classic() + 
        labs(x="Number Features", y="Percent Mitochondial Genes") +    
        xlim(0, NA) + 
        ylim(0, NA) + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nFeature_RNA"), ], aes(x=value, xend=value, color=MAD_metric), y=min(seurat@meta.data$percent.mt), yend=max(seurat@meta.data$percent.mt), linetype="longdash") + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "percent.mt"), ], aes(y=value, yend=value, color=MAD_metric), x=min(seurat@meta.data$nFeature_RNA), xend=max(seurat@meta.data$nFeature_RNA), linetype="longdash") + 
        scale_color_manual(values=colors)
    ggsave(pNfeatures_MTscatter_MAD, filename=paste0(args$out, "/nFeatures_vs_percentMT_QC_scatter_w_MADlines.png"))

    pUMI_MTscatter <- ggplot(seurat@meta.data, aes(x=nCount_RNA, y=percent.mt, color=seurat@meta.data$Pool)) + 
        geom_point(alpha=0.5, size=0.5) + 
        theme_classic() + 
        scale_color_manual(values=rainbow(length(unique(seurat@meta.data$Pool))), name="Pool") + 
        labs(x="Number UMI", y="Percent Mitochondial Genes")     
    ggsave(pUMI_MTscatter, filename=paste0(args$out, "/UMI_vs_percentMT_QC_scatter_colorPool.png"))

    pNfeature_MTscatter <- ggplot(seurat@meta.data, aes(x=nFeature_RNA, y=percent.mt, color=seurat@meta.data$Pool)) + 
        geom_point(alpha=0.5, size=0.5) + 
        theme_classic() + 
        scale_color_manual(values=rainbow(length(unique(seurat@meta.data$Pool))), name="Pool") + 
        labs(x="Number Features", y="Percent Mitochondial Genes")     
    ggsave(pNfeature_MTscatter, filename=paste0(args$out, "/nFeatures_vs_percentMT_QC_scatter_colorPool.png"))

}

if ("nFeature_RNA" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_Genes_scatter_MAD <- ggplot() + 
        geom_point(data=seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA), alpha=0.25, size=0.5) + 
        theme_classic() + 
        labs(x="Number UMI", y="Number Genes")  + 
        xlim(0, NA) + 
        ylim(0, NA) + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nCount_RNA"), ], aes(x=value, xend=value, color=MAD_metric), y=min(seurat@meta.data$nFeature_RNA), yend=max(seurat@meta.data$nFeature_RNA), linetype="longdash") + 
        geom_segment(data=MAD_df_All_long[which(MAD_df_All_long$QC_Metric == "nFeature_RNA"), ], aes(y=value, yend=value, color=MAD_metric), x=min(seurat@meta.data$nCount_RNA), xend=max(seurat@meta.data$nCount_RNA), linetype="longdash") + 
        scale_color_manual(values=colors)
    ggsave(pUMI_Genes_scatter_MAD, filename=paste0(args$out, "/UMI_vs_Genes_QC_scatter_w_MADlines.png"))

    pUMI_Genes_scatter<- ggplot() + 
        geom_point(data=seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=as.factor(seurat@meta.data$Pool)), alpha=0.25, size=0.5) + 
        theme_classic() + 
        scale_color_manual(values=rainbow(length(unique(seurat@meta.data$Pool))), name="Pool") + 
        labs(x="Number UMI", y="Number Genes") 
        
    ggsave(pUMI_Genes_scatter, filename=paste0(args$out, "/UMI_vs_Genes_QC_scatter.png"))
}

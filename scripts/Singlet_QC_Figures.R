.libPaths("/usr/local/lib/R/site-library")
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Seurat)
# library(RColorBrewer)
# library(ggpubr)
# library(ggforce)
library(ggnewscale)

##### Set up arguemtns #####
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
dir <- arguments[1,]
pools_file <- arguments[2,]
datadir <- arguments[3,]
dirs_10x <- arguments[4,]
out <- arguments[5,]
RB_genes_file <- arguments[,8]
MT_genes_file <- arguments[,7]

pools <- read_delim(as.character(pools_file[1]), delim = "\t", col_names = c("Pool"))
pools_list <- pools$Pool

##### Read in the mt and rb gene lists #####
RB_genes <- read_delim(RB_genes_file, delim = "\t")
MT_genes <- read_delim(MT_genes_file, delim = "\t")

##### Read in and format the software assignments #####
assignments_list <- lapply(pools_list, function(x){
    read_delim(paste0(as.character(dir),"/",x,"/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"), delim = "\t")
})
names(assignments_list) <- pools_list

assignments_list <- lapply(pools_list, function(x){
    assignments_list[[x]]$Pool <- x
    assignments_list[[x]]$Barcode <- gsub("-1", "", assignments_list[[x]]$Barcode) 
    assignments_list[[x]]$Barcode  <- paste0(assignments_list[[x]]$Barcode, "_", x)
    return(assignments_list[[x]])
})
assignments <- do.call(rbind,assignments_list)
assignments <- as.data.frame(assignments)
rownames(assignments) <- assignments$Barcode

##### Read in and format the final assignments #####
assignments_final_list <- lapply(pools_list, function(x){
    read_delim(paste0(as.character(dir),"/",x,"/CombinedResults/Final_Assignments_demultiplexing_doublets.txt"), delim = "\t")
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
counts_list <- lapply(pools_list, function(x){
    Read10X(dirs_10x[grep(x, dirs_10x)], gene.column = 1)
})
names(counts_list) <- pools_list

counts_list <- lapply(pools_list, function(x){
    colnames(counts_list[[x]]) <- gsub("-1","", colnames(counts_list[[x]]))
    colnames(counts_list[[x]]) <- paste0(colnames(counts_list[[x]]),"_",x)
    return(counts_list[[x]])
})
counts <- do.call(cbind, counts_list)

##### Read in files - seurat with all metadata #####
seurat <- CreateSeuratObject(counts, meta.data = assignments)
saveRDS(seurat,paste0(out,"/seurat_object_all_pools_all_barcodes_all_metadata.rds"))


##### create seurat object - seurat with metadata #####
seurat <- CreateSeuratObject(counts, meta.data = assignments_final)
saveRDS(seurat,paste0(out,"/seurat_object_all_pools_all_barcodes_final_assignments.rds"))


##### subset for just singlets #####
seurat <- subset(seurat, subset = DropletType == "singlet")
saveRDS(seurat,paste0(out,"/seurat_object_all_pools_singlet_barcodes_final_assignments.rds"))

if ((sum(which(MT_genes$GeneID %in% rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat)))) & (sum(which(MT_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = MT_genes$GeneID)
} else if ((sum(which(MT_genes$ENSG %in% rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (sum(which(MT_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = MT_genes$ENSG)
} else if ((length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(MT_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(MT_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have mitochondrial genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

if ((sum(which(RB_genes$GeneID %in% rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat)))) & (sum(which(RB_genes$GeneID %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = RB_genes$GeneID)
} else if ((sum(which(RB_genes$ENSG %in% rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (sum(which(RB_genes$ENSG %in% rownames(seurat))) > length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = RB_genes$ENSG)
} else if ((length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$GeneID %in% rownames(seurat)))) & (length(grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))) > sum(which(RB_genes$ENSG %in% rownames(seurat))))){
    seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, features = rownames(seurat)[grep(paste0(RB_genes$ENSG, "\\.", collapse = "|"), rownames(seurat))])
} else {
    message("Either you do not have ribosomal genes in your dataset or they are not labeled with ENSG IDs or Gene IDs")
}

##### Create a dataframe of MADs #####
MAD_df_list <- lapply(pools_list, function(x){
    as.data.frame(matrix(ncol = 2, nrow = length((which(c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA") %in% colnames(seurat@meta.data))))))
})
names(MAD_df_list) <- pools_list

MAD_df_list <- lapply(names(MAD_df_list), function(x){
    rownames(MAD_df_list[[x]]) <- colnames(seurat@meta.data)[which(colnames(seurat@meta.data) %in% c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA"))]
    colnames(MAD_df_list[[x]]) <- c("Median","MAD")
    for (QC in rownames(MAD_df_list[[x]])){
        MAD_df_list[[x]][QC, "Median"] <- median(seurat@meta.data[which(seurat@meta.data$Pool == x),QC], na.rm = TRUE)
        MAD_df_list[[x]][QC, "MAD"] <- mad(seurat@meta.data[which(seurat@meta.data$Pool == x),QC], center = MAD_df_list[[x]][QC, "Median"],  constant = 1.4826, na.rm = TRUE,low = FALSE, high = FALSE)
    }
    MAD_df_list[[x]]$Pool <- x
    MAD_df_list[[x]]$QC_Metric <- rownames(MAD_df_list[[x]])
    return(MAD_df_list[[x]])
})
MAD_df <- do.call(rbind,MAD_df_list)
MAD_df$Pool <- as.factor(MAD_df$Pool)

MAD_df_All <- as.data.frame(matrix(ncol = 2, nrow = length((which(c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA") %in% colnames(seurat@meta.data))))))
rownames(MAD_df_All) <- colnames(seurat@meta.data)[which(colnames(seurat@meta.data) %in% c("percent.rb","percent.mt","nCount_RNA","nFeature_RNA"))]
colnames(MAD_df_All) <- c("Median","MAD")
for (QC in rownames(MAD_df_All)){
    MAD_df_All[QC, "Median"] <- median(seurat@meta.data[,QC], na.rm = TRUE)
    MAD_df_All[QC, "MAD"] <- mad(seurat@meta.data[,QC], center = MAD_df_All[QC, "Median"],  constant = 1.4826, na.rm = TRUE,low = FALSE, high = FALSE)
}
MAD_df_All$QC_Metric <- rownames(MAD_df_All)



# MAD_df_long <- pivot_longer(MAD_df, names_to = "QC_metric", cols = c("mt.percent","nUMI","nFeature"), values_to = "value")

##### Make figures #####
### Mt % ###
## Add multiple MAD lines ##
violins_MADperPOOL <- list()
violins_MAD_ALL <- list()
violins <- list()

for (QC in unique(MAD_df$QC_Metric)){
    violins_MADperPOOL[[QC]] <- ggplot(seurat@meta.data, aes( x = Pool, y = seurat@meta.data[,QC])) +
                        geom_violin() +
                        geom_sina(size = 1, alpha = 0.6) +
                        theme_classic() +
                        labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df[which(MAD_df$QC_Metric == QC),], aes(x=as.numeric(as.factor(Pool))-0.5, xend=as.numeric(as.factor(Pool))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                        scale_color_manual("MAD", values=c("Median"="grey","1 MAD"="blue3","2 MAD"="darkviolet", "3 MAD" = "firebrick1", "4 MAD" = "darkorange1", "5 MAD" = "gold1")) +
                        theme(text = element_text(size=14),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            plot.title = element_text(hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5))

    violins_MAD_ALL[[QC]] <- ggplot(seurat@meta.data, aes( x = Pool, y = seurat@meta.data[,QC])) +
                        geom_violin() +
                        geom_sina(size = 1, alpha = 0.6) +
                        theme_classic() +
                        labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median, yend=Median, col = "Median"), size = 1) +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+MAD, yend=Median+MAD,col = "1 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-MAD, yend=Median-MAD, col = "1 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+2*MAD, yend=Median+2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-2*MAD, yend=Median-2*MAD, col = "2 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+3*MAD, yend=Median+3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-3*MAD, yend=Median-3*MAD, col = "3 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+4*MAD, yend=Median+4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-4*MAD, yend=Median-4*MAD, col = "4 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median+5*MAD, yend=Median+5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == QC),], aes(x=min(as.numeric(as.factor(MAD_df_All$Pool)))-0.5, xend=max(as.numeric(as.factor(MAD_df_All$Pool)))+0.5, y=Median-5*MAD, yend=Median-5*MAD, col = "5 MAD"), size = 1, linetype = "longdash") +
                        scale_color_manual("MAD", values=c("Median"="grey","1 MAD"="blue3","2 MAD"="darkviolet", "3 MAD" = "firebrick1", "4 MAD" = "darkorange1", "5 MAD" = "gold1")) +
                        theme(text = element_text(size=14),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            plot.title = element_text(hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5))
    
    violins[[QC]] <- ggplot(seurat@meta.data, aes( x = Pool, y = seurat@meta.data[,QC])) +
                    geom_violin() +
                    geom_sina(size = 1, alpha = 0.6) +
                    theme_classic() +
                    labs(title = paste0(QC, " Violin Plot"), y = QC, x = "Pool") +
                    theme(text = element_text(size=14),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))
}
 
for (QC in names(violins_MAD_ALL)){
    ggsave(violins_MADperPOOL[[QC]], filename = paste0(out, "/", QC, "_violin_MADper_Pool.png"), width = 29.7, height = 21 ,units = c("cm"))
    ggsave(violins_MAD_ALL[[QC]], filename = paste0(out, "/", QC, "_violin_MAD_All.png"), width = 29.7, height = 21 ,units = c("cm"))
    ggsave(violins[[QC]], filename = paste0(out, "/", QC, "_violin_noMADlines.png"), width = 29.7, height = 21 ,units = c("cm"))
}


### Number Genes vs Number UMI per Cell ###
## Color by Individual ##
## Add multiple MAD lines ##
if ("percent.mt" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_MTscatter_MAD <- ggplot(seurat@meta.data, aes(x = nCount_RNA, y = percent.mt, color = Assignment)) +
        geom_scatter(alpha = 0.5) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(seurat@meta.data$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Percent Mitochondial Genes") +   
        xlim(0, NA) +
        ylim(0, NA) +
        new_scale_color() +     
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median, xend = Median, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "Median"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-1*MAD, xend = Median-1*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+1*MAD, xend = Median+1*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-2*MAD, xend = Median-2*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+2*MAD, xend = Median+2*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-3*MAD, xend = Median-3*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+3*MAD, xend = Median+3*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-4*MAD, xend = Median-4*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+4*MAD, xend = Median+4*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-5*MAD, xend = Median-5*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+5*MAD, xend = Median+5*MAD, y = min(seurat@meta.data$percent.mt), yend = max(seurat@meta.data$percent.mt), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median, yend = Median, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "Median"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median-1*MAD, yend = Median-1*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median+1*MAD, yend = Median+1*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median-2*MAD, yend = Median-2*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median+2*MAD, yend = Median+2*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median-3*MAD, yend = Median-3*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median+3*MAD, yend = Median+3*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median-4*MAD, yend = Median-4*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median+4*MAD, yend = Median+4*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median-5*MAD, yend = Median-5*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "percent.mt"),], aes(y=Median+5*MAD, yend = Median+5*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "5 MAD"), linetype = "longdash") +
        scale_color_manual("MAD", values=c("Median"="grey","1 MAD"="blue3","2 MAD"="darkviolet", "3 MAD" = "firebrick1", "4 MAD" = "darkorange1", "5 MAD" = "gold1")) 
    ggsave(pUMI_MTscatter_MAD, filename = paste0(out, "UMI_vs_percentMT_QC_scatter_w_MADlines.png"))

    pUMI_MTscatter <- ggplot(seurat@meta.data, aes(x = nCount_RNA, y = percent.mt, color = Assignment)) +
        geom_scatter(alpha = 0.5) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(seurat@meta.data$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Percent Mitochondial Genes")     
    ggsave(pUMI_MTscatter, filename = paste0(out, "UMI_vs_percentMT_QC_scatter.png"))

}

if ("nFeature_RNA" %in% unique(MAD_df_All$QC_Metric) & ("nCount_RNA" %in% MAD_df_All$QC_Metric)){
    pUMI_Genes_scatter_MAD <- ggplot() +
        geom_point(data = seurat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = as.factor(seurat@meta.data$Pool)), alpha = 0.5, size = 1) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(seurat@meta.data$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Number Genes")  +
        xlim(0, NA) +
        ylim(0, NA) +
        new_scale_color() +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median, xend = Median, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "Median"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-1*MAD, xend = Median-1*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+1*MAD, xend = Median+1*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-2*MAD, xend = Median-2*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+2*MAD, xend = Median+2*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-3*MAD, xend = Median-3*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+3*MAD, xend = Median+3*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-4*MAD, xend = Median-4*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+4*MAD, xend = Median+4*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median-5*MAD, xend = Median-5*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nCount_RNA"),], aes(x=Median+5*MAD, xend = Median+5*MAD, y = min(seurat@meta.data$nFeature_RNA), yend = max(seurat@meta.data$nFeature_RNA), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median, yend = Median, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "Median"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median-1*MAD, yend = Median-1*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median+1*MAD, yend = Median+1*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "1 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median-2*MAD, yend = Median-2*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median+2*MAD, yend = Median+2*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "2 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median-3*MAD, yend = Median-3*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median+3*MAD, yend = Median+3*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "3 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median-4*MAD, yend = Median-4*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median+4*MAD, yend = Median+4*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "4 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median-5*MAD, yend = Median-5*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "5 MAD"), linetype = "longdash") +
        geom_segment(data=MAD_df_All[which(MAD_df_All$QC_Metric == "nFeature_RNA"),], aes(y=Median+5*MAD, yend = Median+5*MAD, x = min(seurat@meta.data$nCount_RNA), xend = max(seurat@meta.data$nCount_RNA), col = "5 MAD"), linetype = "longdash") +
        scale_color_manual(name = "MAD",values=c("Median"="grey","1 MAD"="blue3","2 MAD"="darkviolet", "3 MAD" = "firebrick1", "4 MAD" = "darkorange1", "5 MAD" = "gold1")) 

    ggsave(pUMI_Genes_scatter_MAD, filename = paste0(out, "UMI_vs_Genes_QC_scatter_w_MADlines.png"))

    pUMI_Genes_scatter_MAD <- ggplot() +
        geom_point(data = seurat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = as.factor(seurat@meta.data$Pool)), alpha = 0.5, size = 1) +
        theme_classic() +
        scale_color_manual(values = rainbow(length(unique(seurat@meta.data$Pool))), name = "Pool")+
        labs(x = "Number UMI", y = "Number Genes") 
        
    ggsave(pUMI_Genes_scatter_MAD, filename = paste0(out, "UMI_vs_Genes_QC_scatter.png"))
}
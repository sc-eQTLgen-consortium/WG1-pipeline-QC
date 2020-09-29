#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")
library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library(ggpubr)
library(RColorBrewer)

##### Set directories #####
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
dir <- arguments[1,]
pool <- arguments[2,]

##### Read in Files #####
results <- read_delim(paste0(dir,"/",pool,"/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"), delim = "\t")

# Demuxlet may have not run for some cells, so fill the entries
results[is.na(results$demuxlet_nSNP), "demuxlet_nSNP"] <- 0
results[is.na(results$demuxlet_DropletType), "demuxlet_DropletType"] <- "unassigned"
results[is.na(results$demuxlet_Assignment), "demuxlet_Assignment"] <- "unassigned"
results[is.na(results$demuxlet_SingletLLK), "demuxlet_SingletLLK"] <- 0
results[is.na(results$demuxlet_DoulbetLLK), "demuxlet_DoulbetLLK"] <- 1
results[is.na(results$demuxlet_DiffLLK), "demuxlet_DiffLLK"] <- 0

message("Creating List of Softwares")
##### make a list of all softwares #####
softwares <- c("demuxlet","souporcell","DoubletDetection","scds","scrublet")
demultiplexing_combn <- t(combn(softwares, 5, simplify = TRUE)) %>% apply(. , 1 , paste , collapse = "_" )

##### Make a dataframe to provide assignments for intersection of doublets #####
intersection_doublet_demultiplex <- results[,c("Barcode")]


message("Creating Intersection and Union Assignments for each combination")
### Get a df of the droplet type ###
temp_DropletType <- results[,paste0(softwares,"_DropletType")]
temp_Assignment <- results[,paste0(c("demuxlet","souporcell"),"_Assignment")]

# Create column called "DropletType_temp"
# In this column, if number of singlet assignments per cell are all singlets, then label as singlet
# Otherwise, label as doublets
intersection_doublet_demultiplex <- intersection_doublet_demultiplex %>% mutate(DropletType_temp = if_else(rowSums(temp_DropletType == "singlet") == ncol(temp_DropletType), "singlet", "doublet"))

# For assignments, if DropletType_temp is a singlet, get temp_Assignment for that cell and get the first assignment.
# If all assignments agree, then set it to the first one
intersection_doublet_demultiplex[,"Assignment"] <- ifelse((intersection_doublet_demultiplex[,"DropletType_temp"] == "singlet" & apply(temp_Assignment, 1, function(y) all(y == y[1]))),
  pull(temp_Assignment,1), 
    ifelse(intersection_doublet_demultiplex[,"DropletType_temp"] == "doublet", "doublet","unassigned"))
for (row in 1:nrow(intersection_doublet_demultiplex[,"Assignment"])){
    if (intersection_doublet_demultiplex[row,"Assignment"] == "unassigned"){
        intersection_doublet_demultiplex[row,"DropletType"] <- "unassigned"
    } else {
        intersection_doublet_demultiplex[row,"DropletType"] <- intersection_doublet_demultiplex[row,"DropletType_temp"]
    }
}
intersection_doublet_demultiplex[,"DropletType_temp"] <- NULL

write_delim(intersection_doublet_demultiplex, paste0(dir,"/",pool,"/CombinedResults/Final_Assignments_demultiplexing_doublets.txt"), delim = "\t")

##### Make a figure of the number of cells per assignment per pool #####
pDropletType <- ggplot(intersection_doublet_demultiplex, aes(x = DropletType, fill = DropletType)) +
    geom_bar(position = "dodge", stat = "count") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    labs(title=pool, subtitle = "Number of Each Droplet Type") +
    theme(text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

pAssignment <- ggplot(intersection_doublet_demultiplex, aes(x = Assignment, fill = DropletType)) +
    geom_bar(position = "dodge", stat = "count") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    labs(subtitle = "Number of Each Droplet Assignment") +
    theme(text = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


pCombined <- ggarrange(pDropletType, pAssignment, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
ggexport(pCombined, filename = paste0(dir,"/",pool, "/CombinedResults/DropletType_Assignment_BarPlot.png"), width = 2000, height = 2000, res = 300)
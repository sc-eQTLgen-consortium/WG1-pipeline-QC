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
souporcell_genotype_key_dir <- arguments[3,]


##### Read in Files #####
results <- read_delim(paste0(dir,"/",pool,"/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"), delim = "\t")

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

intersection_doublet_demultiplex[,"DropletType"] <- ifelse(rowSums(temp_DropletType == "singlet") == ncol(temp_DropletType), "singlet", "doublet")
intersection_doublet_demultiplex[,"Assignment"] <- ifelse((intersection_doublet_demultiplex[,"DropletType"] == "singlet" & apply(temp_Assignment, 1, function(y) all(y == y[1]))),
  pull(temp_Assignment,1), 
    ifelse(intersection_doublet_demultiplex[,"DropletType"] == "doublet", "doublet","unassigned"))

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
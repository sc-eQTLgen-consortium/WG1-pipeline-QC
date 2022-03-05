.libPaths("/usr/local/lib/R/site-library")
library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(cowplot)



##### Set directories #####
args <- commandArgs(trailingOnly=TRUE)

sample_sheet <- args[1]
outdir <- args[2]
basedir <- args[3]



##### Read in sample sheet #####
sample_sheet_dt <- fread(sample_sheet, sep = "\t")



##### Read in the CombinedResults #####
results_list <- lapply(sample_sheet_dt$Pool, function(pool){
    fread(paste0(basedir,"/",pool, "/CombinedResults/Final_Assignments_demultiplexing_doublets.txt"))
})
names(results_list) <- sample_sheet_dt$Pool


##### Make datatable that has # individuals, # singlets, # doublets per pool expected and observed
expected_observed_individuals_list <- lapply(names(results_list), function(pool){
    data.table(Pool = rep(pool, 2),
                N_Individuals = c(sample_sheet_dt[Pool == pool]$N_Individuals, length(unique(results_list[[pool]][!(Assignment %in% c("doublet", "unassigned"))]$Assignment))),
                Expected_Observed = c("Expected", "Observed"))
})

                
expected_observed_classification_list <- lapply(names(results_list), function(pool){
    data.table(Pool = rep(pool, 4),
                N_Droplets = c(nrow(results_list[[pool]]) - nrow(results_list[[pool]])^2*0.008/1000, nrow(results_list[[pool]][DropletType == "singlet"]), nrow(results_list[[pool]])^2*0.008/1000, nrow(results_list[[pool]][DropletType == "doublet"])),
                Classification = c("singlet", "singlet", "doublet", "doublet"),
                Expected_Observed = rep(c("Expected", "Observed"), 2))
})


##### Combine results to a single datatable #####
expected_observed_individuals_dt <- rbindlist(expected_observed_individuals_list)
expected_observed_classification_dt <- rbindlist(expected_observed_classification_list)


##### Make plots #####
### 1. Number of individuals expected vs detected ###
plot1 <- ggplot(expected_observed_individuals_dt, aes(Expected_Observed, N_Individuals)) +
    geom_bar(stat = "identity", position="stack") +
    facet_grid(~ Pool)+
    theme_classic() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(angle = 90)) +
    ylab("Number of Individuals") 
    



### 2. Number of singlets and doublets explected vs detected ###
plot2 <- ggplot(expected_observed_classification_dt, aes(Expected_Observed, N_Droplets, fill = Classification)) +
    geom_bar(stat = "identity") +
    facet_grid(~ Pool) +
    theme_classic() +
    theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Number of Droplets") +
    scale_fill_manual(values = c("#85929E","#73C6B6"))


ggsave(plot_grid(plot1, plot2, rel_heights = c(2.5, 2), ncol = 1, nrow = 2, align = "v", axis = "lr"), filename = paste0(outdir, "/expected_observed_individuals_classifications.png"), width = 1.1 * length(unique(expected_observed_classification_dt$Pool)) + 1.1)




library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("ComplexHeatmap")
library(awtools)
library(circlize)
library(ggpubr)
library(viridis)
library(ggnewscale)
library(RColorBrewer)

set.seed(79)

##### Set directories #####
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/SimulatedPools/"
out <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/SimulatedPools/SimulatedOverlap/"
sample_list <- dir(path = dir, pattern = "Size")

result_files <- paste0(dir,sample_list,"/CombinedResults/CombinedDropletAssignments.tsv")
names(result_files) <- sample_list
meta <- read_delim("/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/data/ONEK1K/Pool_Name_N_Individuals.txt", delim = "\t")

colors <- c(ppalette,"#006d43")

# ##### Read in the Pool Metadata #####
# sim_pool_meta <- read_delim(paste0(dir,"../Round1Overlap/SimulatedPoolFirst100.tsv"), delim = "\t")
# sim_pool_meta <- sim_pool_meta[,c("Individuals","dir" )]
# sim_pool_meta_list <- list()
# for (pool in sim_pool_meta$dir){
#   sim_pool_meta_list[[pool]] <- sim_pool_meta[which(sim_pool_meta$dir == pool), "Individuals"]
# }

# sim_pool_meta_list <- lapply(sim_pool_meta_list, function(x){
#   separate(x, col = Individuals, into = paste0("Individual_",1:(str_count(x$Individuals[1], pattern = ",") + 1)), sep = ",") 
# })

# sim_pool_meta_list <- lapply(sim_pool_meta_list, function(x){
#   pivot_longer(x, cols = colnames(x), names_to = "Individual_Number", values_to = "Individual")
# })

# sim_pool_meta_list <- lapply(sim_pool_meta_list, function(x){
#   x$Individual_Number <- as.numeric(as.character(gsub("Individual_","",x$Individual_Number)))
#   return(x)
# })

# ##### Read in the cell_info barcode files #####
# simulated_barcodes <- lapply(sample_list, function(x){
#   read_delim(paste0(dir,x,"/cell_info.tsv"), delim = "\t")
# })
# simulated_barcodes <- lapply(simulated_barcodes, function(x){
#   x$DropletType <- ifelse(grepl("[A,C,G,T]+-[0-9]+D", x$CB_pool), "doublet", ifelse(grepl("[A,C,G,T]+-[0-9]+S", x$CB_pool), "doublet","singlet"))
#   x$DoubletType <- ifelse(grepl("[A,C,G,T]+-[0-9]+D", x$CB_pool), "doublet_2individuals", ifelse(grepl("[A,C,G,T]+-[0-9]+S", x$CB_pool), "doublet_1individual","singlet"))
#   return(x)
# })

# simulated_barcodes <- lapply(simulated_barcodes, function(x){
#   colnames(x) <- gsub("CB_pool","Barcode", colnames(x)) %>% gsub("Sample_id","Individual_Number",.)
#   return(x)
# })
# names(simulated_barcodes) <- sample_list

# ##### Join barcodes with individual identifiers #####
# simulated_barcodes_ref <- lapply(names(simulated_barcodes), function(x){
#   left_join(simulated_barcodes[[x]], sim_pool_meta_list[[x]], by = c("Individual_Number"))
# })
# names(simulated_barcodes_ref) <- names(simulated_barcodes)


# ##### Read in Files #####
# results_list <- lapply(result_files, function(x){
#   read_delim(x, delim = "\t")
# })
# names(results_list) <- sample_list

# ##### Read in the key tables for freemuxlet, scSplit and souporcell #####
# key <- readRDS(paste0(out,"PoolKeys.rds"))
# key <- lapply(key, function(x){
#   x$Cluster_ID <- gsub("CLUST","",x$Cluster_ID)
#   x$Software <- paste0(x$Software, "_Assignment")
#   return(x)
# })

# ##### Pivot the dataframes longer for joining#####
# results_list_long <- lapply(results_list, function(x){
#   pivot_longer(x,cols = c("demuxlet_Assignment","freemuxlet_Assignment","scSplit_Assignment","souporcell_Assignment","vireo_Assignment"), names_to = "Software")
# })

# ##### Left_join the common assignments to the dataframe #####
# results_list_long <- lapply(names(key), function(x){
#   temp <- left_join(results_list_long[[x]],key[[x]], by = c("value" = "Cluster_ID", "Software" = "Software"))
#   temp$Genotype_ID <- ifelse(temp$Software == "demuxlet_Assignment", temp$value, temp$Genotype_ID)
#   temp$Genotype_ID <- ifelse(temp$Software == "vireo_Assignment", temp$value, temp$Genotype_ID)
#   temp$Genotype_ID <- ifelse(temp$value == "doublet", "doublet", temp$Genotype_ID)
#   temp$Genotype_ID <- ifelse(is.na(temp$Genotype_ID), "unassigned", temp$Genotype_ID)
#   return(temp)
# })
# names(results_list_long) <- names(key)

# unsure_unassigned_numbers <- as.data.frame(matrix(nrow = length(results_list_long), ncol = 3))
# colnames(unsure_unassigned_numbers) <- c("Pool","Number_unassigned_unsure","Percent_unassigned_unsure")
# unsure_unassigned_numbers$Pool <- names(results_list_long)
# for (pool in names(results_list_long)){
#   unsure_unassigned_numbers[which(unsure_unassigned_numbers$Pool == pool),"Number_unassigned_unsure"] <- length(which(results_list_long[[pool]]$Genotype_ID == "unassigned")) + length(which(results_list_long[[pool]]$Genotype_ID == "unsure"))
#   unsure_unassigned_numbers[which(unsure_unassigned_numbers$Pool == pool),"Percent_unassigned_unsure"] <- 100*((length(which(results_list_long[[pool]]$Genotype_ID == "unassigned")) + length(which(results_list_long[[pool]]$Genotype_ID == "unsure")))/nrow(results_list_long[[pool]]))
# }
# print(max(unsure_unassigned_numbers$Percent_unassigned_unsure))

# ##### Pivot wider to get the counts of singlet, doublets... and intersection across softwares #####
# results_list_wide <- lapply(results_list_long, function(x){
#   x$value <- NULL
#   x$Correlation <- NULL
#   x <- pivot_wider(x,names_from = "Software", values_from = "Genotype_ID")
#   return(x)
# })

# message("Creating List of Softwares")
# ##### make a list of all softwares #####
# demultiplexing_list <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo")
# doublet_detection_list <- c("DoubletDetection","DoubletFinder","scds","scrublet")

# demultiplexing_deoublet_detection_list <- c(demultiplexing_list,doublet_detection_list)

# ones <- demultiplexing_deoublet_detection_list
# twos_df <- t(combn(demultiplexing_deoublet_detection_list, 2, simplify = TRUE))
# twos <- apply( twos_df , 1 , paste , collapse = "_" )
# threes_df <- t(combn(demultiplexing_deoublet_detection_list, 3, simplify = TRUE))
# threes <- apply( threes_df , 1 , paste , collapse = "_" )
# fours_df <- t(combn(demultiplexing_deoublet_detection_list, 4, simplify = TRUE))
# fours <- apply( fours_df , 1 , paste , collapse = "_" )
# fives_df <- t(combn(demultiplexing_deoublet_detection_list, 5, simplify = TRUE))
# fives <- apply( fives_df , 1 , paste , collapse = "_" )
# sixes_df <- t(combn(demultiplexing_deoublet_detection_list, 6, simplify = TRUE))
# sixes <- apply( sixes_df , 1 , paste , collapse = "_" )
# sevens_df <- t(combn(demultiplexing_deoublet_detection_list, 7, simplify = TRUE))
# sevens <- apply( sevens_df , 1 , paste , collapse = "_" )
# eights_df <-t(combn(demultiplexing_deoublet_detection_list, 8, simplify = TRUE))
# eights <- apply( eights_df , 1 , paste , collapse = "_" )
# nines_df <- t(combn(demultiplexing_deoublet_detection_list, 9, simplify = TRUE))
# nines <- apply( nines_df , 1 , paste , collapse = "_" )

# all_combn <- c(ones,twos,threes,fours,fives,sixes,sevens,eights,nines)

# ##### Make a loop to make a dataframe to provide assignments for intersection of doublets #####
# union_intersection_demultiplex <- lapply(results_list_wide,function(x){
#   x[,c("Barcode")]
# })

# union_intersection_doublet_demultiplex <- lapply(results_list_wide,function(x){
#   x[,c("Barcode")]
# })

# names <- names(union_intersection_demultiplex)

# doublet_only_softwares <- c()
# demultiplex_software <- c()
# doublet_only_softwares_half <- c()
# demultiplex_software_half <- c()

# message("Creating Intersection and Union Assignments for each combination")
# for (comparison in 1:length(all_combn)){
#   print(comparison)
#   combination <- all_combn[comparison]
#   softwares <- strsplit(all_combn[comparison],"_")[[1]]
#   ### Get a df of the droplet type ###
#   temp_DropletType <- lapply(results_list_wide, function(x){
#     x[,paste0(softwares,"_DropletType")]
#   })
#   ### Get a df of the assignment ###
#   if (any(softwares %in% demultiplexing_list)){
#   temp_Assignment <- lapply(results_list_wide, function(x){
#     x[,paste0(softwares[(softwares %in% demultiplexing_list)], "_Assignment")]
#   })
#   }
#   ### Add columns for union and intersection assignments if all softwares are demultiplexing softwares
#   if(all(softwares %in% doublet_detection_list)){
#     message("DoubletsOnly")
#     # union_intersection_doublet_demultiplex <- lapply(names(union_intersection_doublet_demultiplex), function(x){
#     #     union_intersection_doublet_demultiplex[[x]][,paste0(combination,"_DropletType_Intersection")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") == ncol(temp_DropletType[[x]]),"doublet",
#     #       ifelse(rowSums(temp_DropletType[[x]] == "singlet") == ncol(temp_DropletType[[x]]), "singlet","unassigned"))
#     #   union_intersection_doublet_demultiplex[[x]][,paste0(combination,"_DropletType_Union")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") > 0,"doublet",
#     #     ifelse(rowSums(temp_DropletType[[x]] == "singlet") == ncol(temp_DropletType[[x]]), "singlet","unassigned"))
#     #   if(length(softwares) > 2){
#     #     half <- ceiling(length(softwares)/2)
#     #     union_intersection_doublet_demultiplex[[x]][,paste0(combination,"_DropletType_Union")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") >= half ,"doublet",
#     #     ifelse(rowSums(temp_DropletType[[x]] == "singlet") >= half, "singlet","unassigned"))
#     #   }
#     #   return(union_intersection_doublet_demultiplex[[x]])
#     # })
#     # names(union_intersection_doublet_demultiplex) <- names
#     # doublet_only_softwares <- c(doublet_only_softwares,combination)
#   } else {
#     union_intersection_demultiplex <- lapply(names(union_intersection_demultiplex), function(x){
#       union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Intersection")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") == ncol(temp_DropletType[[x]]),"doublet",
#         ifelse(rowSums(temp_DropletType[[x]] == "singlet") == ncol(temp_DropletType[[x]]), "singlet","unassigned"))
#       union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") > 0,"doublet",
#         ifelse(rowSums(temp_DropletType[[x]] == "singlet") == ncol(temp_DropletType[[x]]), "singlet","unassigned"))
#       union_intersection_demultiplex[[x]][,paste0(combination,"_Assignment_Intersection")] <- ifelse((union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Intersection")] == "singlet" & 
#         apply(temp_Assignment[[x]], 1, function(y) all(y == y[1]))),pull(temp_Assignment[[x]],1), pull(union_intersection_demultiplex[[x]],paste0(combination,"_DropletType_Intersection")))
#       union_intersection_demultiplex[[x]][,paste0(combination,"_Assignment_Union")] <- ifelse((union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union")] == "singlet" & 
#         apply(temp_Assignment[[x]], 1, function(y) all(y == y[1]))), pull(temp_Assignment[[x]],1), pull(union_intersection_demultiplex[[x]],paste0(combination,"_DropletType_Union")))
#       if(length(softwares) > 2){
#         half <- ceiling(length(softwares)/2)
#         union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union_Half")] <- ifelse(rowSums(temp_DropletType[[x]] == "doublet") >= half ,"doublet",
#           ifelse(rowSums(temp_DropletType[[x]] == "singlet") >= half, "singlet","unassigned"))
#           if(length(softwares %in% demultiplexing_list) > 2){
#             half <- ceiling(length(which(softwares %in% demultiplexing_list))/2)
#             max_temp <- apply(temp_Assignment[[x]],1,function(y) names(which.max(table(y))))
#             union_intersection_demultiplex[[x]][,paste0(combination,"_Assignment_Union_Half")] <- ifelse(union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union_Half")] == "doublet", "doublet",
#               ifelse((union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union_Half")] == "unassigned"),"unassigned",
#                 ifelse((union_intersection_demultiplex[[x]][,paste0(combination,"_DropletType_Union_Half")] == "singlet" & 
#                   rowSums(temp_Assignment[[x]] == max_temp)>=half), max_temp, "unassigned")))
#       }
#       }
#       return(union_intersection_demultiplex[[x]])
#     })
#     if(length(softwares) > 2){
#      doublet_only_softwares_half <- c(doublet_only_softwares_half[!(doublet_only_softwares_half %in% demultiplex_software_half)],combination)
#     }
#     if(length(softwares %in% demultiplexing_list) > 2){
#       demultiplex_software_half <- c(demultiplex_software_half,combination)
#     }
#     names(union_intersection_demultiplex) <- names
#     demultiplex_software <- c(demultiplex_software,combination)
#   }
# }
# # saveRDS(union_intersection_doublet_demultiplex,paste0(out,"Intersection_Union_Assignments_doublets.rds"))
# saveRDS(union_intersection_demultiplex,paste0(out,"Intersection_Union_Assignments_demultiplexing.rds"))

# ##### Left join to add in the reference individual assignments #####
# message("Add reference information to dataframes")
# union_intersection_demultiplex <- lapply(names(union_intersection_demultiplex), function(x){
#   left_join(union_intersection_demultiplex[[x]], simulated_barcodes_ref[[x]], by = c("Barcode"))
# })

# union_intersection_demultiplex <- lapply(union_intersection_demultiplex, function(x){
#   x$Individual <- ifelse(x$DropletType == "singlet", x$Individual, x$DropletType)
#   return(x)
# })
# names(union_intersection_demultiplex) <- sample_list

# union_intersection_doublet_demultiplex <- lapply(names(union_intersection_doublet_demultiplex), function(x){
#   left_join(union_intersection_doublet_demultiplex[[x]], simulated_barcodes_ref[[x]], by = c("Barcode"))
# })

# union_intersection_doublet_demultiplex <- lapply(union_intersection_doublet_demultiplex, function(x){
#   x$Individual <- ifelse(x$DropletType == "singlet", x$Individual, x$DropletType)
#   return(x)
# })
# names(union_intersection_doublet_demultiplex) <- sample_list

# saveRDS(union_intersection_doublet_demultiplex,paste0(out,"Intersection_Union_Assignments_doublets.rds"))
# saveRDS(union_intersection_demultiplex,paste0(out,"Intersection_Union_Assignments_demultiplexing.rds"))

union_intersection_demultiplex <- readRDS(paste0(out,"Intersection_Union_Assignments_demultiplexing.rds"))
union_intersection_doublet_demultiplex <- readRDS(paste0(out,"Intersection_Union_Assignments_doublets.rds"))

##### Create dataframes for benchmarking results
message("Create")
# Benchmarking_demultiplex <- lapply(union_intersection_demultiplex, function(x){
#   df <- as.data.frame(matrix(nrow = 6, ncol = (2*length(demultiplex_software) + length(demultiplex_software_half) +1)))
#   rownames(df) <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy")
#   if (length(demultiplex_software_half) > 0){
#     colnames(df) <- c("Assessment",paste0(demultiplex_software,"_Intersection"),paste0(demultiplex_software,"_Union"), paste0(demultiplex_software_half,"_Union_Half"))
#   } else {
#     colnames(df) <- c("Assessment",paste0(demultiplex_software,"_Intersection"),paste0(demultiplex_software,"_Union"))
#   }
#   df$Assessment <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy")
#   return(df)
# })

# Benchmarking_doublet <- lapply(union_intersection_doublet_demultiplex, function(x){
#   df <- as.data.frame(matrix(nrow = 6, ncol = (2*length(doublet_only_softwares) + length(doublet_only_softwares_half) +1)))
#   rownames(df) <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy")
#   if (length(doublet_only_softwares_half) > 0 ){
#     colnames(df) <- c("Assessment",paste0(doublet_only_softwares,"_Intersection"),paste0(doublet_only_softwares,"_Union"), paste0(doublet_only_softwares_half,"_Union_Half"))
#   } else {
#     colnames(df) <- c("Assessment",paste0(doublet_only_softwares,"_Intersection"),paste0(doublet_only_softwares,"_Union"))
#   }
#   df$Assessment <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy")
#   return(df)
# })

## Try new method to get colnames ##
Benchmarking_demultiplex <- lapply(union_intersection_demultiplex, function(x){
  df <- as.data.frame(matrix(nrow = 11, ncol = 1))
  rownames(df) <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy","TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total")
  colnames(df) <- c("Assessment")
  df$Assessment <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy","TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total")
  return(df)
})

demultiplex_combinations_list <- unique(gsub("_DropletType_Intersection","",grep("_DropletType", value = TRUE, colnames(union_intersection_demultiplex[[1]]))) %>% gsub("_DropletType_Union_Half", "",.) %>% gsub("_DropletType_Union","", .))

Benchmarking_doublet <- lapply(union_intersection_doublet_demultiplex, function(x){
  df <- as.data.frame(matrix(nrow = 11, ncol = 1))
  rownames(df) <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy","TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total")
  colnames(df) <- c("Assessment")
  df$Assessment <- c("TrueSingletRate","FalseSingletRate","TrueDoubletRate","FalseDoubletRate","MCC","Accuracy","TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total")
  return(df)
})

doublet_combinations_list <- unique(gsub("_DropletType_Intersection","",grep("_DropletType", value = TRUE, colnames(union_intersection_doublet_demultiplex[[1]]))) %>% gsub("_DropletType_Union_Half", "",.) %>% gsub("_DropletType_Union","", .))

##### Calculate TrueSingletRate, FalseSingletRate, TrueDoubletRate and FalseDoubletRate #####
message("Calculating quality measures")
Benchmarking_demultiplex <- lapply(names(Benchmarking_demultiplex), function(x){
  for (y in demultiplex_combinations_list){
    print(y)
    ### For the intersection method ###
    TP <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "singlet" &  union_intersection_demultiplex[[x]][,paste0(y,"_Assignment_Intersection")] == union_intersection_demultiplex[[x]][,"Individual"] & union_intersection_demultiplex[[x]][,"DropletType"] == "singlet"))
    TN <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "doublet" & union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == union_intersection_demultiplex[[x]][,"DropletType"] & union_intersection_demultiplex[[x]][,"DropletType"] == "doublet"))
    FN <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "doublet")) - TN

    Benchmarking_demultiplex[[x]]["TrueSinglet",paste0(y,"_Intersection")] <- TP
    Benchmarking_demultiplex[[x]]["TrueDoublet",paste0(y,"_Intersection")] <- TN
    Benchmarking_demultiplex[[x]]["FalseSinglet",paste0(y,"_Intersection")] <- FP
    Benchmarking_demultiplex[[x]]["FalseDoublet",paste0(y,"_Intersection")] <- FN
    Benchmarking_demultiplex[[x]]["Total",paste0(y,"_Intersection")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_Intersection")] <- TP/(TP+FN) ## TPR
    Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_Intersection")] <- TN/(TN + FP) ## TNR
    Benchmarking_demultiplex[[x]]["FalseSingletRate",paste0(y,"_Intersection")] <- 1 - Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_Intersection")] ## FPR
    Benchmarking_demultiplex[[x]]["FalseDoubletRate",paste0(y,"_Intersection")] <- 1 - Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_Intersection")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    Benchmarking_demultiplex[[x]]["MCC",paste0(y,"_Intersection")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    Benchmarking_demultiplex[[x]]["Accuracy",paste0(y,"_Intersection")] <- (TP+TN)/(TP+TN+FP+FN)

    ### For the doublet union method ###
    TP <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == "singlet" &  union_intersection_demultiplex[[x]][,paste0(y,"_Assignment_Union")] == union_intersection_demultiplex[[x]][,"Individual"] & union_intersection_demultiplex[[x]][,"DropletType"] == "singlet"))
    TN <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == "doublet" & union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == union_intersection_demultiplex[[x]][,"DropletType"] & union_intersection_demultiplex[[x]][,"DropletType"] == "doublet"))
    FN <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "doublet" )) - TN

    Benchmarking_demultiplex[[x]]["TrueSinglet",paste0(y,"_DoubletUnion")] <- TP
    Benchmarking_demultiplex[[x]]["TrueDoublet",paste0(y,"_DoubletUnion")] <- TN
    Benchmarking_demultiplex[[x]]["FalseSinglet",paste0(y,"_DoubletUnion")] <- FP
    Benchmarking_demultiplex[[x]]["FalseDoublet",paste0(y,"_DoubletUnion")] <- FN
    Benchmarking_demultiplex[[x]]["Total",paste0(y,"_DoubletUnion")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_DoubletUnion")] <- TP/(TP+FN) ## TPR
    Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_DoubletUnion")] <- TN/(TN + FP) ## TNR
    Benchmarking_demultiplex[[x]]["FalseSingletRate",paste0(y,"_DoubletUnion")] <- 1 - Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_DoubletUnion")] ## FPR
    Benchmarking_demultiplex[[x]]["FalseDoubletRate",paste0(y,"_DoubletUnion")] <- 1 - Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_DoubletUnion")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    Benchmarking_demultiplex[[x]]["MCC",paste0(y,"_DoubletUnion")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    Benchmarking_demultiplex[[x]]["Accuracy",paste0(y,"_DoubletUnion")] <- (TP+TN)/(TP+TN+FP+FN)

    ### For the doublet remainder method ###
    TP <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "singlet" &  union_intersection_demultiplex[[x]][,paste0(y,"_Assignment_Intersection")] == union_intersection_demultiplex[[x]][,"Individual"] & union_intersection_demultiplex[[x]][,"DropletType"] == "singlet"))
    TN <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "doublet" & union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] != "singlet" & union_intersection_demultiplex[[x]][,"DropletType"] == "doublet"))
    FN <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "doublet")) - TN

    Benchmarking_demultiplex[[x]]["TrueSinglet",paste0(y,"_DoubletRemainder")] <- TP
    Benchmarking_demultiplex[[x]]["TrueDoublet",paste0(y,"_DoubletRemainder")] <- TN
    Benchmarking_demultiplex[[x]]["FalseSinglet",paste0(y,"_DoubletRemainder")] <- FP
    Benchmarking_demultiplex[[x]]["FalseDoublet",paste0(y,"_DoubletRemainder")] <- FN
    Benchmarking_demultiplex[[x]]["Total",paste0(y,"_DoubletRemainder")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_DoubletRemainder")] <- TP/(TP+FN) ## TPR
    Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_DoubletRemainder")] <- TN/(TN + FP) ## TNR
    Benchmarking_demultiplex[[x]]["FalseSingletRate",paste0(y,"_DoubletRemainder")] <- 1 - Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_DoubletRemainder")] ## FPR
    Benchmarking_demultiplex[[x]]["FalseDoubletRate",paste0(y,"_DoubletRemainder")] <- 1 - Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_DoubletRemainder")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    Benchmarking_demultiplex[[x]]["MCC",paste0(y,"_DoubletRemainder")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    Benchmarking_demultiplex[[x]]["Accuracy",paste0(y,"_DoubletRemainder")] <- (TP+TN)/(TP+TN+FP+FN)

    ### For the half union/intersection method ###
    if (y %in% gsub("_DropletType_Union_Half","",grep("_Union_Half",colnames(union_intersection_demultiplex[[x]]), value = TRUE))){
      print("HalfIntersection")
      TP <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union_Half")] == "singlet" & union_intersection_demultiplex[[x]][,paste0(y,"_Assignment_Union_Half")] == union_intersection_demultiplex[[x]][,"Individual"] & union_intersection_demultiplex[[x]][,"DropletType"] == "singlet"))
      TN <- length(which(union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union_Half")] == "doublet" & union_intersection_demultiplex[[x]][,paste0(y,"_DropletType_Union_Half")] == union_intersection_demultiplex[[x]][,"DropletType"] & union_intersection_demultiplex[[x]][,"DropletType"] != "singlet"))
      FN <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
      FP <- length(which(union_intersection_demultiplex[[x]][,"DropletType"] == "doublet" )) - TN

    Benchmarking_demultiplex[[x]]["TrueSinglet",paste0(y,"_HalfIntersection")] <- TP
    Benchmarking_demultiplex[[x]]["TrueDoublet",paste0(y,"_HalfIntersection")] <- TN
    Benchmarking_demultiplex[[x]]["FalseSinglet",paste0(y,"_HalfIntersection")] <- FP
    Benchmarking_demultiplex[[x]]["FalseDoublet",paste0(y,"_HalfIntersection")] <- FN
    Benchmarking_demultiplex[[x]]["Total",paste0(y,"_HalfIntersection")] <- nrow(union_intersection_demultiplex[[x]])

      Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_HalfIntersection")] <- TP/(TP+FN) ## TPR
      Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_HalfIntersection")] <- TN/(TN + FP) ## TNR
      Benchmarking_demultiplex[[x]]["FalseSingletRate",paste0(y,"_HalfIntersection")] <- 1 - Benchmarking_demultiplex[[x]]["TrueDoubletRate",paste0(y,"_HalfIntersection")] ## FPR
      Benchmarking_demultiplex[[x]]["FalseDoubletRate",paste0(y,"_HalfIntersection")] <- 1 - Benchmarking_demultiplex[[x]]["TrueSingletRate",paste0(y,"_HalfIntersection")] ## FNR
      PPV <- TP/(TP+FP)
      TPR <- TP/(TP+FN)
      TNR <- TN/(TN+FN)
      NPV <- TN/(TN+FN)
      FDR <- FP/(FP+TP)
      FNR <- FN/(FN+TP)
      FPR <- FP/(FP+TN)
      FOR <- FN/(FN+TN)
      Benchmarking_demultiplex[[x]]["MCC",paste0(y,"_HalfIntersection")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
      Benchmarking_demultiplex[[x]]["Accuracy",paste0(y,"_HalfIntersection")] <- (TP+TN)/(TP+TN+FP+FN)
    }
  }
  Benchmarking_demultiplex[[x]]$Pool <- x
  return(Benchmarking_demultiplex[[x]])
})


Benchmarking_doublet <- lapply(names(Benchmarking_doublet), function(x){
  for (y in doublet_combinations_list){
    print(y)
    ### For the intersection method ###
    TP <- length(which(union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "singlet" &  union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == union_intersection_doublet_demultiplex[[x]][,"DropletType"]))
    TN <- length(which(union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "doublet" & union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == union_intersection_doublet_demultiplex[[x]][,"DropletType"]))
    FN <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "doublet")) - TN

    Benchmarking_doublet[[x]]["TrueSinglet",paste0(y,"_Intersection")] <- TP
    Benchmarking_doublet[[x]]["TrueDoublet",paste0(y,"_Intersection")] <- TN
    Benchmarking_doublet[[x]]["FalseSinglet",paste0(y,"_Intersection")] <- FP
    Benchmarking_doublet[[x]]["FalseDoublet",paste0(y,"_Intersection")] <- FN
    Benchmarking_doublet[[x]]["Total",paste0(y,"_Intersection")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_Intersection")] <- TP/(TP+FN) ## TPR
    Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_Intersection")] <- TN/(TN + FP) ## TNR
    Benchmarking_doublet[[x]]["FalseSingletRate",paste0(y,"_Intersection")] <- 1 - Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_Intersection")] ## FPR
    Benchmarking_doublet[[x]]["FalseDoubletRate",paste0(y,"_Intersection")] <- 1 - Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_Intersection")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    if (PPV > 0 & TPR > 0 & TNR > 0 & NPV > 0 & FDR > 0 & FNR > 0 & FPR > 0 & FOR > 0) {
      Benchmarking_doublet[[x]]["MCC",paste0(y,"_Intersection")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    } else {
      Benchmarking_doublet[[x]]["MCC",paste0(y,"_Intersection")] <- NA
    }
    Benchmarking_doublet[[x]]["Accuracy",paste0(y,"_Intersection")] <- (TP+TN)/(TP+TN+FP+FN)

    ### For the union method ###
    TP <- length(which(union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == "singlet" &  union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == union_intersection_doublet_demultiplex[[x]][,"DropletType"]))
    TN <- length(which(union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == "doublet" & union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Union")] == union_intersection_doublet_demultiplex[[x]][,"DropletType"]))
    FN <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "doublet" )) - TN

    Benchmarking_doublet[[x]]["TrueSinglet",paste0(y,"_DoubletUnion")] <- TP
    Benchmarking_doublet[[x]]["TrueDoublet",paste0(y,"_DoubletUnion")] <- TN
    Benchmarking_doublet[[x]]["FalseSinglet",paste0(y,"_DoubletUnion")] <- FP
    Benchmarking_doublet[[x]]["FalseDoublet",paste0(y,"_DoubletUnion")] <- FN
    Benchmarking_doublet[[x]]["Total",paste0(y,"_DoubletUnion")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_DoubletUnion")] <- TP/(TP+FN) ## TPR
    Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_DoubletUnion")] <- TN/(TN + FP) ## TNR
    Benchmarking_doublet[[x]]["FalseSingletRate",paste0(y,"_DoubletUnion")] <- 1 - Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_DoubletUnion")] ## FPR
    Benchmarking_doublet[[x]]["FalseDoubletRate",paste0(y,"_DoubletUnion")] <- 1 - Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_DoubletUnion")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    Benchmarking_doublet[[x]]["MCC",paste0(y,"_DoubletUnion")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    Benchmarking_doublet[[x]]["Accuracy",paste0(y,"_DoubletUnion")] <- (TP+TN)/(TP+TN+FP+FN)

    ### For the doublet remainder method ###
    TP <- length(which(union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == "singlet" &  union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] == union_intersection_doublet_demultiplex[[x]][,"DropletType"]))
    TN <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "doublet" & union_intersection_doublet_demultiplex[[x]][,paste0(y,"_DropletType_Intersection")] != "singlet"))
    FN <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "singlet")) - TP
    FP <- length(which(union_intersection_doublet_demultiplex[[x]][,"DropletType"] == "doublet")) - TN

    Benchmarking_doublet[[x]]["TrueSinglet",paste0(y,"_DoubletRemainder")] <- TP
    Benchmarking_doublet[[x]]["TrueDoublet",paste0(y,"_DoubletRemainder")] <- TN
    Benchmarking_doublet[[x]]["FalseSinglet",paste0(y,"_DoubletRemainder")] <- FP
    Benchmarking_doublet[[x]]["FalseDoublet",paste0(y,"_DoubletRemainder")] <- FN
    Benchmarking_doublet[[x]]["Total",paste0(y,"_DoubletRemainder")] <- nrow(union_intersection_demultiplex[[x]])

    Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_DoubletRemainder")] <- TP/(TP+FN) ## TPR
    Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_DoubletRemainder")] <- TN/(TN + FP) ## TNR
    Benchmarking_doublet[[x]]["FalseSingletRate",paste0(y,"_DoubletRemainder")] <- 1 - Benchmarking_doublet[[x]]["TrueDoubletRate",paste0(y,"_Intersection")] ## FPR
    Benchmarking_doublet[[x]]["FalseDoubletRate",paste0(y,"_DoubletRemainder")] <- 1 - Benchmarking_doublet[[x]]["TrueSingletRate",paste0(y,"_Intersection")] ## FNR
    PPV <- TP/(TP+FP)
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FN)
    NPV <- TN/(TN+FN)
    FDR <- FP/(FP+TP)
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FOR <- FN/(FN+TN)
    Benchmarking_doublet[[x]]["MCC",paste0(y,"_DoubletRemainder")] <- sqrt(PPV*TPR*TNR*NPV) - sqrt(FDR*FNR*FPR*FOR)
    Benchmarking_doublet[[x]]["Accuracy",paste0(y,"_DoubletRemainder")] <- (TP+TN)/(TP+TN+FP+FN)
  }
  Benchmarking_doublet[[x]]$Pool <- x
  return(Benchmarking_doublet[[x]])
})

saveRDS(Benchmarking_doublet, paste0(out,"doublet_simulated_assessment_results.rds"))
saveRDS(Benchmarking_demultiplex, paste0(out,"demultiplex_simulated_assessment_results.rds"))
message("Combining dataframes")

Benchmarking_doublet <- readRDS(paste0(out,"doublet_simulated_assessment_results.rds"))
Benchmarking_demultiplex <- readRDS(paste0(out,"demultiplex_simulated_assessment_results.rds"))

##### Separate into two dataframes - one with the metrics and one with the absolute numbers #####
Benchmarking_doublet_numbers <- lapply(Benchmarking_doublet, function(x){
  x[c("TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total"),]
})

Benchmarking_doublet <- lapply(Benchmarking_doublet, function(x){
  x[which(!c(rownames(x) %in% c("TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total"))),]
})

Benchmarking_demultiplex_numbers <- lapply(Benchmarking_demultiplex, function(x){
  x[c("TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total"),]
})

Benchmarking_demultiplex <- lapply(Benchmarking_demultiplex, function(x){
  x[which(!c(rownames(x) %in% c("TrueSinglet","TrueDoublet","FalseSinglet","FalseDoublet","Total"))),]
})

##### Pivot longer #####
Benchmarking_demultiplex <- lapply(Benchmarking_demultiplex, function(x){
  pivot_longer(x,-c(Assessment,Pool), names_to="Comparison",values_to="Value")
})

Benchmarking_demultiplex_numbers <- lapply(Benchmarking_demultiplex_numbers, function(x){
  pivot_longer(x,-c(Assessment,Pool), names_to="Comparison",values_to="Value")
})

Benchmarking_doublet <- lapply(Benchmarking_doublet, function(x){
  pivot_longer(x,-c(Assessment,Pool), names_to="Comparison",values_to="Value")
})

Benchmarking_doublet_numbers <- lapply(Benchmarking_doublet_numbers, function(x){
  pivot_longer(x,-c(Assessment,Pool), names_to="Comparison",values_to="Value")
})

Benchmarking_demultiplex_df <- do.call(rbind, Benchmarking_demultiplex)
Benchmarking_demultiplex_numbers_df <- do.call(rbind, Benchmarking_demultiplex_numbers)
Benchmarking_doublet_df <- do.call(rbind, Benchmarking_doublet)
Benchmarking_doublet_numbers_df <- do.call(rbind, Benchmarking_doublet_numbers)
Benchmarking <- rbind(Benchmarking_demultiplex_df,Benchmarking_doublet_df)
Benchmarking_counts <- rbind(Benchmarking_demultiplex_numbers_df, Benchmarking_doublet_numbers_df)

message("Writing output")
write_delim(Benchmarking, delim = "\t", path = paste0(out,"Benchmarking_simulated_results.tsv"))
write_delim(Benchmarking_counts, delim = "\t", path = paste0(out,"Benchmarking_simulated_results_counts.tsv"))
Benchmarking <- read_delim(delim = "\t", paste0(out,"Benchmarking_simulated_results.tsv"))
Benchmarking_counts <- read_delim(delim = "\t", paste0(out,"Benchmarking_simulated_results_counts.tsv"))

Benchmarking_wide <- pivot_wider(Benchmarking, names_from="Assessment", values_from = "Value")
Benchmarking_counts_wide <- pivot_wider(Benchmarking_counts, names_from="Assessment", values_from = "Value")
Benchmarking_wide <- Benchmarking_wide[which(Benchmarking_wide$Comparison != "demuxlet_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "freemuxlet_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "scSplit_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "souporcell_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "vireo_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "DoubletDetection_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "DoubletFinder_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "scds_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "scrublet_DoubletUnion" &
                                              Benchmarking_wide$Comparison != "demuxlet_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "freemuxlet_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "scSplit_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "souporcell_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "vireo_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "DoubletDetection_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "DoubletFinder_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "scds_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "scrublet_DoubletRemainder" &
                                              Benchmarking_wide$Comparison != "demuxlet_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "freemuxlet_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "scSplit_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "souporcell_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "vireo_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "DoubletDetection_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "DoubletFinder_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "scds_HalfIntersection" &
                                              Benchmarking_wide$Comparison != "scrublet_HalfIntersection" ),]
Benchmarking_counts_wide <- Benchmarking_counts_wide[which(Benchmarking_counts_wide$Comparison != "demuxlet_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "freemuxlet_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "scSplit_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "souporcell_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "vireo_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "DoubletDetection_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "DoubletFinder_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "scds_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "scrublet_DoubletUnion" &
                                              Benchmarking_counts_wide$Comparison != "demuxlet_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "freemuxlet_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "scSplit_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "souporcell_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "vireo_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "DoubletDetection_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "DoubletFinder_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "scds_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "scrublet_DoubletRemainder" &
                                              Benchmarking_counts_wide$Comparison != "demuxlet_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "freemuxlet_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "scSplit_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "souporcell_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "vireo_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "DoubletDetection_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "DoubletFinder_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "scds_HalfIntersection" &
                                              Benchmarking_counts_wide$Comparison != "scrublet_HalfIntersection" ),]
Benchmarking_wide$N_softwares <- (str_count(gsub("scrublet_DoubletUnion","", Benchmarking_wide$Comparison) %>% gsub("_Intersection","",.) %>% gsub("_HalfIntersection","", .) %>% gsub("_DoubletRemainder","", .), pattern = "_") + 1)
Benchmarking_counts_wide$N_softwares <- (str_count(gsub("scrublet_DoubletUnion","", Benchmarking_counts_wide$Comparison) %>% gsub("_Intersection","",.) %>% gsub("_HalfIntersection","", .) %>% gsub("_DoubletRemainder","", .), pattern = "_") + 1)
Benchmarking_wide$Combination_Type <- gsub("demuxlet_","",Benchmarking_wide$Comparison) %>%
                                        gsub("freemuxlet_","",.) %>%
                                        gsub("scSplit_","",.) %>%
                                        gsub("souporcell_","",.) %>%
                                        gsub("vireo_","",.) %>%
                                        gsub("DoubletDetection_","",.) %>%
                                        gsub("DoubletFinder_","",.) %>%
                                        gsub("scds_","",.) %>%
                                        gsub("scrublet_","",.)
Benchmarking_counts_wide$Combination_Type <- gsub("demuxlet_","",Benchmarking_counts_wide$Comparison) %>%
                                        gsub("freemuxlet_","",.) %>%
                                        gsub("scSplit_","",.) %>%
                                        gsub("souporcell_","",.) %>%
                                        gsub("vireo_","",.) %>%
                                        gsub("DoubletDetection_","",.) %>%
                                        gsub("DoubletFinder_","",.) %>%
                                        gsub("scds_","",.) %>%
                                        gsub("scrublet_","",.)
Benchmarking_wide$N <- gsub("Size","",Benchmarking_wide$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
Benchmarking_counts_wide$N <- gsub("Size","",Benchmarking_counts_wide$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
Benchmarking_wide$N <- as.numeric(as.character(Benchmarking_wide$N))
Benchmarking_counts_wide$N <- as.numeric(as.character(Benchmarking_counts_wide$N))
Benchmarking_wide$Combination_Type <- factor(Benchmarking_wide$Combination_Type, levels = c("Intersection","DoubletUnion","DoubletRemainder","HalfIntersection"))
Benchmarking_counts_wide$Combination_Type <- factor(Benchmarking_counts_wide$Combination_Type, levels = c("Intersection","DoubletUnion","DoubletRemainder","HalfIntersection"))
write_delim(Benchmarking_wide, paste0(out, "Benchmarking_simulated_results_wide.tsv"), delim = "\t")
write_delim(Benchmarking_counts_wide, paste0(out, "Benchmarking_simulated_results_wide_counts.tsv"), delim = "\t")

Benchmarking_wide <- read_delim(paste0(out, "Benchmarking_simulated_results_wide.tsv"), delim = "\t")
Benchmarking_counts_wide <- read_delim(paste0(out, "Benchmarking_simulated_results_wide_counts.tsv"), delim = "\t")
# Benchmarking_wide$Combination_Groups <- NA
# for (row in 1:nrow(Benchmarking_wide)){
#   if (all(strsplit(gsub("_Union","", Benchmarking_wide$Comparison) %>% gsub("_Intersection","",.) %>% gsub("_Union_Half","", .),"_")[[row]] %in% demultiplexing_list)){
#     Benchmarking_wide$Combination_Groups[row] <- "All Demultiplexing"
#   }else if (all(strsplit(gsub("_Union","", Benchmarking_wide$Comparison) %>% gsub("_Intersection","",.) %>% gsub("_Union_Half","", .),"_")[[row]] %in% doublet_detection_list)){
#      Benchmarking_wide$Combination_Groups[row] <- "All Doublet Detection"
#   } else{
#     Benchmarking_wide$Combination_Groups[row] <- "Combination Demultiplexing and Doublet Detection"
#   }
# }
Benchmarking_wide$TrueDoubletRate_FalseSingletRate <- Benchmarking_wide$TrueDoubletRate/Benchmarking_wide$FalseSingletRate


p_TrueDoubletRate_FalseSingletRate_MCC <- ggplot(Benchmarking_wide, aes(x=TrueDoubletRate_FalseSingletRate,y=MCC,color=as.factor(N_softwares))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~Combination_Type, nrow = 1) +
  scale_color_brewer(name = "Number Softwares",palette = "Set1") +
  theme(text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(title="MCC and TrueDoubletRate/FalseSingletRate for Each Pool and Software Combination",
      subtitle="Faceted by Number of Softwares in Combination")

ggsave(p_TrueDoubletRate_FalseSingletRate_MCC, filename = paste0(out,"MCC_TrueDoubletRate_FalseSingletRate.png"), height = 6, width = 16)

p_TrueDoubletRate_FalseSingletRate_MCC_numberind <- ggplot(Benchmarking_wide, aes(x=TrueDoubletRate_FalseSingletRate,y=MCC,color=factor(N))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~Combination_Type, nrow = 1) +
  scale_color_brewer(name = "Number Individuals in Pool", palette = "Set1")+
  theme(text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "MCC and TrueDoubletRate/FalseSingletRate for Each Pool and Software Combination", 
      subtitle = "Faceted by Number of Softwares in Combination")

ggsave(p_TrueDoubletRate_FalseSingletRate_MCC_numberind, filename = paste0(out,"MCC_TrueDoubletRate_FalseSingletRate_number_individuals.png"), height = 6, width = 16)

p_Accuracy_MCC_combGroup <- ggplot(Benchmarking_wide, aes(x=Accuracy,y=MCC,color=factor(Combination_Groups))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~Combination_Type, nrow = 1) +
  scale_color_brewer(name = "Demultiplexing Doublet Detection Combination", palette = "Set1")+
  theme(text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("MCC and Accuracy for Each Pool and Software Combination")

ggsave(p_Accuracy_MCC_combGroup, filename = paste0(out,"MCC_accuracy_combination_demultiplex_doublet_group_type.png"), height = 6, width = 16)

Benchmarking_wide_singlets <- Benchmarking_wide[which(Benchmarking_wide$N_softwares == 1),]
Benchmarking_wide_singlets <- Benchmarking_wide_singlets[grep(".+_Intersection",Benchmarking_wide_singlets$Comparison),]
Benchmarking_wide_singlets$Comparison <- gsub("^demuxlet_Intersection$","demuxlet", Benchmarking_wide_singlets$Comparison) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.)

p_Accuracy_MCC_SingleSoftware <- ggplot(Benchmarking_wide_singlets, aes(x=Accuracy,y=MCC,color=factor(Comparison))) +
  geom_point(size = 1, alpha = 1)+
  theme_classic() +
  scale_color_brewer(name = "Software", palette = "Set1")+
  theme(text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("MCC and Accuracy for Each Pool and Individual Software")

ggsave(p_Accuracy_MCC_SingleSoftware, filename = paste0(out,"MCC_accuracy_single_software.png"), width = 8, height = 7)


p_TrueDoubletRate_FalseSingletRate_MCC_SingleSoftware <- ggplot(Benchmarking_wide_singlets, aes(x=TrueDoubletRate_FalseSingletRate,y=MCC,color=factor(Comparison))) +
  geom_point(size = 1, alpha = 1)+
  theme_classic() +
  scale_color_brewer(name = "Software", palette = "Set1")+
  theme(text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("MCC and TrueDoubletRate/FalseSingletRate for Each Pool and Each Individual Software")

ggsave(p_TrueDoubletRate_FalseSingletRate_MCC_SingleSoftware, filename = paste0(out,"TrueDoubletRate_FalseSingletRate_MCC_SingleSoftware.png"), width = 8, height = 7)

### Same plots but with TPR (sensitivity) and TNR (specificity)
p_TPR_TNR <- ggplot(Benchmarking_wide, aes(x=TrueSingletRate,y=TrueDoubletRate,color=factor(N_softwares))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~Combination_Type, nrow = 1) +
  scale_color_brewer(name = "Number Softwares",palette = "Set1") +
  theme(text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "True Singlet Rate and True Positive Rate for Each Pool and Software Combination", 
      subtitle = "Faceted by Number of Softwares in Combination")

ggsave(p_TPR_TNR, filename = paste0(out,"TPR_TNR.png"), width = 16, height = 6)

p_TPR_TNR_numberind <- ggplot(Benchmarking_wide, aes(x=TrueSingletRate,y=TrueDoubletRate,color=factor(N))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~Combination_Type, nrow = 1) +
  scale_color_brewer(name = "Number Individuals in Pool", palette = "Set1")+
  theme(text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "True Singlet Rate and True Positive Rate for Each Pool and Software Combination", 
      subtitle = "Faceted by Number of Softwares in Combination")

ggsave(p_TPR_TNR_numberind, filename = paste0(out,"TPR_TNR_number_individuals.png"), width = 16, height = 6)

p_TPR_TNR_SingleSoftware <- ggplot(Benchmarking_wide_singlets, aes(x=TrueSingletRate,y=TrueDoubletRate,color=Comparison)) +
  geom_point(size = 1, alpha = 1)+
  theme_classic() +
  scale_color_brewer(name = "Software", palette = "Set1")+
  theme(text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "True Singlet Rate and True Positive Rate for Each Pool and Individual Software")

ggsave(p_TPR_TNR_SingleSoftware, filename = paste0(out,"TPR_TNR_SingelSoftware.png"), width = 8, height = 7)


p_TPR_TNR_combGroup <- ggplot(Benchmarking_wide, aes(x=TrueSingletRate,y=TrueDoubletRate,color=factor(Combination_Groups))) +
  geom_point(size = 0.5, alpha = 0.5)+
  theme_classic() +
  facet_wrap(~N_softwares ) +
  scale_color_brewer(name = "Demultiplexing Doublet Detection Combination", palette = "Set1")+
  theme(text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="top", 
      legend.box = "horizontal",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("True Singlet Rate and True Positive Rate for Each Pool and Software Combination")

ggsave(p_TPR_TNR_combGroup, filename = paste0(out,"TPR_TNR_combination_demultiplex_doublet_group_type.png"))

pca <- prcomp(Benchmarking_wide[,3:ncol(Benchmarking_wide)][!rowSums(is.na(Benchmarking_wide[,3:ncol(Benchmarking_wide)])) > 0,])
pca_figure <- ggbiplot(pca)
pca_figure <- autoplot(pca, data = Benchmarking_wide[,3:ncol(Benchmarking_wide)][!rowSums(is.na(Benchmarking_wide[,3:ncol(Benchmarking_wide)])) > 0,], colour = 'MCC')
ggsave(pca_figure, filename = paste0(out,"Assessment_PCA.png"))



### Summarize to average for each assessment in each comparison across all pools ###
Summary_Benchmarking_mean <- Benchmarking_wide %>%
  group_by(Comparison) %>%
  summarize(TrueSingletRate = mean(TrueSingletRate, na.rm = TRUE),
            FalseSingletRate = mean(FalseSingletRate, na.rm = TRUE),
            TrueDoubletRate = mean(TrueDoubletRate, na.rm = TRUE),
            FalseDoubletRate = mean(FalseDoubletRate, na.rm = TRUE),
            MCC = mean(MCC, na.rm = TRUE),
            Accuracy = mean(Accuracy, na.rm = TRUE))

Summary_Benchmarking_sd <- Benchmarking_wide %>%
  group_by(Comparison) %>%
  summarize(TrueSingletRate = sd(TrueSingletRate, na.rm = TRUE),
            FalseSingletRate = sd(FalseSingletRate, na.rm = TRUE),
            TrueDoubletRate = sd(TrueDoubletRate, na.rm = TRUE),
            FalseDoubletRate = sd(FalseDoubletRate, na.rm = TRUE),
            MCC = sd(MCC, na.rm = TRUE),
            Accuracy = sd(Accuracy, na.rm = TRUE))

Summary_Benchmarking_mean_long <- pivot_longer(Summary_Benchmarking_mean, -c(Comparison), names_to = "Metric", values_to="Mean")
Summary_Benchmarking_sd_long <- pivot_longer(Summary_Benchmarking_sd, -c(Comparison), names_to = "Metric", values_to="SD")

Summary_Benchmarking_long <- left_join(Summary_Benchmarking_mean_long,Summary_Benchmarking_sd_long)

Summary_Benchmarking_mean_top100MCC <- Summary_Benchmarking_mean[order(-Summary_Benchmarking_mean$MCC),]
Summary_Benchmarking_mean_top100MCC <- Summary_Benchmarking_mean_top100MCC[1:100,]
top100MCC_order <- Summary_Benchmarking_mean_top100MCC$Comparison
Summary_Benchmarking_mean_top100MCC_long <- Summary_Benchmarking_mean_top100MCC[,"Comparison"]
Summary_Benchmarking_mean_top100MCC_long <- inner_join(Summary_Benchmarking_mean_top100MCC_long, Summary_Benchmarking_long)
Summary_Benchmarking_mean_top100MCC_long$Metric <- factor(Summary_Benchmarking_mean_top100MCC_long$Metric, levels = rev(c("MCC","Accuracy","TrueSingletRate","TrueDoubletRate","FalseSingletRate","FalseDoubletRate")))
Summary_Benchmarking_mean_top100MCC_long$Comparison <- factor(Summary_Benchmarking_mean_top100MCC_long$Comparison, levels = unique(Summary_Benchmarking_mean_top100MCC_long$Comparison))
Summary_Benchmarking_mean_top100MCC_long$Comparison <- gsub("^demuxlet_Intersection$","demuxlet", Summary_Benchmarking_mean_top100MCC_long$Comparison) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.)
Summary_Benchmarking_mean_top100MCC_long$Comparison <- factor(Summary_Benchmarking_mean_top100MCC_long$Comparison ,levels = gsub("^demuxlet_Intersection$","demuxlet", top100MCC_order) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.))


Summary_Benchmarking_mean_bottom100MCC <- Summary_Benchmarking_mean[order(Summary_Benchmarking_mean$MCC),]
Summary_Benchmarking_mean_bottom100MCC <- Summary_Benchmarking_mean_bottom100MCC[1:100,]
bottom100MCC_order <- Summary_Benchmarking_mean_bottom100MCC$Comparison
Summary_Benchmarking_mean_bottom100MCC <- Summary_Benchmarking_mean_bottom100MCC[,"Comparison"]
Summary_Benchmarking_mean_bottom100MCC <- inner_join(Summary_Benchmarking_mean_bottom100MCC, Summary_Benchmarking_long)
Summary_Benchmarking_mean_bottom100MCC$Metric <- factor(Summary_Benchmarking_mean_bottom100MCC$Metric, levels = rev(c("MCC","Accuracy","TrueSingletRate","TrueDoubletRate","FalseSingletRate","FalseDoubletRate")))
Summary_Benchmarking_mean_bottom100MCC$Comparison <- factor(Summary_Benchmarking_mean_bottom100MCC$Comparison, levels = unique(Summary_Benchmarking_mean_bottom100MCC$Comparison))
Summary_Benchmarking_mean_bottom100MCC$Comparison <- gsub("^demuxlet_Intersection$","demuxlet", Summary_Benchmarking_mean_bottom100MCC$Comparison) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.)
Summary_Benchmarking_mean_bottom100MCC$Comparison <- factor(Summary_Benchmarking_mean_bottom100MCC$Comparison ,levels = gsub("^demuxlet_Intersection$","demuxlet", bottom100MCC_order) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.))

Summary_Benchmarking_mean_top100MCC_long$End <- "Top"
Summary_Benchmarking_mean_bottom100MCC$End <- "Bottom"

Summary_Benchmarking_mean_top_bottom <- rbind(Summary_Benchmarking_mean_top100MCC_long[which(Summary_Benchmarking_mean_top100MCC_long$Comparison %in% (gsub("^demuxlet_Intersection$","demuxlet", top100MCC_order[1:50]) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.))),], Summary_Benchmarking_mean_bottom100MCC[which(Summary_Benchmarking_mean_bottom100MCC$Comparison %in% c(gsub("^demuxlet_Intersection$","demuxlet", bottom100MCC_order[51:100]) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.))),])
Summary_Benchmarking_mean_top_bottom$Comparison <- factor(Summary_Benchmarking_mean_top_bottom$Comparison, levels = c((gsub("^demuxlet_Intersection$","demuxlet", top100MCC_order[1:50]) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.)),(gsub("^demuxlet_Intersection$","demuxlet", bottom100MCC_order[51:100]) %>%
                                                          gsub("^freemuxlet_Intersection$","freemuxlet",.) %>%
                                                          gsub("^scSplit_Intersection$","scSplit",.) %>%
                                                          gsub("^souporcell_Intersection$","souporcell",.) %>%
                                                          gsub("^vireo_Intersection$","vireo",.) %>%
                                                          gsub("^DoubletDetection_Intersection$","DoubletDetection",.) %>%
                                                          gsub("^DoubletFinder_Intersection$","DoubletFinder",.) %>%
                                                          gsub("^scds_Intersection$","scds",.) %>%
                                                          gsub("^scrublet_Intersection$","scrublet",.))))
Summary_Benchmarking_mean_top_bottom$End <- factor(Summary_Benchmarking_mean_top_bottom$End, levels = c("Top","Bottom"))

mid <- (max(Summary_Benchmarking_mean_top_bottom$Mean) + min(Summary_Benchmarking_mean_top_bottom$Mean))/2
p_CircleHeat <- ggplot(Summary_Benchmarking_mean_top_bottom, aes(x=Comparison, y = Metric)) +
  geom_point(aes(colour = Mean,  size =SD)) +
  theme_classic() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank()) +
  scale_color_gradient2(midpoint = 0.5, low = "blue", mid = "white",
                            high = "red", space = "Lab" ) +
  facet_wrap(~End, ncol = 2, scales = "free_x")

ggsave(p_CircleHeat, filename = paste0(out,"CircleHeatmap_Top100.png"), width = 16, height = 8)


##### Make a faceted boxplot for the top 100 combinations #####
Benchmarking_top100 <- Summary_Benchmarking_mean_top100MCC[,"Comparison"]
Benchmarking_top100 <- inner_join(Benchmarking_top100, Benchmarking)
Benchmarking_top100$N_softwares <- (str_count(gsub("scrublet_DoubletUnion","", Benchmarking_top100$Comparison) %>% gsub("_Intersection","",.) %>% gsub("_HalfIntersection","", .) %>% gsub("_DoubletRemainder","", .), pattern = "_") + 1)
Benchmarking_top100$Combination_Type <- gsub("demuxlet_","",Benchmarking_top100$Comparison) %>%
                                        gsub("freemuxlet_","",.) %>%
                                        gsub("scSplit_","",.) %>%
                                        gsub("souporcell_","",.) %>%
                                        gsub("vireo_","",.) %>%
                                        gsub("DoubletDetection_","",.) %>%
                                        gsub("DoubletFinder_","",.) %>%
                                        gsub("scds_","",.) %>%
                                        gsub("scrublet_","",.)
Benchmarking_top100$N <- gsub("Size","",Benchmarking_top100$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
Benchmarking_top100$N <- as.numeric(as.character(Benchmarking_top100$N))
Benchmarking_top100$Combination_Type <- factor(Benchmarking_top100$Combination_Type, levels = c("Intersection","DoubletUnion","DoubletRemainder","HalfIntersection"))
Benchmarking_top100$Assessment <- factor(Benchmarking_top100$Assessment, levels = c("MCC","Accuracy","TrueSingletRate","TrueDoubletRate","FalseSingletRate","FalseDoubletRate"))
levels(Benchmarking_top100$Assessment) <- c("MCC","Accuracy","True Singlet\nRate","True Doublet\nRate","False Singlet\nRate","False Doublet\nRate")
Benchmarking_top20 <- Benchmarking_top100[which(Benchmarking_top100$Comparison %in% top100MCC_order[1:20]),]
Benchmarking_top20$Comparison <- factor(Benchmarking_top20$Comparison, levels = top100MCC_order[1:20])

pMetricsBox <- ggplot(Benchmarking_top20, aes(x = Comparison, y = Value, color = factor(N))) +
  geom_jitter(size = 0.5, alpha = 0.5) +
  facet_wrap(~Assessment, ncol = 1,strip.position="left", scales="free_y") +
  theme_linedraw() +
  scale_color_viridis(option = "plasma",discrete = TRUE, name = "Number\nIndividuals") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
      axis.title.x = element_blank(),
      text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="left", 
      legend.box = "vertical",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
      plot.margin=unit(c(5.5, 5.5, 5.5, 40), "points")) +
  labs(title = "Top 20 Software Combinations")

ggsave(pMetricsBox, filename = paste0(out,"MetricsJitterTop100.png"), width = 16, height = 8)

Summary_Benchmarking_mean$TrueDoubletRate_FalseSingletRate <- Summary_Benchmarking_mean$TrueDoubletRate/Summary_Benchmarking_mean$FalseSingletRate
Summary_Benchmarking_lowFalseSinglet_highTrueDoublet20 <- Summary_Benchmarking_mean[order(-Summary_Benchmarking_mean$TrueDoubletRate_FalseSingletRate)[1:50],]
Benchmarking_lowFalseSinglet20 <- Summary_Benchmarking_lowFalseSinglet_highTrueDoublet20[,"Comparison"]
Benchmarking_lowFalseSinglet20 <- inner_join(Benchmarking_lowFalseSinglet20, Benchmarking)
Benchmarking_lowFalseSinglet20$N <- gsub("Size","",Benchmarking_lowFalseSinglet20$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
Benchmarking_lowFalseSinglet20$Comparison <- factor(Benchmarking_lowFalseSinglet20$Comparison, levels = Summary_Benchmarking_lowFalseSinglet_highTrueDoublet20$Comparison)

pMetricsLowFalseSingletRate <- ggplot(Benchmarking_lowFalseSinglet20, aes(x = Comparison, y = Value, color = factor(N))) +
  geom_jitter(size = 0.5, alpha = 0.5) +
  facet_wrap(~Assessment, ncol = 1,strip.position="left", scales="free_y") +
  theme_linedraw() +
  scale_color_viridis(option = "plasma",discrete = TRUE, name = "Number\nIndividuals") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
      axis.title.x = element_blank(),
      text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="left", 
      legend.box = "vertical",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
      plot.margin=unit(c(5.5, 5.5, 5.5, 40), "points")) +
  labs(title = "Top 50 Software Combinations by the (True Doublet Rate)/(False Singlet Rate)")

ggsave(pMetricsLowFalseSingletRate, filename = paste0(out,"MetricsJitterLowSinglet20.png"), width = 16, height = 8)

##### Figures for numbers per pool for each classification #####
Benchmarking_counts_wide_prop <- Benchmarking_counts_wide[,c("Pool","Comparison","N_softwares","Combination_Type")]
Benchmarking_counts_wide_prop$TrueSingletProp <- Benchmarking_counts_wide$TrueSinglet/Benchmarking_counts_wide$Total
Benchmarking_counts_wide_prop$TrueDoubletProp <- Benchmarking_counts_wide$TrueDoublet/Benchmarking_counts_wide$Total
Benchmarking_counts_wide_prop$FalseSingletProp <- Benchmarking_counts_wide$FalseSinglet/Benchmarking_counts_wide$Total
Benchmarking_counts_wide_prop$FalseDoubletProp <- Benchmarking_counts_wide$FalseDoublet/Benchmarking_counts_wide$Total
Benchmarking_counts_wide_prop$UnassignedProp <- (Benchmarking_counts_wide$Total - Benchmarking_counts_wide$TrueSinglet - Benchmarking_counts_wide$TrueDoublet - Benchmarking_counts_wide$FalseSinglet - Benchmarking_counts_wide$FalseDoublet)/Benchmarking_counts_wide$Total


SummaryBenchmarking_counts_mean <- Benchmarking_counts_wide_prop %>%
  group_by(Comparison) %>%
  summarize(TrueSingletProp = mean(TrueSingletProp, na.rm = TRUE),
            TrueDoubletProp = mean(TrueDoubletProp, na.rm = TRUE),
            FalseSingletProp = mean(FalseSingletProp, na.rm = TRUE),
            FalseDoubletProp = mean(FalseDoubletProp, na.rm = TRUE),
            UnassignedProp = mean(UnassignedProp, na.rm = TRUE))

SummaryBenchmarking_counts_sd <- Benchmarking_counts_wide_prop %>%
  group_by(Comparison) %>%
  summarize(TrueSingletProp = sd(TrueSingletProp, na.rm = TRUE),
            TrueDoubletProp = sd(TrueDoubletProp, na.rm = TRUE),
            FalseSingletProp = sd(FalseSingletProp, na.rm = TRUE),
            FalseDoubletProp = sd(FalseDoubletProp, na.rm = TRUE),
            UnassignedProp = sd(UnassignedProp, na.rm = TRUE))

SummaryBenchmarking_counts_mean$True <- SummaryBenchmarking_counts_mean$TrueDoubletProp + SummaryBenchmarking_counts_mean$TrueSingletProp
SummaryBenchmarking_counts_mean$False <- SummaryBenchmarking_counts_mean$FalseDoubletProp + SummaryBenchmarking_counts_mean$FalseSingletProp

SummaryBenchmarking_counts_mean_long <- pivot_longer(SummaryBenchmarking_counts_mean, -c(Comparison, True), names_to = "Metric", values_to="Mean")
SummaryBenchmarking_counts_sd_long <- pivot_longer(SummaryBenchmarking_counts_sd, -c(Comparison), names_to = "Metric", values_to="SD")

SummaryBenchmarking_counts_long <- left_join(SummaryBenchmarking_counts_mean_long,SummaryBenchmarking_counts_sd_long)
SummaryBenchmarking_counts_long <- inner_join(SummaryBenchmarking_counts_long,Benchmarking_counts_wide_prop[,c("Comparison","N_softwares","Combination_Type")])
SummaryBenchmarking_counts_long <- unique(SummaryBenchmarking_counts_long)

SummaryBenchmarking_counts_long <- SummaryBenchmarking_counts_long[order(-SummaryBenchmarking_counts_long$True),]

Benchmarking_counts_long_prop_top20 <- SummaryBenchmarking_counts_long[1:100,]

pProportion <- ggplot(Benchmarking_counts_long_prop_top20, aes(Comparison, Mean, fill = Metric)) +
  geom_bar(stat = "identity", position = "stack") 

ggsave(pProportion, filename = paste0(out,"Proportion_Bar.png"))

SummaryBenchmarking_counts_maxTrue_True <- SummaryBenchmarking_counts_long[order(-SummaryBenchmarking_counts_long$True)[1:120],]
# SummaryBenchmarking_counts_maxTrue_True<- SummaryBenchmarking_counts_maxTrue[,c("Comparison","Metric","Mean","SD")]
SummaryBenchmarking_counts_maxTrue_True$Facet <- "Proportion"
SummaryBenchmarking_counts_maxTrue_MCC <- inner_join(SummaryBenchmarking_counts_maxTrue_True[,"Comparison"],Benchmarking, by = "Comparison")
SummaryBenchmarking_counts_maxTrue_MCC$Facet <- "Proportion"


SummaryBenchmarking_counts_maxMCC_MCC <- inner_join(Benchmarking_top20[,"Comparison"],Benchmarking, by = "Comparison")
SummaryBenchmarking_counts_maxMCC_MCC$Facet <- "MCC"
SummaryBenchmarking_counts_maxMCC_True <- inner_join(Benchmarking_top20[,"Comparison"],SummaryBenchmarking_counts_long, by = "Comparison")
SummaryBenchmarking_counts_maxMCC_True$Facet <- "MCC"

order <- c(as.character(unique(SummaryBenchmarking_counts_long$Comparison[order(-SummaryBenchmarking_counts_long$True)[1:120]])),as.character(unique(Benchmarking_top20$Comparison)))
SummaryBenchmarking_counts_True <- rbind(SummaryBenchmarking_counts_maxTrue_True, SummaryBenchmarking_counts_maxMCC_True)
SummaryBenchmarking_counts_True <- SummaryBenchmarking_counts_True[which(SummaryBenchmarking_counts_True$Metric == "TrueSingletProp" | 
                                                                          SummaryBenchmarking_counts_True$Metric == "TrueDoubletProp" |
                                                                          SummaryBenchmarking_counts_True$Metric == "FalseSingletProp" |
                                                                          SummaryBenchmarking_counts_True$Metric == "FalseDoubletProp"),]
SummaryBenchmarking_counts_True <- unique(SummaryBenchmarking_counts_True)
SummaryBenchmarking_counts_True$Facet <- "MCC"
SummaryBenchmarking_counts_True$Comparison <- factor(SummaryBenchmarking_counts_True$Comparison, levels = unique(order))
SummaryBenchmarking_counts_MCC <- rbind(SummaryBenchmarking_counts_maxTrue_MCC, SummaryBenchmarking_counts_maxMCC_MCC)
SummaryBenchmarking_counts_MCC <- SummaryBenchmarking_counts_MCC[which(SummaryBenchmarking_counts_MCC$Assessment == "MCC" | SummaryBenchmarking_counts_MCC$Assessment == "Accuracy"),]
SummaryBenchmarking_counts_MCC <- unique(SummaryBenchmarking_counts_MCC)
SummaryBenchmarking_counts_MCC$Comparison <- factor(SummaryBenchmarking_counts_MCC$Comparison, levels = order)

pMCC <- ggplot(SummaryBenchmarking_counts_MCC, aes(Comparison, Value, color = (gsub("Size","", Pool) %>% gsub("SimulatedPools[0-9]+","",.)))) +
  geom_jitter(size = 0.5, alpha = 0.5) +
  facet_grid(Assessment~Facet, scales="free") +
  theme_linedraw() +
  scale_color_viridis(option = "plasma",discrete = TRUE, name = "Number\nIndividuals") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
      axis.title.x = element_blank(),
      text = element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="left", 
      legend.box = "vertical",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
      plot.margin=unit(c(5.5, 5.5, 5.5, 40), "points")) +
  labs(title = "Top 20 Software Combinations")

ggsave(pMCC, filename = paste0(out,"MCC_part_figure.png"), width = 20, height = 8)


pProp <- ggplot(SummaryBenchmarking_counts_True, aes(Comparison, Mean)) +
  geom_bar(stat = "identity", position = "stack",aes(fill = Metric)) +
  geom_errorbar(aes(ymax=Mean-SD, ymin=Mean+SD), width=.2,stat = "identity",position = "identity") +
  theme_linedraw() +
  facet_wrap(~Facet, ncol = 2, scales="free_x") +
  scale_fill_brewer(name = "Software", palette = "Set1")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
  text = element_text(size=14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position="left", 
      legend.box = "verticle",
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        plot.margin=unit(c(5.5, 5.5, 5.5, 40), "points")) 

ggsave(pProp, filename = paste0(out,"Proportion_Bar_part_figure.png"), width = 20, height = 8)





# ##### Calculate accuracy of each software and combinations of softwares #####
# accuracy_list <- lapply(results_list_wide, function(x){
#   df <- as.data.frame(matrix(nrow=1, ncol=6))
#   colnames(df) <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet")
#   df[1,"demuxlet"] <- length(which(x$demuxlet_Assignment == x$Individual))/nrow(x)
#   df[1,"freemuxlet"] <- length(which(x$freemuxlet_Assignment == x$Individual))/nrow(x)
#   df[1,"scSplit"] <- length(which(x$scSplit_Assignment == x$Individual))/nrow(x)
#   df[1,"souporcell"] <- length(which(x$souporcell_Assignment == x$Individual))/nrow(x)
#   df[1,"vireo"] <- length(which(x$vireo_Assignment == x$Individual))/nrow(x)
#   df[1,"scrublet"] <- length(which(x$scrublet_Assignment == x$DropletType))/nrow(x)
#   return(df)
# }) 
# names(accuracy_list) <- names(key)

# accuracy_list <- lapply(names(accuracy_list), function(x){
#   accuracy_list[[x]]$Pool <- x
#   return(accuracy_list[[x]])
# })
# names(accuracy_list) <- names(key)

# accuracy_df <- do.call(rbind,accuracy_list)
# accuracy_df_long <- pivot_longer(accuracy_df, cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"), values_to = "Accuracy", names_to = "Software")
# accuracy_df_long$Size <- gsub("_SimulatedPools[0-9]+","",accuracy_df_long$Pool)
# accuracy_df_long$Size <- factor(accuracy_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))
# accuracy_df_long$Software <- factor(accuracy_df_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"))

# p_accuracy <- ggplot(accuracy_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   a_primary_fill() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Accuracy of Demultiplexing Softwares") +
#   ylab("Accuracy") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_accuracy, filename=paste0(out,"SingleSoftwareAccuracies.png"),width = 16, height = 9)


# ##### Calculate TPR and TNR of each software #####
# ### TPR ###
# TPR_list <- lapply(results_list_wide, function(x){
#   df <- as.data.frame(matrix(nrow=1, ncol=5))
#   colnames(df) <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo")
#   df[1,"demuxlet"] <- length(which(x$demuxlet_Assignment == x$Individual & x$demuxlet_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"freemuxlet"] <- length(which(x$freemuxlet_Assignment == x$Individual & x$freemuxlet_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"scSplit"] <- length(which(x$scSplit_Assignment == x$Individual & x$scSplit_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"souporcell"] <- length(which(x$souporcell_Assignment == x$Individual & x$souporcell_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"vireo"] <- length(which(x$vireo_Assignment == x$Individual & x$vireo_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"scrublet"] <- length(which(x$scrublet_Assignment == x$DropletType & x$scrublet_Assignment == "singlet"))/length(which(x$DropletType == "singlet"))
#   return(df)
# })
# names(TPR_list) <- names(key)

# TPR_list <- lapply(names(TPR_list), function(x){
#   TPR_list[[x]]$Pool <- x
#   return(TPR_list[[x]])
# })

# TPR_df <- do.call(rbind, TPR_list)
# TPR_df_long <- pivot_longer(TPR_df, cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"), values_to = "Accuracy", names_to = "Software")
# TPR_df_long$Size <- gsub("_SimulatedPools[0-9]+","",TPR_df_long$Pool)
# TPR_df_long$Size <- factor(TPR_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))
# TPR_df_long$Software <- factor(TPR_df_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"))

# p_TPR <- ggplot(TPR_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   a_primary_fill() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Sensitivity of Demultiplexing Softwares") +
#   ylab("TPR") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_TPR, filename=paste0(out,"SingleSoftwareTPR.png"),width = 16, height = 9)

# ### TNR ###
# TNR_list <- lapply(results_list_wide, function(x){
#   df <- as.data.frame(matrix(nrow=1, ncol=5))
#   colnames(df) <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo")
#   df[1,"demuxlet"] <- length(which(x$demuxlet_Assignment == x$Individual & x$demuxlet_DropletType == "doublet"))/length(which(x$DropletType == "doublet"))
#   df[1,"freemuxlet"] <- length(which(x$freemuxlet_Assignment == x$Individual & x$freemuxlet_DropletType == "doublet"))/length(which(x$DropletType == "doublet"))
#   df[1,"scSplit"] <- length(which(x$scSplit_Assignment == x$Individual & x$scSplit_DropletType == "doublet"))/length(which(x$DropletType == "doublet"))
#   df[1,"souporcell"] <- length(which(x$souporcell_Assignment == x$Individual & x$souporcell_DropletType == "doublet"))/length(which(x$DropletType == "doublet"))
#   df[1,"vireo"] <- length(which(x$vireo_Assignment == x$Individual & x$vireo_DropletType == "doublet"))/length(which(x$DropletType == "doublet"))
#   df[1,"scrublet"] <- length(which(x$scrublet_Assignment == x$DropletType & x$scrublet_Assignment == "doublet"))/length(which(x$DropletType == "doublet"))
#   return(df)
# })
# names(TNR_list) <- names(key)

# TNR_list <- lapply(names(TNR_list), function(x){
#   TNR_list[[x]]$Pool <- x
#   return(TNR_list[[x]])
# })

# TNR_df <- do.call(rbind, TNR_list)
# TNR_df_long <- pivot_longer(TNR_df, cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"), values_to = "Accuracy", names_to = "Software")
# TNR_df_long$Size <- gsub("_SimulatedPools[0-9]+","",TNR_df_long$Pool)
# TNR_df_long$Size <- factor(TNR_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))
# TNR_df_long$Software <- factor(TNR_df_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"))

# p_TNR <- ggplot(TNR_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   a_primary_fill() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Specificity of Demultiplexing Softwares") +
#   ylab("TNR") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_TNR, filename=paste0(out,"SingleSoftwareTNR.png"),width = 16, height = 9)



# # ##### Calculate accuracy of each software considering unassigned as a "removal" #####
# # results_reomove_list_wide <- lapply(results_list_wide, function(x){
# #   x$demuxlet_Assignment <- gsub("unassigned","doublet", x$demuxlet_Assignment) %>% gsub("unsure","doublet",.)
# #   x$freemuxlet_Assignment <- gsub("unassigned","doublet", x$freemuxlet_Assignment) %>% gsub("unsure","doublet",.)
# #   x$scSplit_Assignment <- gsub("unassigned","doublet", x$scSplit_Assignment) %>% gsub("unsure","doublet",.)
# #   x$souporcell_Assignment <- gsub("unassigned","doublet", x$souporcell_Assignment) %>% gsub("unsure","doublet",.)
# #   x$vireo_Assignment <- gsub("unassigned","doublet", x$vireo_Assignment) %>% gsub("unsure","doublet",.)
# #   return(x)
# # })
# # names(results_reomove_list_wide) <- names(key)

# # accuracy_remove_list <- lapply(results_reomove_list_wide, function(x){
# #   df <- as.data.frame(matrix(nrow=1, ncol=5))
# #   colnames(df) <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo")
# #   df[1,"demuxlet"] <- length(which(x$demuxlet_Assignment == x$Individual))/nrow(x)
# #   df[1,"freemuxlet"] <- length(which(x$freemuxlet_Assignment == x$Individual))/nrow(x)
# #   df[1,"scSplit"] <- length(which(x$scSplit_Assignment == x$Individual))/nrow(x)
# #   df[1,"souporcell"] <- length(which(x$souporcell_Assignment == x$Individual))/nrow(x)
# #   df[1,"vireo"] <- length(which(x$vireo_Assignment == x$Individual))/nrow(x)
# #   return(df)
# # }) 
# # names(accuracy_remove_list) <- names(key)

# # accuracy_remove_list <- lapply(names(accuracy_remove_list), function(x){
# #   accuracy_remove_list[[x]]$Pool <- x
# #   return(accuracy_remove_list[[x]])
# # })
# # names(accuracy_remove_list) <- names(key)

# # accuracy_remove_df <- do.call(rbind,accuracy_remove_list)
# # accuracy_remove_df_long <- pivot_longer(accuracy_remove_df, cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo"), values_to = "Accuracy", names_to = "Software")
# # accuracy_remove_df_long$Size <- gsub("_SimulatedPools[0-9]+","",accuracy_remove_df_long$Pool)
# # accuracy_remove_df_long$Size <- factor(accuracy_remove_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))

# # p_accuracy_remove <- ggplot(accuracy_remove_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
# #   geom_violin() +
# #   a_primary_fill() +
# #   theme_classic() +
# #     theme(text = element_text(size=18),
# #       plot.title = element_text(hjust = 0.5),
# #       legend.background = element_blank(),
# #         legend.box.background = element_rect(colour = "black")) +
# #   ggtitle("Accuracy of Demultiplexing Softwares") +
# #   ylab("Accuracy") +
# #   xlab("Size of Simulated Pool")
# # ggsave(plot = p_accuracy_remove, filename=paste0(out,"SingleSoftwareAccuracies_unassigneddoublets.png"),width = 16, height = 9)


# ##### Calculate the accuracy of the intersection of two of the softwares #####
# paired_df_list <- lapply(results_list_wide, function(x){
#   temp <- x[,c("Barcode","DropletType","Individual","DoubletType")]
#   temp[,"demuxlet_freemuxlet"] <- ifelse(x$demuxlet_Assignment ==x$freemuxlet_Assignment, x$demuxlet_Assignment,"discordant")
#   temp[,"demuxlet_scSplit"] <- ifelse(x$demuxlet_Assignment ==x$scSplit_Assignment, x$demuxlet_Assignment,"discordant")
#   temp[,"demuxlet_souporcell"] <- ifelse(x$demuxlet_Assignment ==x$souporcell_Assignment, x$demuxlet_Assignment,"discordant")
#   temp[,"demuxlet_vireo"] <- ifelse(x$demuxlet_Assignment ==x$vireo_Assignment, x$demuxlet_Assignment,"discordant")
#   temp[,"demuxlet_scrublet"] <- ifelse(x$demuxlet_DropletType == x$scrublet_Assignment, x$demuxlet_Assignment,"discordant")
#   temp[,"freemuxlet_scSplit"] <- ifelse(x$freemuxlet_Assignment ==x$scSplit_Assignment, x$freemuxlet_Assignment,"discordant")
#   temp[,"freemuxlet_souporcell"] <- ifelse(x$freemuxlet_Assignment ==x$souporcell_Assignment, x$freemuxlet_Assignment,"discordant")
#   temp[,"freemuxlet_vireo"] <- ifelse(x$freemuxlet_Assignment ==x$vireo_Assignment, x$freemuxlet_Assignment,"discordant")
#   temp[,"freemuxlet_scrublet"] <- ifelse(x$freemuxlet_DropletType ==x$scrublet_Assignment, x$freemuxlet_Assignment,"discordant")
#   temp[,"scSplit_souporcell"] <- ifelse(x$scSplit_Assignment ==x$souporcell_Assignment, x$scSplit_Assignment,"discordant")
#   temp[,"scSplit_vireo"] <- ifelse(x$scSplit_Assignment ==x$vireo_Assignment, x$scSplit_Assignment,"discordant")
#   temp[,"scSplit_scrublet"] <- ifelse(x$scSplit_DropletType ==x$scrublet_Assignment, x$scSplit_Assignment,"discordant")
#   temp[,"souporcell_vireo"] <- ifelse(x$vireo_Assignment ==x$souporcell_Assignment, x$vireo_Assignment,"discordant")
#   temp[,"souporcell_scrublet"] <- ifelse(x$souporcell_DropletType == x$scrublet_Assignment, x$souporcell_Assignment,"discordant")
#   temp[,"vireo_scrublet"] <- ifelse(x$vireo_DropletType == x$scrublet_Assignment, x$vireo_Assignment,"discordant")
#   return(temp)
# })


# accuracy2_list <- lapply(paired_df_list, function(x){
#   df <- as.data.frame(matrix(nrow=1, ncol=15))
#   colnames(df) <- c("demuxlet_freemuxlet","demuxlet_scSplit","demuxlet_souporcell","demuxlet_vireo","demuxlet_scrublet","freemuxlet_scSplit","freemuxlet_souporcell","freemuxlet_vireo","freemuxlet_scrublet","scSplit_souporcell","scSplit_vireo","scSplit_scrublet","souporcell_vireo","souporcell_scrublet","vireo_scrublet")
#   df[1,"demuxlet_freemuxlet"] <- length(which(x$demuxlet_freemuxlet == x$Individual))/nrow(x)
#   df[1,"demuxlet_scSplit"] <- length(which(x$demuxlet_scSplit == x$Individual))/nrow(x)
#   df[1,"demuxlet_souporcell"] <- length(which(x$demuxlet_souporcell == x$Individual))/nrow(x)
#   df[1,"demuxlet_vireo"] <- length(which(x$demuxlet_vireo == x$Individual))/nrow(x)
#   df[1,"demuxlet_scrublet"] <- length(which(x$demuxlet_scrublet == x$Individual))/nrow(x)
#   df[1,"freemuxlet_scSplit"] <- length(which(x$freemuxlet_scSplit == x$Individual))/nrow(x)
#   df[1,"freemuxlet_souporcell"] <- length(which(x$freemuxlet_souporcell == x$Individual))/nrow(x)
#   df[1,"freemuxlet_vireo"] <- length(which(x$freemuxlet_vireo == x$Individual))/nrow(x)
#   df[1,"freemuxlet_scrublet"] <- length(which(x$freemuxlet_scrublet == x$Individual))/nrow(x)
#   df[1,"scSplit_souporcell"] <- length(which(x$scSplit_souporcell == x$Individual))/nrow(x)
#   df[1,"scSplit_vireo"] <- length(which(x$scSplit_vireo == x$Individual))/nrow(x)
#   df[1,"scSplit_scrublet"] <- length(which(x$scSplit_scrublet == x$Individual))/nrow(x)
#   df[1,"souporcell_vireo"] <- length(which(x$souporcell_vireo == x$Individual))/nrow(x)
#   df[1,"souporcell_scrublet"] <- length(which(x$souporcell_scrublet == x$Individual))/nrow(x)
#   df[1,"vireo_scrublet"] <- length(which(x$vireo_scrublet == x$Individual))/nrow(x)
#   return(df)
# }) 
# names(accuracy2_list) <- names(key)

# accuracy2_list <- lapply(names(accuracy2_list), function(x){
#   accuracy2_list[[x]]$Pool <- x
#   return(accuracy2_list[[x]])
# })
# names(accuracy2_list) <- names(key)

# accuracy2_df <- do.call(rbind,accuracy2_list)
# accuracy2_df_long <- pivot_longer(accuracy2_df, cols = c("demuxlet_freemuxlet","demuxlet_scSplit","demuxlet_souporcell","demuxlet_vireo","demuxlet_scrublet","freemuxlet_scSplit","freemuxlet_souporcell","freemuxlet_vireo","freemuxlet_scrublet","scSplit_souporcell","scSplit_vireo","scSplit_scrublet","souporcell_vireo","souporcell_scrublet","vireo_scrublet"), values_to = "Accuracy", names_to = "Software")
# accuracy2_df_long$Size <- gsub("_SimulatedPools[0-9]+","",accuracy2_df_long$Pool)
# accuracy2_df_long$Size <- factor(accuracy2_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))
# accuracy2_df_long$Software <- factor(accuracy2_df_long$Software, levels = c("demuxlet_freemuxlet","demuxlet_scSplit","demuxlet_souporcell","demuxlet_vireo","demuxlet_scrublet","freemuxlet_scSplit","freemuxlet_souporcell","freemuxlet_vireo","freemuxlet_scrublet","scSplit_souporcell","scSplit_vireo","scSplit_scrublet","souporcell_vireo","souporcell_scrublet","vireo_scrublet"))

# p_accuracy2 <- ggplot(accuracy2_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Accuracy of Demultiplexing Softwares") +
#   ylab("Accuracy") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_accuracy2, filename=paste0(out,"Accuracy_paired.png"),width = 16, height = 9)

# accruacy_comparison <- rbind(accuracy_df_long,accuracy2_df_long)

# p_accuracy_comparison <- ggplot(accruacy_comparison, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Accuracy of Demultiplexing Softwares") +
#   ylab("Accuracy") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_accuracy_comparison, filename=paste0(out,"Accuracy_paired_single.png"),width = 16, height = 9)

# ##### Calculate TPR and TNR of each software #####
# ### TPR ###
# TPR2_list <- lapply(paired_df_list, function(x){
#   df <- as.data.frame(matrix(nrow=1, ncol=5))
#   colnames(df) <- c("demuxlet","freemuxlet","scSplit","souporcell","vireo")
#   df[1,"demuxlet"] <- length(which(x$demuxlet_Assignment == x$Individual & x$demuxlet_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"freemuxlet"] <- length(which(x$freemuxlet_Assignment == x$Individual & x$freemuxlet_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"scSplit"] <- length(which(x$scSplit_Assignment == x$Individual & x$scSplit_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"souporcell"] <- length(which(x$souporcell_Assignment == x$Individual & x$souporcell_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"vireo"] <- length(which(x$vireo_Assignment == x$Individual & x$vireo_DropletType == "singlet"))/length(which(x$DropletType == "singlet"))
#   df[1,"scrublet"] <- length(which(x$scrublet_Assignment == x$DropletType & x$scrublet_Assignment == "singlet"))/length(which(x$DropletType == "singlet"))
#   return(df)
# })
# names(TPR_list) <- names(key)

# TPR_list <- lapply(names(TPR_list), function(x){
#   TPR_list[[x]]$Pool <- x
#   return(TPR_list[[x]])
# })

# TPR_df <- do.call(rbind, TPR_list)
# TPR_df_long <- pivot_longer(TPR_df, cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"), values_to = "Accuracy", names_to = "Software")
# TPR_df_long$Size <- gsub("_SimulatedPools[0-9]+","",TPR_df_long$Pool)
# TPR_df_long$Size <- factor(TPR_df_long$Size, levels = c("Size8","Size10","Size12","Size14","Size16"))
# TPR_df_long$Software <- factor(TPR_df_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","scrublet"))

# p_TPR <- ggplot(TPR_df_long, aes(x=Size, y = Accuracy, fill = Software)) +
#   geom_violin() +
#   a_primary_fill() +
#   theme_classic() +
#     theme(text = element_text(size=18),
#       plot.title = element_text(hjust = 0.5),
#       legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   ggtitle("Sensitivity of Demultiplexing Softwares") +
#   ylab("TPR") +
#   xlab("Size of Simulated Pool")
# ggsave(plot = p_TPR, filename=paste0(out,"SingleSoftwareTPR.png"),width = 16, height = 9)




# ##### Make a dataframe that has the number of singlets, doublets, unassigned in each pool for each software #####
# results_singlet_doulbet_unnassigned <- lapply(results_list_wide, function(x){
#   temp <- as.data.frame(matrix(nrow=7, ncol=4))
#   colnames(temp) <- c("software","singlets","doublets","unassigned")
#   i <- 1
#   for (software in c("demuxlet","freemuxlet","scSplit","souporcell","vireo")){
#     temp[i,"software"] <- software
#     temp[i,"singlets"] <- nrow(x[which(x[,paste0(software,"_DropletType")] == "singlet"),])
#     temp[i,"doublets"] <- nrow(x[which(x[,paste0(software,"_DropletType")] == "doublet"),])
#     temp[i,"unassigned"] <- nrow(x[which(x[,paste0(software,"_DropletType")] == "unassigned" | is.na(x[,paste0(software,"_DropletType")])),])
#     i <- i + 1
#   }
#   temp[i,"software"] <- "intersection"
#   temp[i,"singlets"] <- length(which(x$demuxlet_DropletType == "singlet" & 
#                                        x$freemuxlet_DropletType == "singlet" &
#                                        x$scSplit_DropletType == "singlet" &
#                                        x$souporcell_DropletType == "singlet" &
#                                        x$vireo_DropletType == "singlet" &
#                                  x$demuxlet_Assignment == x$freemuxlet_Assignment & 
#                                  x$demuxlet_Assignment == x$scSplit_Assignment &
#                                    x$demuxlet_Assignment == x$souporcell_Assignment &
#                                    x$demuxlet_Assignment == x$vireo_Assignment))
#   temp[i,"doublets"] <- length(which(x$demuxlet_DropletType == "doublet" &
#                                        x$freemuxlet_DropletType == "doublet" &
#                                        x$scSplit_DropletType == "doublet" &
#                                        x$souporcell_DropletType == "doublet" &
#                                        x$vireo_DropletType == "doublet"))
#   temp[i,"unassigned"] <- nrow(x) - temp[i,"doublets"] -temp[i,"singlets"]
#   i <- i+1
#   temp[i,"software"] <- "intersection_no_scSplit"
#   temp[i,"singlets"] <- length(which(x$demuxlet_DropletType == "singlet" & 
#                                        x$freemuxlet_DropletType == "singlet" &
#                                        x$souporcell_DropletType == "singlet" &
#                                        x$vireo_DropletType == "singlet" &
#                                        x$demuxlet_Assignment == x$freemuxlet_Assignment & 
#                                        x$demuxlet_Assignment == x$souporcell_Assignment &
#                                        x$demuxlet_Assignment == x$vireo_Assignment))
#   temp[i,"doublets"] <- length(which(x$demuxlet_DropletType == "doublet" &
#                                        x$freemuxlet_DropletType == "doublet" &
#                                        x$souporcell_DropletType == "doublet" &
#                                        x$vireo_DropletType == "doublet"))
#   temp[i,"unassigned"] <- nrow(x) - temp[i,"doublets"] -temp[i,"singlets"]
#   return(temp)
# })

# results_singlet_doulbet_unnassigned <- lapply(names(results_singlet_doulbet_unnassigned), function(x){
#   results_singlet_doulbet_unnassigned[[x]]$Pool <- x
#   return(results_singlet_doulbet_unnassigned[[x]])
# })
# names(results_singlet_doulbet_unnassigned) <- sample_list

# results_singlet_doulbet_unnassigned_df <- do.call(rbind, results_singlet_doulbet_unnassigned)
# results_singlet_doulbet_unnassigned_df <- pivot_longer(results_singlet_doulbet_unnassigned_df, cols = c("singlets","doublets","unassigned"), names_to = "DropletType", values_to = "Count")

# results_singlet_doulbet_unnassigned_means <- results_singlet_doulbet_unnassigned_df %>%
#   group_by(DropletType, software) %>%
#   summarise( 
#     n=n(),
#     mean=mean(Count),
#     sd=sd(Count)
#   ) %>%
#   mutate( se=sd/sqrt(n))  %>%
#   mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# results_singlet_doulbet_unnassigned_means$software <- factor(results_singlet_doulbet_unnassigned_means$software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","intersection","intersection_no_scSplit"))
# results_singlet_doulbet_unnassigned_df$software <- factor(results_singlet_doulbet_unnassigned_df$software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","intersection","intersection_no_scSplit"))


# pCellNumberPerDropletType <- ggplot(results_singlet_doulbet_unnassigned_df, aes(x = DropletType, y = Count, fill = factor(software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","intersection")))) +
#   geom_bar(alpha = 0.9, stat="identity", position = "dodge") +
#   theme_classic() +
#   facet_wrap(~Pool, ncol = 8, nrow = 10) +
#   labs(fill = "software") +
#   a_primary_fill() +
#   ylim(0, 22000) 
#   # geom_text(data = meta, aes(x = -Inf, y = Inf, label = N, group = DropletType),
#   #           size = 3.5,
#   #           hjust = -0.025,
#   #           vjust = 1.4,
#   #           inherit.aes = FALSE) +
#   # geom_text(data = meta, aes(x = -Inf, y = Inf, label = Genotyped, group = DropletType),
#   #           size = 3.5,
#   #           hjust =-0.025,
#   #           vjust = 3,
#   #           inherit.aes = FALSE)
# ggsave(filename = paste0(out, "CellNumberPerDropletType.png"), plot = pCellNumberPerDropletType, width = 16, height = 9)

# pCellNumberPerDropletType_woscSplit <- ggplot(results_singlet_doulbet_unnassigned_df, aes(x = DropletType, y = Count, fill = software)) +
#   geom_bar(alpha = 0.9, stat="identity", position = "dodge") +
#   theme_classic() +
#   facet_wrap(~Pool, nrow = 10, ncol = 8) +
#   a_primary_fill() +
#   ylim(0, 22000) 
#   # geom_text(data = annotation, aes(x = -Inf, y = Inf, label = N, group = DropletType),
#   #           size = 3.5,
#   #           hjust = -0.025,
#   #           vjust = 1.4,
#   #           inherit.aes = FALSE) +
#   # geom_text(data = annotation, aes(x = -Inf, y = Inf, label = Genotyped, group = DropletType),
#   #           size = 3.5,
#   #           hjust =-0.025,
#   #           vjust = 3,
#   #           inherit.aes = FALSE)
# ggsave(filename = paste0(out, "CellNumberPerDropletType_woscSplit.png"), plot = pCellNumberPerDropletType_woscSplit, width = 16, height = 9)


##### Make a dataframe of singlets #####
## Barcode, genotype ID ##
results_singlets <- lapply(results_list_wide, function(x){
  temp <- x[which(x$demuxlet_DropletType == "singlet" & 
  x$freemuxlet_DropletType == "singlet" &
  x$scSplit_DropletType == "singlet" &
  x$souporcell_DropletType == "singlet" &
  x$vireo_DropletType == "singlet" &
  x$scds_DropletType == "singlet" &
  x$scrublet_DropletType == "singlet" &
  x$DoubletDetection_DropletType == "singlet" &
  x$DoubletFinder_DropletType == "singlet" &
  x$demuxlet_Assignment == x$freemuxlet_Assignment & 
  x$demuxlet_Assignment == x$scSplit_Assignment &
  x$demuxlet_Assignment == x$souporcell_Assignment &
  x$demuxlet_Assignment == x$vireo_Assignment),c("Barcode","demuxlet_Assignment")]
  colnames(temp) <- c("Barcode","Individual")
  return(temp)
})

results_singlets_demultiplex <- lapply(results_list_wide, function(x){
  temp <- x[which(x$demuxlet_DropletType == "singlet" & 
  x$freemuxlet_DropletType == "singlet" &
  x$scSplit_DropletType == "singlet" &
  x$souporcell_DropletType == "singlet" &
  x$vireo_DropletType == "singlet" &
  x$demuxlet_Assignment == x$freemuxlet_Assignment & 
  x$demuxlet_Assignment == x$scSplit_Assignment &
  x$demuxlet_Assignment == x$souporcell_Assignment &
  x$demuxlet_Assignment == x$vireo_Assignment),c("Barcode","demuxlet_Assignment")]
  colnames(temp) <- c("Barcode","Individual")
  return(temp)
})

## Make a list of individuas in list of pools that has just barcode IDs as DF ##
results_singlets_individuals_list <- lapply(results_singlets, function(x){
  temp <- list()
  for (ind in unique(x$Individual)){
    temp[[ind]] <- x[which(x$Individual == ind),"Barcode"]
  }
  return(temp)
})

# lapply(names(results_singlets_individuals_list), function(x){
#   lapply(names(results_singlets_individuals_list[[x]]), function(y){
#     write_delim(results_singlets_individuals_list[[x]][[y]], delim = "\t", path = paste0(out,x,"/",x,"_",y,"_droplet_barcodes.tsv"), col_names = FALSE)
#   })
# })

##### Make a boxplot figure of the number of singlets per individual per pool #####
results_singlets_4df <- lapply(names(results_singlets), function(x){
  results_singlets[[x]] <- as.data.frame(table(results_singlets[[x]][,c("Individual")]))
  colnames(results_singlets[[x]]) <- c("Individual","Count")
  results_singlets[[x]]$Pool <- x
  return(results_singlets[[x]])
})
results_singlets_df <- do.call(rbind,results_singlets_4df)
results_singlets_df <- left_join(results_singlets_df, meta, by = c("Pool"))
results_singlets_df$N <- gsub("Size","",results_singlets_df$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
results_singlets_df$N <- as.numeric(as.character(results_singlets_df$N))
results_singlets_df <- results_singlets_df[order(results_singlets_df$N),]
results_singlets_df$Pool <- gsub("SimulatedPools","",results_singlets_df$Pool)
results_singlets_df$Pool <- factor(results_singlets_df$Pool, levels = unique(results_singlets_df$Pool))
results_singlets_df$Condition <- "All_Softwares"


results_singlets_demultiplex_4df <- lapply(names(results_singlets_demultiplex), function(x){
  results_singlets_demultiplex[[x]] <- as.data.frame(table(results_singlets_demultiplex[[x]][,c("Individual")]))
  colnames(results_singlets_demultiplex[[x]]) <- c("Individual","Count")
  results_singlets_demultiplex[[x]]$Pool <- x
  return(results_singlets_demultiplex[[x]])
})
results_singlets_demultiplex_df <- do.call(rbind,results_singlets_demultiplex_4df)
results_singlets_demultiplex_df <- left_join(results_singlets_demultiplex_df, meta, by = c("Pool"))
results_singlets_demultiplex_df$N <- gsub("Size","",results_singlets_demultiplex_df$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
results_singlets_demultiplex_df$N <- as.numeric(as.character(results_singlets_demultiplex_df$N))
results_singlets_demultiplex_df <- results_singlets_demultiplex_df[order(results_singlets_demultiplex_df$N),]
results_singlets_demultiplex_df$Pool <- gsub("SimulatedPools","",results_singlets_demultiplex_df$Pool)
results_singlets_demultiplex_df$Pool <- factor(results_singlets_demultiplex_df$Pool, levels = unique(results_singlets_demultiplex_df$Pool))
results_singlets_demultiplex_df$Condition <- "Only_Demultiplexing_Softwares"

# pSingletsBox <- ggplot(results_singlets_df, aes(x=factor(Pool), y = Count)) +
#   geom_boxplot() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90))
# ggsave(filename = paste0(out, "SingletsBoxPlot.png"), plot = pSingletsBox, width = 16, height = 9)

##### Make a dataframe of doublets #####
results_doublets <- lapply(results_list_wide, function(x){
  temp <- x[which(x$demuxlet_DropletType == "doublet" &
                                       x$freemuxlet_DropletType == "doublet" &
                                       x$scSplit_DropletType == "doublet" &
                                       x$souporcell_DropletType == "doublet" &
                                       x$vireo_DropletType == "doublet" &
                                       x$scds_DropletType == "doublet" &
                                       x$scrublet_DropletType == "doublet" &
                                       x$DoubletDetection_DropletType == "doublet" &
                                       x$DoubletFinder_DropletType == "doublet"
                                       ),c("Barcode","demuxlet_Assignment")]
  colnames(temp) <- c("Barcode","Individual")
  return(temp)
})

results_doublets_demultiplex <- lapply(results_list_wide, function(x){
  temp <- x[which(x$demuxlet_DropletType == "doublet" &
                                       x$freemuxlet_DropletType == "doublet" &
                                       x$scSplit_DropletType == "doublet" &
                                       x$souporcell_DropletType == "doublet" &
                                       x$vireo_DropletType == "doublet"
                                       ),c("Barcode","demuxlet_Assignment")]
  colnames(temp) <- c("Barcode","Individual")
  return(temp)
})

## Make a list of individuas in list of pools that has just barcode IDs as DF ##
results_doublets_individuals_list <- lapply(results_doublets, function(x){
  temp <- list()
  for (ind in unique(x$Individual)){
    temp[[ind]] <- x[which(x$Individual == ind),"Barcode"]
  }
  return(temp)
})

results_doublets_demultiplex_individuals_list <- lapply(results_doublets_demultiplex, function(x){
  temp <- list()
  for (ind in unique(x$Individual)){
    temp[[ind]] <- x[which(x$Individual == ind),"Barcode"]
  }
  return(temp)
})

# lapply(names(results_doublets_individuals_list), function(x){
#   lapply(names(results_doublets_individuals_list[[x]]), function(y){
#     write_delim(results_doublets_individuals_list[[x]][[y]], delim = "\t", path = paste0(out,x,"/",x,"_",y,"_droplet_barcodes.tsv"), col_names = FALSE)
#   })
# })

##### Make a boxplot figure of the number of singlets per individual per pool #####
results_doublets_4df <- lapply(names(results_doublets), function(x){
  results_doublets[[x]] <- as.data.frame(table(results_doublets[[x]][,c("Individual")]))
  if (nrow(results_doublets[[x]]) == 0){
    results_doublets[[x]] <- data.frame("Individual" = c("doublet"), "Count" = c(0))
  }
  colnames(results_doublets[[x]]) <- c("Individual","Count")
  results_doublets[[x]]$Pool <- x
  return(results_doublets[[x]])
})

results_doublets_df <- do.call(rbind,results_doublets_4df)
results_doublets_df$N <- gsub("Size","",results_doublets_df$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
results_doublets_df$N <- as.numeric(as.character(results_doublets_df$N))
results_doublets_df$Pool <- gsub("SimulatedPools","",results_doublets_df$Pool)
results_doublets_df <- results_doublets_df[order(results_doublets_df$N),]
results_doublets_df$Pool <- factor(results_doublets_df$Pool, levels = unique(results_doublets_df$Pool))
results_doublets_df$Condition <- "All_Softwares"

### With just the demultiplexing Softwares ###
results_doublets_demultiplex_4df <- lapply(names(results_doublets_demultiplex), function(x){
  results_doublets_demultiplex[[x]] <- as.data.frame(table(results_doublets_demultiplex[[x]][,c("Individual")]))
  if (nrow(results_doublets[[x]]) == 0){
    results_doublets[[x]] <- data.frame("Individual" = c("doublet"), "Count" = c(0))
  }
  colnames(results_doublets_demultiplex[[x]]) <- c("Individual","Count")
  results_doublets_demultiplex[[x]]$Pool <- x
  return(results_doublets_demultiplex[[x]])
})
results_doublets_demultiplex_df <- do.call(rbind,results_doublets_demultiplex_4df)
results_doublets_demultiplex_df$N <- gsub("Size","",results_doublets_demultiplex_df$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
results_doublets_demultiplex_df$N <- as.numeric(as.character(results_doublets_demultiplex_df$N))
results_doublets_demultiplex_df$Pool <- gsub("SimulatedPools","",results_doublets_demultiplex_df$Pool)
results_doublets_demultiplex_df <- results_doublets_demultiplex_df[order(results_doublets_demultiplex_df$N),]
results_doublets_demultiplex_df$Pool <- factor(results_doublets_demultiplex_df$Pool, levels = unique(results_doublets_demultiplex_df$Pool))
results_doublets_demultiplex_df$Condition <- "Only_Demultiplexing_Softwares"
### Combine with and without doublet detection dataframes ###
results_doublets_combined_df <- rbind(results_doublets_df,results_doublets_demultiplex_df)
results_doublets_combined_df$N <- as.numeric(as.character(results_doublets_combined_df$N))
results_doublets_combined_df <- results_doublets_combined_df[order(results_doublets_combined_df$N),]
results_doublets_combined_df$Pool <- factor(results_doublets_combined_df$Pool, levels = unique(results_doublets_combined_df$Pool))
results_doublets_combined_df$Condition <- factor(results_doublets_combined_df$Condition, levels=c("Only_Demultiplexing_Softwares","All_Softwares"))

results_singlets_combined_df <- rbind(results_singlets_df,results_singlets_demultiplex_df)
results_singlets_combined_df$N <- as.numeric(as.character(results_singlets_combined_df$N))
results_singlets_combined_df <- results_singlets_combined_df[order(results_singlets_combined_df$N),]
results_singlets_combined_df$Pool <- factor(results_singlets_combined_df$Pool, levels = unique(results_singlets_combined_df$Pool))
results_singlets_combined_df$N_Condition <- paste0(results_singlets_combined_df$N, results_singlets_combined_df$Condition)
results_singlets_combined_df$Condition <- factor(results_singlets_combined_df$Condition, levels=c("Only_Demultiplexing_Softwares","All_Softwares"))

pSingletsBox_doubletBar <- ggplot() +
  geom_bar(data = results_doublets_combined_df, aes(x=Pool, y = Count,  fill =  ""), stat="identity", position = "dodge") +
  scale_fill_manual("Doublets", values = "grey") +
  new_scale("fill") +
  geom_boxplot(data = results_singlets_combined_df, aes(x=Pool, y = Count, fill = factor(N), color = "")) +
  scale_color_manual("Singlets Per Individual", values=1) +
  scale_fill_viridis(option = "plasma",discrete = TRUE, name = "Number Individuals in Pool") +
  # scale_fill_manual(values = c(viridis(option="plasma", n = 5)[1],alpha(viridis(option="plasma", n = 5))[1],
  #                   viridis(option="plasma", n = 5)[2],alpha(viridis(option="plasma", n = 5))[2],
  #                   viridis(option="plasma", n = 5)[3],alpha(viridis(option="plasma", n = 5))[3],
  #                   viridis(option="plasma", n = 5)[4],alpha(viridis(option="plasma", n = 5))[4],
  #                   viridis(option="plasma", n = 5)[5],alpha(viridis(option="plasma", n = 5))[5])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
      axis.title.x = element_blank(),
      text = element_text(size=18),
      legend.position="top", 
      legend.box = "horizontal",
      plot.title = element_text(hjust = 0.5),
      legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggtitle("Number of Singlets and Doublets Common Across All Demultiplexing Softwares") +
  ylab("Number Droplets") +
  facet_wrap( ~ Condition, ncol=1) +
  guides(color = guide_legend(order = 1, nrow = 1),
         fill = guide_legend(order = 2, nrow = 1))


ggsave(filename = paste0(out, "SingletsBoxPlot_DoubletBar_genotyped.png"), plot = pSingletsBox_doubletBar, width = 16, height = 7.5)

# ##### Make a dataframe that has the number of unique inidividuals identified per software per pool #####
# Indiv_identified <- as.data.frame(matrix(nrow = length(results_list_wide), ncol = 6))
# colnames(Indiv_identified) <- c("Pool","demuxlet","freemuxlet","scSplit","souporcell","vireo")
# Indiv_identified$Pool <- names(results_list_wide)

# for (pool in Indiv_identified$Pool){
#   Indiv_identified[which(Indiv_identified$Pool == pool),"demuxlet"] <- length(unique(results_list_wide[[pool]]$demuxlet_Assignment)[!(unique(results_list_wide[[pool]]$demuxlet_Assignment) %in% c("doublet","unsure","unassigned"))])
#   Indiv_identified[which(Indiv_identified$Pool == pool),"freemuxlet"] <- length(unique(results_list_wide[[pool]]$freemuxlet_Assignment)[!(unique(results_list_wide[[pool]]$freemuxlet_Assignment) %in% c("doublet","unsure","unassigned"))])
#   Indiv_identified[which(Indiv_identified$Pool == pool),"scSplit"] <- length(unique(results_list_wide[[pool]]$scSplit_Assignment)[!(unique(results_list_wide[[pool]]$scSplit_Assignment) %in% c("doublet","unsure","unassigned"))])
#   Indiv_identified[which(Indiv_identified$Pool == pool),"souporcell"] <- length(unique(results_list_wide[[pool]]$souporcell_Assignment)[!(unique(results_list_wide[[pool]]$souporcell_Assignment) %in% c("doublet","unsure","unassigned"))])
#   Indiv_identified[which(Indiv_identified$Pool == pool),"vireo"] <- length(unique(results_list_wide[[pool]]$vireo_Assignment)[!(unique(results_list_wide[[pool]]$vireo_Assignment) %in% c("doublet","unsure","unassigned"))][!str_detect(unique(results_list_wide[[pool]]$vireo_Assignment)[!(unique(results_list_wide[[pool]]$vireo_Assignment) %in% c("doublet","unsure","unassigned"))],pattern="donor")])
# }

# Indiv_identified$Genotyped <- gsub("Size","",Indiv_identified$Pool) %>% gsub("_SimulatedPools[0-9]+","",.)
# Indiv_identified$Genotyped <- as.numeric(as.character(Indiv_identified$Genotyped))

# Indiv_identified_long <- pivot_longer(Indiv_identified,cols = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","Genotyped"), names_to = "Software")
# Indiv_identified_long$Software <- factor(Indiv_identified_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo","Genotyped"))
# Indiv_identified_long$Pool  <- gsub("OneK1K_scRNA_Sample","Pool", Indiv_identified_long$Pool)

# pIndividualsIdentified <- ggplot(Indiv_identified_long, aes(x=1, y = value, fill = Software)) +
#   geom_bar(stat="identity", position = "dodge") +
#   theme_classic() +
#   facet_wrap(~Pool, nrow = 4, ncol = 20) +
#   a_primary_fill()

# ggsave(filename = paste0(out, "NumberIndividualsIdentified.png"), plot = pIndividualsIdentified, width = 16, height = 9)



# Indiv_identified_percent <- Indiv_identified
# Indiv_identified_percent$demuxlet_percentage <- Indiv_identified$demuxlet/Indiv_identified$Genotyped
# Indiv_identified_percent$freemuxlet_percentage <- Indiv_identified$freemuxlet/Indiv_identified$Genotyped
# Indiv_identified_percent$scSplit_percentage <- Indiv_identified$scSplit/Indiv_identified$Genotyped
# Indiv_identified_percent$souporcell_percentage <- Indiv_identified$souporcell/Indiv_identified$Genotyped
# Indiv_identified_percent$vireo_percentage <- Indiv_identified$vireo/Indiv_identified$Genotyped
# Indiv_identified_percent$demuxlet <- NULL
# Indiv_identified_percent$freemuxlet <- NULL
# Indiv_identified_percent$scSplit <- NULL
# Indiv_identified_percent$souporcell <- NULL
# Indiv_identified_percent$vireo <- NULL
# # Indiv_identified_percent$Genotyped <- NULL
# Indiv_identified_percent_long <- pivot_longer(Indiv_identified_percent,cols = c("demuxlet_percentage","freemuxlet_percentage","scSplit_percentage","souporcell_percentage","vireo_percentage"), names_to = "Software")
# Indiv_identified_percent_long$Software <- factor(Indiv_identified_percent_long$Software, levels = c("demuxlet_percentage","freemuxlet_percentage","scSplit_percentage","souporcell_percentage","vireo_percentage"))

# pIndividualsIdentified_percnt <- ggplot(Indiv_identified_percent_long, aes(x=factor(Genotyped), y = value*100, fill = gsub("_percentage","",Software))) +
#   geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.9, jitter.width = 0.12), shape = 21, color = "black", size = 5, stroke = 1.25) +
#   theme_classic() +
#   a_primary_fill(name = "Software") +
#   theme(text = element_text(size=30)) +
#   ylab("Percent of Individuals Identified") +
#   xlab("Number of Individuals in Pool")

# ggsave(filename = paste0(out, "NumberIndividualsIdentified_percent_violin.png"), plot = pIndividualsIdentified_percnt, width = 16, height = 9)

# pIndividualsIdentified_percnt_heat <- ggplot(Indiv_identified_percent_long, aes(y=gsub("_percentage","",Software), x= gsub("_SimulatedPools",".",Pool), fill = value*100)) +
#   geom_tile(color = "white") +
#   theme_classic() +
#   theme(text = element_text(size=20),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     legend.position="bottom", 
#     legend.box = "horizontal") +
#   ylab("Software") +
#   xlab("Simulated Pool") +
#   scale_fill_viridis(direction = -1, name = "Percent Individuals Identified")

# ggsave(filename = paste0(out, "NumberIndividualsIdentified_percent_heatmap.png"), plot = pIndividualsIdentified_percnt_heat, width = 20, height = 3.5)


# Indiv_missed <- Indiv_identified
# Indiv_missed$demuxlet_percentage <- Indiv_identified$Genotyped- Indiv_identified$demuxlet
# Indiv_missed$freemuxlet_percentage <- Indiv_identified$Genotyped- Indiv_identified$freemuxlet
# Indiv_missed$scSplit_percentage <- Indiv_identified$Genotyped- Indiv_identified$scSplit
# Indiv_missed$souporcell_percentage <- Indiv_identified$Genotyped- Indiv_identified$souporcell
# Indiv_missed$vireo_percentage <- Indiv_identified$Genotyped- Indiv_identified$vireo
# Indiv_missed$TotalMissing <- Indiv_missed$demuxlet_percentage + Indiv_missed$freemuxlet_percentage + Indiv_missed$scSplit_percentage + Indiv_missed$souporcell_percentage + Indiv_missed$vireo_percentage
# Indiv_missed$demuxlet <- NULL
# Indiv_missed$freemuxlet <- NULL
# Indiv_missed$scSplit <- NULL
# Indiv_missed$souporcell <- NULL
# Indiv_missed$vireo <- NULL
# colnames(Indiv_missed) <- gsub("_percentage","",colnames(Indiv_missed))
# Indiv_missed <- Indiv_missed[order( Indiv_missed[,"Genotyped"], Indiv_missed[,"TotalMissing"] ),]

# ht_opt(legend_title_position = "topleft", 
#   legend_labels_gp = gpar(fontsize = 18),
#   legend_title_gp = gpar(fontsize = 18, fontface = "bold"),
#   heatmap_row_names_gp = gpar(fontsize = 14),
#   heatmap_column_title_gp = gpar(fontsize = 18))
# colors <- list("Number Individuals in Pool" = plasma(length(unique(Indiv_missed$Genotyped))))
# colors[["Number Individuals in Pool"]] <- factor(colors[["Number Individuals in Pool"]], levels = colors[["Number Individuals in Pool"]])
# names(colors[["Number Individuals in Pool"]]) <- unique(Indiv_missed$Genotyped)
# column_ha = HeatmapAnnotation("Number Individuals in Pool" = Indiv_missed$Genotyped, col = colors,
#   annotation_name_side = "left",
#   annotation_legend_param = list(
#         "Number Individuals in Pool" = list(
#                 title = "Number Individuals in Pool",
#                 at = names(colors[["Number Individuals in Pool"]]),
#                 nrow = 1
#             )))

# matrix <- as.matrix(t(Indiv_missed[,c("demuxlet","freemuxlet","scSplit","souporcell","vireo")]))
# colnames(matrix) <- NULL
# png(paste0(out, "NumberIndividualsMissed_heatmap.png"), width = 16, height = 2.25, units = "in", res = 300)
# draw(Heatmap(matrix, name = "Number Individuals Missed", 
#   col=viridis(2),
#   row_order = c("demuxlet","freemuxlet","scSplit","souporcell","vireo"),
#   column_order = 1:ncol(matrix),
#   bottom_annotation = column_ha,
#   column_split = Indiv_missed$Genotyped,
#   rect_gp = gpar(col = "grey", lwd = 1),
#   column_title = "Pools", 
#   column_title_side = "bottom",
#   row_names_side = "left",
#   heatmap_legend_param = list(nrow = 1)
#   ), merge_legend = TRUE,annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
# dev.off()



# # Indiv_missed$Genotyped <- NULL
# Indiv_missed_long <- pivot_longer(Indiv_missed,cols = c("demuxlet_percentage","freemuxlet_percentage","scSplit_percentage","souporcell_percentage","vireo_percentage"), names_to = "Software")
# Indiv_missed_long$Software <- gsub("_percentage","",Indiv_missed_long$Software)
# Indiv_missed_long$Software <- factor(Indiv_missed_long$Software, levels = c("demuxlet","freemuxlet","scSplit","souporcell","vireo"))
# Indiv_missed_long$Size <- gsub("_SimulatedPools[0-9]+","",Indiv_missed_long$Pool) %>% gsub("Size","",.)
# Indiv_missed_long$Size <- as.numeric(Indiv_missed_long$Size)
# Indiv_missed_long <- Indiv_missed_long[order( Indiv_missed_long[,"Size"], Indiv_missed_long[,"value"] ),]
# Indiv_missed_long$Pool <- factor(Indiv_missed_long$Pool, levels = unique(Indiv_missed_long$Pool))

# pIndividualsMissed_heat <- ggplot(Indiv_missed_long, aes(y=Software, x= factor(gsub("_SimulatedPools",".",Pool), levels=unique(gsub("_SimulatedPools",".",Pool))), fill = value)) +
#   geom_tile(color = "white") +
#   theme_classic() +
#   theme(text = element_text(size=20),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     legend.position="bottom", 
#     legend.box = "horizontal") +
#   ylab("Software") +
#   xlab("Simulated Pool") +
#   scale_fill_viridis( name = "Percent Individuals Identified")

# ggsave(filename = paste0(out, "NumberIndividualsMissed_heatmap.png"), plot = pIndividualsMissed_heat, width = 20, height = 3.5)

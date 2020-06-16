library(tidyr)
library(tidyverse)
library(dplyr)
library(vcfR)
library(lsa)
library(ComplexHeatmap)

########## Set up directories ##########
dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/SimulatedPools/"
out <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/SimulatedPools/SimulatedOverlap/"

########## Set up variables ##########
pools <- dir(path = dir, pattern = "Size")

########## Set up functions ##########
calculate_DS <- function(GP_df){
    columns <- c()
    for (i in 1:ncol(GP_df)){
        columns <- c(columns, paste0(colnames(GP_df)[i],"-0"), paste0(colnames(GP_df)[i],"-1"), paste0(colnames(GP_df)[i],"-2"))
    }
    df <- GP_df
    colnames(df) <- paste0("c", colnames(df))
    colnames_orig <- colnames(df)
    for (i in 1:length(colnames_orig)){
        df <- separate(df, sep = ",", col = colnames_orig[i], into = columns[(1+(3*(i-1))):(3+(3*(i-1)))])
    }
    df <- mutate_all(df, function(x) as.numeric(as.character(x)))
    for (i in 1: ncol(GP_df)){
        GP_df[,i] <- df[,(2+((i-1)*3))] + 2* df[,(3+((i-1)*3))]
    }
    return(GP_df)
}

pearson_correlation <- function(df, ref_df, clust_df){
    for (col in colnames(df)){
        for (row in rownames(df)){
            df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson")
        }
    }
    return(df)
}


########## Read in vcf files for each of three non-reference genotype softwares ##########
ref_geno_list <- lapply(pools, function(x){
    temp <- list()
    temp[["freemuxlet"]] <- read.vcfR(paste0(dir,x,"/popscle/freemuxlet/Individual_genotypes_subset.vcf.gz"))
    temp[["scSplit"]] <- read.vcfR(paste0(dir,x,"/scSplit/Individual_genotypes_subset.vcf.gz"))
    temp[["souporcell"]] <- read.vcfR(paste0(dir,x,"/souporcell/Individual_genotypes_subset.vcf.gz"))
    return(temp)
})
names(ref_geno_list) <- pools

cluster_geno_list <- lapply(pools, function(x){
    temp <- list()
    temp[["freemuxlet"]] <- read.vcfR(paste0(dir,x,"/popscle/freemuxlet/freemuxletOUT.clust1.vcf.gz"))
    temp[["scSplit"]] <- read.vcfR(paste0(dir,x,"/scSplit/scSplit.vcf"))
    temp[["souporcell"]] <- read.vcfR(paste0(dir,x,"/souporcell/cluster_genotypes.vcf"))
    return(temp)
})
names(cluster_geno_list) <- pools

########## Convert to tidy data frame ##########
ref_geno_tidy_list <- lapply(ref_geno_list, function(x){
    temp <- list()
    temp[["freemuxlet"]] <- as_tibble(extract.gt(element = "DS",x[["freemuxlet"]], IDtoRowNames =F))
    temp[["freemuxlet"]]$ID <- paste0(x[["freemuxlet"]]@fix[,'CHROM'],":", x[["freemuxlet"]]@fix[,'POS'],"_", x[["freemuxlet"]]@fix[,'REF'], "_",x[["freemuxlet"]]@fix[,'ALT'])
    temp[["scSplit"]] <- as_tibble(extract.gt(element = "GP", x[["scSplit"]], IDtoRowNames =F))
    temp[["scSplit"]] <- calculate_DS(temp[["scSplit"]])
    temp[["scSplit"]]$ID <- paste0(x[["scSplit"]]@fix[,'CHROM'],":", x[["scSplit"]]@fix[,'POS'])
    temp[["scSplit"]] <- temp[["scSplit"]][!(temp[["scSplit"]]$ID %in% temp[["scSplit"]]$ID[duplicated(temp[["scSplit"]]$ID)]),]
    temp[["souporcell"]] <- as_tibble(extract.gt(element = "DS",x[["souporcell"]], IDtoRowNames =F))
    temp[["souporcell"]]$ID <- paste0(x[["souporcell"]]@fix[,'CHROM'],":", x[["souporcell"]]@fix[,'POS'],"_", x[["souporcell"]]@fix[,'REF'], "_",x[["souporcell"]]@fix[,'ALT'])
    return(temp)
})

cluster_geno_tidy_list <- lapply(cluster_geno_list, function(x){
    temp <- list()
    temp[["freemuxlet"]] <- as_tibble(extract.gt(element = "GP",x[["freemuxlet"]], IDtoRowNames =F))
    temp[["freemuxlet"]] <- calculate_DS(temp[["freemuxlet"]])
    temp[["freemuxlet"]]$ID <- paste0(x[["freemuxlet"]]@fix[,'CHROM'],":", x[["freemuxlet"]]@fix[,'POS'],"_", x[["freemuxlet"]]@fix[,'REF'], "_",x[["freemuxlet"]]@fix[,'ALT'])
    temp[["freemuxlet"]] <- temp[["freemuxlet"]][colSums(!is.na(temp[["freemuxlet"]])) > 0]
    temp[["freemuxlet"]] <- temp[["freemuxlet"]][complete.cases(temp[["freemuxlet"]]),]
    
    temp[["scSplit"]] <- as_tibble(extract.gt(element = "GP",x[["scSplit"]], IDtoRowNames =F))
    temp[["scSplit"]] <- calculate_DS(temp[["scSplit"]])
    temp[["scSplit"]]$ID <- paste0(x[["scSplit"]]@fix[,'CHROM'],":", x[["scSplit"]]@fix[,'POS'])
    temp[["scSplit"]] <- temp[["scSplit"]][colSums(!is.na(temp[["scSplit"]])) > 0]
    temp[["scSplit"]] <- temp[["scSplit"]][complete.cases(temp[["scSplit"]]),]

    temp[["souporcell"]] <- as_tibble(extract.gt(element = "GT",x[["souporcell"]], IDtoRowNames =F))
    temp[["souporcell"]] <- as_tibble(lapply(temp[["souporcell"]], function(x) {gsub("0/0",0, x)}) %>%
                                lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                lapply(., function(x) {gsub("1/1",2, x)}))
    temp[["souporcell"]]$ID <- paste0(x[["souporcell"]]@fix[,'CHROM'],":", x[["souporcell"]]@fix[,'POS'],"_", x[["souporcell"]]@fix[,'REF'], "_",x[["souporcell"]]@fix[,'ALT'])
    temp[["souporcell"]] <- temp[["souporcell"]][colSums(!is.na(temp[["souporcell"]])) > 0]
    temp[["souporcell"]] <- temp[["souporcell"]][complete.cases(temp[["souporcell"]]),]
    return(temp)
})

########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
locations_list <- lapply(names(ref_geno_list), function(x){
    temp <- list()
    temp_freemuxlet <- inner_join(ref_geno_tidy_list[[x]][["freemuxlet"]][,"ID"],cluster_geno_tidy_list[[x]][["freemuxlet"]][,"ID"])
    temp[["freemuxlet"]] <- temp_freemuxlet[!(temp_freemuxlet$ID %in% temp_freemuxlet[duplicated(temp_freemuxlet),"ID"]),]
    temp_scSplit <- inner_join(ref_geno_tidy_list[[x]][["scSplit"]][,"ID"],cluster_geno_tidy_list[[x]][["scSplit"]][,"ID"])
    temp[["scSplit"]] <- temp_scSplit[!(temp_scSplit$ID %in% temp_scSplit[duplicated(temp_scSplit),"ID"]),]
    temp_souporcell <- inner_join(ref_geno_tidy_list[[x]][["souporcell"]][,"ID"],cluster_geno_tidy_list[[x]][["souporcell"]][,"ID"])
    temp[["souporcell"]] <- temp_souporcell[!(temp_souporcell$ID %in% temp_souporcell[duplicated(temp_souporcell),"ID"]),]
    return(temp)
})
names(locations_list) <- pools


########## Keep just the SNPs that overlap ##########
ref_geno_tidy_list <- lapply(names(ref_geno_tidy_list), function(x){
    temp <- list()
    temp[["freemuxlet"]] <- left_join(locations_list[[x]][["freemuxlet"]], ref_geno_tidy_list[[x]][["freemuxlet"]])
    temp[["scSplit"]] <- left_join(locations_list[[x]][["scSplit"]], ref_geno_tidy_list[[x]][["scSplit"]])
    temp[["souporcell"]] <- left_join(locations_list[[x]][["souporcell"]], ref_geno_tidy_list[[x]][["souporcell"]])
    return(temp)
})
names(ref_geno_tidy_list) <- pools

cluster_geno_tidy_list <- lapply(names(cluster_geno_tidy_list), function(x){
    temp <- list()
    temp[["freemuxlet"]] <- left_join(locations_list[[x]][["freemuxlet"]], cluster_geno_tidy_list[[x]][["freemuxlet"]])
    temp[["scSplit"]] <- left_join(locations_list[[x]][["scSplit"]], cluster_geno_tidy_list[[x]][["scSplit"]])
    temp[["souporcell"]] <- left_join(locations_list[[x]][["souporcell"]], cluster_geno_tidy_list[[x]][["souporcell"]])
    return(temp)
})
names(cluster_geno_tidy_list) <- pools


########## Correlate all the cluster genotypes with the individuals genotyped ##########
##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
pearson_correlations <- lapply(names(ref_geno_tidy_list), function(x){
    temp <- list()
    temp[["freemuxlet"]] <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy_list[[x]][["freemuxlet"]]) -1), ncol = (ncol(ref_geno_tidy_list[[x]][["freemuxlet"]]) -1)))
    colnames(temp[["freemuxlet"]]) <- colnames(ref_geno_tidy_list[[x]][["freemuxlet"]])[2:(ncol(ref_geno_tidy_list[[x]][["freemuxlet"]]))]
    rownames(temp[["freemuxlet"]]) <- colnames(cluster_geno_tidy_list[[x]][["freemuxlet"]])[2:(ncol(cluster_geno_tidy_list[[x]][["freemuxlet"]]))]
    temp[["freemuxlet"]] <- pearson_correlation(temp[["freemuxlet"]], ref_geno_tidy_list[[x]][["freemuxlet"]], cluster_geno_tidy_list[[x]][["freemuxlet"]])

    temp[["scSplit"]] <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy_list[[x]][["scSplit"]]) -1), ncol = (ncol(ref_geno_tidy_list[[x]][["scSplit"]]) -1)))
    colnames(temp[["scSplit"]]) <- colnames(ref_geno_tidy_list[[x]][["scSplit"]])[2:(ncol(ref_geno_tidy_list[[x]][["scSplit"]]))]
    rownames(temp[["scSplit"]]) <- colnames(cluster_geno_tidy_list[[x]][["scSplit"]])[2:(ncol(cluster_geno_tidy_list[[x]][["scSplit"]]))]
    temp[["scSplit"]] <- pearson_correlation(temp[["scSplit"]], ref_geno_tidy_list[[x]][["scSplit"]], cluster_geno_tidy_list[[x]][["scSplit"]])

    temp[["souporcell"]] <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy_list[[x]][["souporcell"]]) -1), ncol = (ncol(ref_geno_tidy_list[[x]][["souporcell"]]) -1)))
    colnames(temp[["souporcell"]]) <- colnames(ref_geno_tidy_list[[x]][["souporcell"]])[2:(ncol(ref_geno_tidy_list[[x]][["souporcell"]]))]
    rownames(temp[["souporcell"]]) <- colnames(cluster_geno_tidy_list[[x]][["souporcell"]])[2:(ncol(cluster_geno_tidy_list[[x]][["souporcell"]]))]
    temp[["souporcell"]] <- pearson_correlation(temp[["souporcell"]], ref_geno_tidy_list[[x]][["souporcell"]], cluster_geno_tidy_list[[x]][["souporcell"]])

    return(temp)
})
names(pearson_correlations) <- pools

########## Save the correlation dataframes ##########
saveRDS(pearson_correlations, file= paste0(out,"reference_cluster_genotype_pearson_correlations.rds"))
lapply(names(pearson_correlations), function(x){
    dir.create(paste0(out,x))
    lapply(names(pearson_correlations[[x]]), function(y){
        write_delim(pearson_correlations[[x]][[y]],path = paste0(out,x,"/",x,"_",y,"_pearson_correlations.tsv"), delim = "\t" )
    })
})

pearson_correlations <- readRDS(file= paste0(out,"reference_cluster_genotype_pearson_correlations.rds"))


########## Create correlation figures ##########
col_fun = colorRampPalette(c("white", "red"))
pPearsonCorrelations_list <- lapply(names(pearson_correlations), FUN = function(x){
    lapply(names(pearson_correlations[[x]]), function(y){
        Heatmap(as.matrix(pearson_correlations[[x]][[y]]), cluster_rows = T, col = col_fun(100), column_title = paste0(x,"_",y))
    })
}) 
names(pPearsonCorrelations_list) <- pools
pPearsonCorrelations_list <- lapply(pPearsonCorrelations_list, function(x){
    names(x) <- c("freemuxlet","scSplit","souporcell")
    return(x)
})


########## Save the correlation figures ##########
lapply(names(pPearsonCorrelations_list), FUN = function(x){
    lapply(names(pPearsonCorrelations_list[[x]]), function(y){
        png(filename = paste0(out,x,"/",x,"_",y,"_pearson_correlation.png"), width = 500)
        print(pPearsonCorrelations_list[[x]][[y]])
    dev.off()
    })
})


########## Assign individual to cluster based on highest correlating individual ##########
key_list_no_filtering <- lapply(pearson_correlations, function(x){
    lapply(x, function(y){
        df <- as.data.frame(matrix(nrow = ncol(y), ncol = 3))
        colnames(df) <- c("Genotype_ID","Cluster_ID","Correlation")
        df$Genotype_ID <- colnames(y)
        for (id in df$Genotype_ID){
                df$Cluster_ID[which(df$Genotype_ID == id)] <- rownames(y)[which.max(y[,id])]
                df$Correlation[which(df$Genotype_ID == id)] <- max(y[,id])
        }
        return(df)
    })
})

key_list <- lapply(pearson_correlations, function(x){
    lapply(x, function(y){
        df <- as.data.frame(matrix(nrow = ncol(y), ncol = 3))
        colnames(df) <- c("Genotype_ID","Cluster_ID","Correlation")
        df$Genotype_ID <- colnames(y)
        for (id in df$Genotype_ID){
            if (max(y[,id]) == max(y[rownames(y)[which.max(y[,id])],])){
                df$Cluster_ID[which(df$Genotype_ID == id)] <- rownames(y)[which.max(y[,id])]
                df$Correlation[which(df$Genotype_ID == id)] <- max(y[,id])
            } else {
                df <- df[-c(which(df$Genotype_ID == id)),]
            }
        }
        return(df)
    })
})

key_list_no_filtering <- lapply(key_list_no_filtering, function(x){
    lapply(names(x), function(y){
        x[[y]]$Software <- y
        return(x[[y]])
    })
})

key_list <- lapply(key_list, function(x){
    lapply(names(x), function(y){
        x[[y]]$Software <- y
        return(x[[y]])
    })
})

key_list_no_filtering <- lapply(key_list_no_filtering, function(x){
    names(x) <- c("freemuxlet","scSplit","souporcell")
    return(x)
})

key_list <- lapply(key_list, function(x){
    names(x) <- c("freemuxlet","scSplit","souporcell")
    return(x)
})

key_no_filtering <- lapply(key_list_no_filtering, function(x){
    do.call(rbind, x)
})

key <- lapply(key_list, function(x){
    do.call(rbind, x)
})

saveRDS(key_list_no_filtering, paste0(out,"PoolKeys_preFiltering.rds"))
saveRDS(key, paste0(out,"PoolKeys.rds"))



########## Check that all assignments are unique ##########
unique_assignments_list <- lapply(names(key_list), function(x){
    lapply(names(key_list[[x]]), function(y){
        unique_assignments <- as.data.frame(matrix(ncol=3, nrow = 1))
        colnames(unique_assignments) <- c("Pool","Software","AllClusterIDsUnique")
        unique_assignments$Pool <- x
        unique_assignments$Software <- y
        if(all(!duplicated(key_list[[x]][[y]]$Cluster_ID))){
            unique_assignments$AllClusterIDsUnique <- TRUE
            # print(paste0("All the cluster IDs are unique for pool ",x," and software ",y))
        } else{
            unique_assignments$AllClusterIDsUnique <- FALSE
            # print(paste0("Not all the cluster IDs are unique. You will have to manually check for pool ", x, " for software ", y))
        }
        return(unique_assignments)
    })
})

unique_assignments_list <- lapply(unique_assignments_list, function(x){
    do.call(rbind, x)
})

unique_assignments <- do.call(rbind, unique_assignments_list)
write_delim(unique_assignments, path = paste0(out,"UniqueClusterAssignmentCheck.tsv"), delim = "\t" )

########## Make a list of those that would have failed pre-filtering
unique_assignments_prefilter_list <- lapply(names(key_list_no_filtering), function(x){
    lapply(names(key_list_no_filtering[[x]]), function(y){
        unique_assignments <- as.data.frame(matrix(ncol=3, nrow = 1))
        colnames(unique_assignments) <- c("Pool","Software","AllClusterIDsUnique")
        unique_assignments$Pool <- x
        unique_assignments$Software <- y
        if(all(!duplicated(key_list_no_filtering[[x]][[y]]$Cluster_ID))){
            unique_assignments$AllClusterIDsUnique <- TRUE
            # print(paste0("All the cluster IDs are unique for pool ",x," and software ",y))
        } else{
            unique_assignments$AllClusterIDsUnique <- FALSE
            # print(paste0("Not all the cluster IDs are unique. You will have to manually check for pool ", x, " for software ", y))
        }
        return(unique_assignments)
    })
})

unique_assignments_prefilter_list <- lapply(unique_assignments_prefilter_list, function(x){
    do.call(rbind, x)
})

unique_assignments_prefilter <- do.call(rbind, unique_assignments_prefilter_list)
write_delim(unique_assignments_prefilter, path = paste0(out,"UniqueClusterAssignmentCheck_prefiltering.tsv"), delim = "\t" )

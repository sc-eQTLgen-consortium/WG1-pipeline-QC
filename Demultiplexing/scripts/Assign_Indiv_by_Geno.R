#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/Assign_Indiv_by_Geno.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-b", "--basedir", required = TRUE, help="")
parser$add_argument("-p", "--pool", required = TRUE, type = "character", help="")
parser$add_argument("-r", "--result_file", required = TRUE, type = "character", help="")
parser$add_argument("-c", "--correlation_limit", required = FALSE, default = 0.7, type = "double", help="The minimum correlation between a cluster and provided SNP genotypes to consider that cluster assigned to that individual.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()


suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(vcfR)))
suppressMessages(suppressWarnings(library(lsa)))
suppressMessages(suppressWarnings(library(ComplexHeatmap)))

## make sure the directory exists ###
dir.create(paste0(args$basedir, "/", args$pool, "/souporcell/genotype_correlations/"))


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
            df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson", use = "complete.obs")
        }
    }
    return(df)
}


########## Read in vcf files for each of three non-reference genotype softwares ##########
ref_geno <- read.vcfR(paste0(args$basedir, "/", args$pool, "/souporcell/Individual_genotypes_subset.vcf.gz"))
cluster_geno <- read.vcfR(paste0(args$basedir, "/", args$pool, "/souporcell/cluster_genotypes.vcf"))


########## Convert to tidy data frame ##########
ref_geno_tidy <- as_tibble(extract.gt(element = "DS",ref_geno, IDtoRowNames =F))
ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'], "_",ref_geno@fix[,'ALT'])
ref_geno_tidy <- ref_geno_tidy[!(ref_geno_tidy$ID %in% ref_geno_tidy$ID[duplicated(ref_geno_tidy$ID)]),]

cluster_geno_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno, IDtoRowNames =F))
cluster_geno_tidy <- as_tibble(lapply(cluster_geno_tidy, function(x) {gsub("0/0",0, x)}) %>%
                                lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                lapply(., function(x) {gsub("1/1",2, x)}))
cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'], "_",cluster_geno@fix[,'ALT'])
cluster_geno_tidy <- cluster_geno_tidy[colSums(!is.na(cluster_geno_tidy)) > 0]
# cluster_geno_tidy <- cluster_geno_tidy[complete.cases(cluster_geno_tidy),]
cluster_geno_tidy <- cluster_geno_tidy[!(cluster_geno_tidy$ID %in% cluster_geno_tidy$ID[duplicated(cluster_geno_tidy$ID)]),]


########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
locations  <- inner_join(ref_geno_tidy[,"ID"],cluster_geno_tidy[,"ID"])
locations <- locations[!(locations$ID %in% locations[duplicated(locations),"ID"]),]

########## Keep just the SNPs that overlap ##########
ref_geno_tidy <- left_join(locations, ref_geno_tidy)
cluster_geno_tidy <- left_join(locations, cluster_geno_tidy)

########## Correlate all the cluster genotypes with the individuals genotyped ##########
##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
pearson_correlations <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy) -1), ncol = (ncol(ref_geno_tidy) -1)))
colnames(pearson_correlations) <- colnames(ref_geno_tidy)[2:(ncol(ref_geno_tidy))]
rownames(pearson_correlations) <- colnames(cluster_geno_tidy)[2:(ncol(cluster_geno_tidy))]
pearson_correlations <- pearson_correlation(pearson_correlations, ref_geno_tidy, cluster_geno_tidy)

########## Save the correlation dataframes ##########
write_delim(pearson_correlations, path = paste0(args$basedir, "/", args$pool, "/souporcell/genotype_correlations/pearson_correlations.tsv"), delim = "\t" )


########## Create correlation figures ##########
col_fun = colorRampPalette(c("white", "red"))(101)
pPearsonCorrelations <- Heatmap(as.matrix(pearson_correlations), cluster_rows = T, col = col_fun, column_title = args$pool)

########## Save the correlation figures ##########
png(filename = paste0(args$basedir, "/", args$pool, "/souporcell/genotype_correlations/pearson_correlation.png"), width = 500)
print(pPearsonCorrelations)
dev.off()

########## Assign individual to cluster based on highest correlating individual ##########
key <- as.data.frame(matrix(nrow = ncol(pearson_correlations), ncol = 3))
colnames(key) <- c("Genotype_ID","Cluster_ID","Correlation")
key$Genotype_ID <- colnames(pearson_correlations)
for (id in key$Genotype_ID){
    if (max(pearson_correlations[,id]) == max(pearson_correlations[rownames(pearson_correlations)[which.max(pearson_correlations[,id])],]) & max(pearson_correlations[,id]) > args$correlation_limit){
        key$Cluster_ID[which(key$Genotype_ID == id)] <- rownames(pearson_correlations)[which.max(pearson_correlations[,id])]
        key$Correlation[which(key$Genotype_ID == id)] <- max(pearson_correlations[,id])
    } else {
        key$Cluster_ID[which(key$Genotype_ID == id)] <- "unassigned"
        key$Correlation[which(key$Genotype_ID == id)] <- NA
    }
}

write_delim(key, path =paste0(args$basedir, "/", args$pool, "/souporcell/genotype_correlations/Genotype_ID_key.txt"), delim = "\t")


##### Read in Files #####
print(as.character(args$result_file))
result <- read_delim(file = as.character(args$result_file), delim = "\t")

##### Left_join the common assignments to the dataframe #####
col_order <- colnames(result)
result <- left_join(result,key, by = c("souporcell_Assignment" = "Cluster_ID"))
result$Genotype_ID <- ifelse(result$souporcell_DropletType == "doublet", "doublet", result$Genotype_ID)
result$Genotype_ID <- ifelse(is.na(result$Genotype_ID), "unassigned", result$Genotype_ID)
result$souporcell_Assignment <- NULL
colnames(result) <- gsub("Genotype_ID", "souporcell_Assignment", colnames(result))
result$Correlation <- NULL
result <- result[,col_order]


write_delim(result, path =paste0(args$basedir, "/", args$pool, "/CombinedResults/CombinedDropletAssignments_w_genotypeIDs.tsv"), delim = "\t")

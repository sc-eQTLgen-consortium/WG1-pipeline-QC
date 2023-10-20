#!/usr/bin/env Rscript
# Author: drneavin
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-s", "--scores", required=TRUE, type="character", help="")
parser$add_argument("-os", "--onekg_scores", required=FALSE, type="character", help="")
parser$add_argument("-p", "--psam", required=FALSE, type="character", help="")
parser$add_argument("-op", "--onekg_psam", required=FALSE, type="character", help="")
parser$add_argument("-sc", "--sex_check", required=FALSE, type="character", help="")
parser$add_argument("-o", "--out", required=TRUE, type="character", help="The output directory where results will be saved.")

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(RColorBrewer)))

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

##### Read in Results #####
data_score <- read_delim(as.character(args$scores), delim = "\t", col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "Provided_Ancestry" = "c"))
onekg_score <- read_delim(as.character(args$onekg_scores), delim = "\t", col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "#IID" = "c", "SuperPop" = "c", "Population" = "c"))
data_anc <- read_delim(as.character(args$psam), delim = "\t", col_types = cols(.default = "c", "SEX" = "d"))
onekg_anc <- read_delim(as.character(args$onekg_psam), delim = "\t", col_types = cols(.default = "c", "SEX" = "d"))
sex_check <- read_delim(as.character(args$sex_check), delim = "\t", col_types = cols(.default = "c", "PEDSEX" = "d", "SNPSEX" = "d", "F" = "d"))


##### Set up variables #####


##### Read in PCA Results #####

##### Remove # from colnames #####
colnames(data_score) <- gsub("#", "", colnames(data_score))
colnames(onekg_score) <- gsub("#", "", colnames(onekg_score))
colnames(onekg_anc) <- gsub("#", "", colnames(onekg_anc))
colnames(data_anc) <- gsub("#", "", colnames(data_anc))
colnames(sex_check) <- gsub("#", "", colnames(sex_check))

print("data_score")
print(head(as.data.frame(data_score)))
print("onekg_score")
print(head(as.data.frame(onekg_score)))
print("onekg_anc")
print(head(as.data.frame(onekg_anc)))
print("data_anc")
print(head(as.data.frame(data_anc)))
print("sex_check")
print(head(as.data.frame(sex_check)))

##### Combine PCA Results Together #####
if (any(grepl("FID", colnames(onekg_score)))){
  onekg_score_temp <- onekg_score[, c(colnames(onekg_score)[grep("FID", colnames(onekg_score))], colnames(onekg_score)[grep("IID", colnames(onekg_score))], colnames(onekg_score)[grep("PC", colnames(onekg_score))])]
  colnames(onekg_score_temp) <- c("FID", "IID", paste0("PC", 1:(ncol(onekg_score_temp)-2)))
} else {
  onekg_score_temp <- onekg_score[, c("IID", colnames(onekg_score)[grep("PC", colnames(onekg_score))])]
  colnames(onekg_score_temp) <- c("IID", paste0("PC", 1:(ncol(onekg_score_temp)-1)))
  onekg_score_temp$FID <- NA
  onekg_score_temp <- onekg_score_temp[, c("FID", "IID", paste0("PC", 1:(ncol(onekg_score_temp)-2)))]
}
print(head(as.data.frame(onekg_score_temp)))

if (any(grepl("FID", colnames(data_score)))){
  data_score_temp <- data_score[, c(colnames(data_score)[grep("FID", colnames(data_score))], colnames(data_score)[grep("IID", colnames(data_score))], colnames(data_score)[grep("PC", colnames(data_score))])]
  colnames(data_score_temp) <- c("FID", "IID", paste0("PC", 1:(ncol(data_score_temp)-2)))
} else {
  data_score_temp <- data_score[, c(colnames(data_score)[grep("IID", colnames(data_score))], colnames(data_score)[grep("PC", colnames(data_score))])]
  colnames(data_score_temp) <- c("IID", paste0("PC", 1:(ncol(data_score_temp)-1)))
  data_score_temp$FID <- NA
}
print(head(as.data.frame(data_score_temp)))


scores <- rbind(onekg_score_temp, data_score_temp)
print(head(as.data.frame(scores)))

scores$IID <- as.character(scores$IID)
onekg_anc$IID <- as.character(onekg_anc$IID)
data_anc$IID <- as.character(data_anc$IID)
scores <- left_join(scores, onekg_anc[, c("IID", "SuperPop")])
scores <- left_join(scores, data_anc[, c("IID", "Provided_Ancestry")])
print(head(as.data.frame(scores)))


##### Calculate Medoids and Assign Clusters #####
pam_res <- pam(scores[, 2:11], 6)


##### Assign Ancestries to Individuals #####
scores$Cluster <- factor(pam_res$clustering)
print(head(as.data.frame(scores)))

conversion_table <- table(scores$SuperPop, scores$Cluster)
conversion_key <- data.frame(Cluster = colnames(conversion_table), Assignment = NA)

for (clust in conversion_key$Cluster){
  conversion_key$Assignment[which(conversion_key$Cluster == clust)] <- rownames(conversion_table)[which.max(conversion_table[, clust])]
}
print(conversion_key)


scores <- left_join(scores, conversion_key)

scores$combined_assignment <- ifelse(is.na(scores$SuperPop), scores$Assignment, scores$SuperPop)
scores$combined_assignment <- ifelse(is.na(scores$SuperPop), scores$Assignment, scores$SuperPop)
scores$Changed <- ifelse(is.na(scores$Provided_Ancestry), "Matched", ifelse(scores$Provided_Ancestry == scores$Assignment, "Matched", paste0("Unmatched-", scores$Provided_Ancestry, "->", scores$Assignment)))

##### Plot Results #####
df4plots <- rbind(data.frame(scores[which(is.na(scores$Provided_Ancestry)), ], Plot = "1000G Reference"), data.frame(scores[which(!is.na(scores$Provided_Ancestry)), ], Plot = "Projected Data Assignments"), data.frame(scores[which(!is.na(scores$Provided_Ancestry)), ], Plot = "Projected Data Assignments vs Original Assignments"))
df4plots$Final_Assignment <- ifelse(df4plots$Plot == "Projected Data Assignments", df4plots$Assignment, ifelse(df4plots$Plot == "Projected Data Assignments vs Original Assignments", df4plots$Changed, df4plots$combined_assignment))
df4plots <- arrange(df4plots, Final_Assignment)
df4plots$Final_Assignment <- factor(df4plots$Final_Assignment, levels = c(unique(df4plots$Assignment), unique(df4plots$Changed)))


##### Set up population colors #####
matching_colors <- c("gray88", colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(unique(scores$Changed))-1))
names(matching_colors) <- sort(unique(scores$Changed))

pop_colors <- brewer.pal(length(unique(scores$Assignment)), "Dark2")
names(pop_colors) <- unique(scores$Assignment)

colors <- c(matching_colors, pop_colors)

plot_PCs_medoids <- ggplot(df4plots, aes(PC1, PC2, color = Final_Assignment)) +
  geom_point() +
  theme_bw() +
  facet_wrap(vars(Plot)) +
  scale_color_manual(values = colors)


ggsave(plot_PCs_medoids, filename = paste0(args$out, "Ancestry_PCAs.png"), height = 5, width = 12)



##### Subset the Mismatching Ancestry and Sex data #####
sex_mismatch <- sex_check[which(sex_check$STATUS == "PROBLEM"), ]
print(head(sex_mismatch))
if (nrow(sex_mismatch) > 0){
  sex_mismatch$`UPDATE/REMOVE/KEEP` <- NA
} else {
   sex_mismatch$`UPDATE/REMOVE/KEEP` <- character()
}
colnames(sex_mismatch)[1] <- paste0("#", colnames(sex_mismatch)[1])
write_delim(sex_mismatch, paste0(args$out, "/sex_update_remove.tsv"), na = "", delim = "\t")

anc_mismatch <- df4plots[which(df4plots$Plot == "Projected Data Assignments vs Original Assignments" & df4plots$Final_Assignment != "Matched"), ]
anc_temp <- df4plots[, c("FID", "IID", "Assignment")]
colnames(anc_temp) <- c("FID", "IID", "PCA_Assignment")

anc_mismatch <- left_join(anc_mismatch[, c("FID", "IID")], data_anc)
anc_mismatch <- left_join(anc_mismatch, anc_temp)
anc_mismatch <- unique(anc_mismatch)

if (nrow(anc_mismatch) > 0){
  anc_mismatch$`UPDATE/REMOVE/KEEP` <- NA
} else {
   anc_mismatch$`UPDATE/REMOVE/KEEP` <- character()
}
colnames(anc_mismatch)[1] <- paste0("#", colnames(anc_mismatch)[1])

anc_mismatch <- anc_mismatch[, c("#FID", "IID", "Provided_Ancestry", "PCA_Assignment", "UPDATE/REMOVE/KEEP")]

print("writing acestry_update_remove.tsv file")
write_delim(anc_mismatch, paste0(args$out, "/ancestry_update_remove.tsv"), na = "", delim = "\t")

anc_mafs <- data.frame(Ancestry = c("AFR", "AMR", "EAS", "EUR", "SAS"), MAF = NA)
print("writing ancestry_mafs.tsv file")
write_delim(anc_mafs, paste0(args$out, "/ancestry_mafs.tsv"), na = "", delim = "\t")

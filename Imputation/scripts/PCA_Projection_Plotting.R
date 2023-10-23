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
parser$add_argument("-o", "--out", required=TRUE, type="character", help="The output directory where results will be saved.")

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

args$scores <- "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/2023-10-20-ImputationTestDataset/Step1-Imputation/pca_projection/subset_pruned_data_pcs.sscore"
args$onekg_scores <- "/groups/umcg-biogen/tmp01/output/2022-09-01-scMetaBrainConsortium/2023-09-06-scMetaBrain-WorkGroup1QC/2023-10-20-ImputationTestDataset/Step1-Imputation/pca_projection/subset_pruned_1000g_pcs_projected.sscore"

print("Options in effect:")
paste0("  --scores ", args$scores)
paste0("  --onekg_scores ", args$onekg_scores)
paste0("  --out ", args$out)
print("")

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(RColorBrewer)))

## make sure the directory exists ###
dir.create(args$out, recursive=TRUE)

##### Read in Results #####
data_score <- read_delim(as.character(args$scores), delim = "\t", col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "Provided_Ancestry" = "c"))
onekg_score <- read_delim(as.character(args$onekg_scores), delim = "\t", col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "#IID" = "c", "SuperPop" = "c", "Population" = "c"))

##### Set up variables #####

##### Read in PCA Results #####

##### Remove # from colnames #####
colnames(data_score) <- gsub("#", "", colnames(data_score))
colnames(onekg_score) <- gsub("#", "", colnames(onekg_score))

print("data_score")
print(head(as.data.frame(data_score)))
print("onekg_score")
print(head(as.data.frame(onekg_score)))

##### Combine PCA Results Together #####
if (!"FID" %in% colnames(onekg_score)) {
  onekg_score$FID <- NA
}
onekg_score$Provided_Ancestry <- NA
onekg_score <- onekg_score[, c("FID", "IID", colnames(onekg_score)[grep("PC", colnames(onekg_score))], "SuperPop", "Provided_Ancestry")]
colnames(onekg_score) <- gsub("_AVG", "", colnames(onekg_score))

if (!"FID" %in% colnames(data_score)) {
  data_score$FID <- NA
}
data_score$SuperPop <- NA
data_score <- data_score[, c("FID", "IID", colnames(data_score)[grep("PC", colnames(data_score))], "SuperPop", "Provided_Ancestry")]
colnames(data_score) <- gsub("_AVG", "", colnames(data_score))

scores <- rbind(onekg_score, data_score)
print(head(as.data.frame(scores)))

##### Calculate Medoids and Assign Clusters #####
pam_res <- pam(scores[, grep("PC", colnames(scores))], 6)


##### Assign Ancestries to Individuals #####
scores$Cluster <- factor(pam_res$clustering)
print(head(as.data.frame(scores)))

conversion_table <- table(scores$SuperPop, scores$Cluster)
conversion_key <- data.frame(Cluster = colnames(conversion_table), Assignment = rownames(conversion_table)[apply(conversion_table,2, which.max)])
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

##### Subset the Mismatching Ancestry #####
anc_mismatch <- df4plots[which(df4plots$Plot == "Projected Data Assignments vs Original Assignments" & df4plots$Final_Assignment != "Matched"), c("FID", "IID", "Provided_Ancestry", "Assignment")]
colnames(anc_mismatch) <- c("#FID", "IID", "Provided_Ancestry", "PCA_Assignment")
if (nrow(anc_mismatch) > 0){
  anc_mismatch$`UPDATE/REMOVE/KEEP` <- NA
} else {
   anc_mismatch$`UPDATE/REMOVE/KEEP` <- character()
}

print("writing acestry_update_remove.tsv file")
write_delim(anc_mismatch, paste0(args$out, "/ancestry_update_remove.tsv"), na = "", delim = "\t")

anc_mafs <- data.frame(Ancestry = rownames(conversion_table), MAF = NA)
print("writing ancestry_mafs.tsv file")
write_delim(anc_mafs, paste0(args$out, "/ancestry_mafs.tsv"), na = "", delim = "\t")

.libPaths("/usr/local/lib/R/site-library")
library(data.table)
library(ggplot2)
library(cowplot)
library(RColorBrewer)



##### Set directories #####
args <- commandArgs(trailingOnly=TRUE)

dir <- args[1]
outdir <- args[2]




##### Read in Results #####
sex_check <- fread(paste0(dir,"/check_sex/check_sex.sexcheck.tsv"), sep = "\t", colClasses=list(character = c("FID", "IID", "STATUS"), double=c("PEDSEX","SNPSEX","F")))
colnames(sex_check)[1] <- "#FID"
sex_decisions <- fread(paste0(dir,"/pca_sex_checks/check_sex_update_remove.tsv"), sep = "\t", colClasses=list(character = c("#FID", "IID", "STATUS", "UPDATE/REMOVE/KEEP"), double=c("PEDSEX","SNPSEX","F")))

psam <- fread(paste0(dir,"/indiv_missingness/indiv_missingness.psam"), colClasses=list(character = c("#FID", "IID","Provided_Ancestry")))
ancestry_decisions <- fread(paste0(dir,"/pca_sex_checks/ancestry_update_remove.tsv"), colClasses=list(character = c("#FID", "IID", "Provided_Ancestry", "PCA_Assignment", "UPDATE/REMOVE/KEEP")))




##### Set up Dataframes #####
### Sex ###
sex_dt1 <- sex_check[,c('#FID', 'IID', 'PEDSEX')]
sex_dt1$ProvidedSNP <- "Provided"
colnames(sex_dt1)[3] <- "Sex"

sex_dt2 <- sex_check[,c('#FID', 'IID', 'SNPSEX')]
sex_dt2$ProvidedSNP <- "Final"
colnames(sex_dt2)[3] <- "Sex"

sex_dt <- rbind(sex_dt1, sex_dt2)

sex_dt$Sex <- ifelse(sex_dt$Sex == 0, "unknown", ifelse(sex_dt$Sex == 1, "male","female"))
sex_dt$Sex <- factor(sex_dt$Sex, levels = c("female", "male"))
sex_dt$ProvidedSNP <- factor(sex_dt$ProvidedSNP, levels = c("Provided", "Final"))

sex_dt <- sex_decisions[,c('#FID', 'IID','UPDATE/REMOVE/KEEP')][sex_dt, on = c('#FID', 'IID')]

sex_dt <- sex_dt[(is.na(`UPDATE/REMOVE/KEEP`) | `UPDATE/REMOVE/KEEP` == "UPDATE" | `UPDATE/REMOVE/KEEP` == "KEEP") | (`UPDATE/REMOVE/KEEP` == "REMOVE" & ProvidedSNP == "Provided")]



### Ancestry ###
ancestry <- psam[,c('#FID', 'IID', 'Provided_Ancestry')]
ancestry_decisions <- ancestry_decisions[,c('#FID', 'IID', 'PCA_Assignment', 'UPDATE/REMOVE/KEEP')]

ancestry_combined <- ancestry_decisions[ancestry, on = c('#FID', 'IID')]
ancestry_combined$PCA_Assignment <- ifelse(is.na(ancestry_combined$PCA_Assignment), ancestry_combined$Provided_Ancestry, ancestry_combined$PCA_Assignment)

ancestry_combined_long <- melt(ancestry_combined, id.vars = c('#FID', 'IID', 'UPDATE/REMOVE/KEEP'),
                measure.vars = c("Provided_Ancestry", "PCA_Assignment"), variable.name = "ProvidedSNP", value.name = "Ancestry")

ancestry_combined_long <- ancestry_combined_long[(is.na(`UPDATE/REMOVE/KEEP`) | `UPDATE/REMOVE/KEEP` == "UPDATE" | `UPDATE/REMOVE/KEEP` == "KEEP") | (`UPDATE/REMOVE/KEEP` == "REMOVE" & ProvidedSNP == "Provided")]



##### Make Figures #####
### Sex ###
p_sex1 <- ggplot(sex_dt, aes(x = Sex, fill = ProvidedSNP)) +
    geom_bar( position=position_dodge()) +
    theme_classic() +
    theme(axis.title.x=element_blank()) +
    ylab("Number Individuals") +
    scale_fill_manual(values = c("#dea9cd", "#7dc5c6"))

if(nrow(sex_decisions) > 0){
    p_sex2 <- ggplot(sex_decisions, aes(x = `UPDATE/REMOVE/KEEP`)) +
        geom_bar(fill = "grey60") +
        theme_classic() +
        theme(axis.title.x=element_blank()) +
        ylab("Number Individuals") +
        ggtitle(paste0("Nonmatching Sex\n(",nrow(sex_decisions), " of ", nrow(sex_check), " individuals)")) +
        theme(plot.title = element_text(hjust = 0.5)) 

    ggsave(plot_grid(p_sex1, p_sex2, rel_heights = c(2.5, 2), ncol = 1, nrow = 2, align = "v", axis = "lr"), filename = paste0(outdir, "/sex_summary.png"), width = 4)


} else {
    ggsave(p_sex1, filename = paste0(outdir, "/sex_summary.png"))
}


### Ancestry ###
p_ancestry1 <- ggplot(ancestry_combined_long, aes(x = Ancestry, fill = ProvidedSNP)) +
    geom_bar(position=position_dodge()) +
    theme_classic() +
    theme(axis.title.x=element_blank()) +
    ylab("Number Individuals")+
    scale_fill_brewer(palette = "Set2")


if(nrow(ancestry_decisions) > 0){
    p_ancestry2 <- ggplot(ancestry_combined_long[!is.na(`UPDATE/REMOVE/KEEP`) & ProvidedSNP == "PCA_Assignment"], aes(x = `UPDATE/REMOVE/KEEP`)) +
        geom_bar(fill = "grey60") +
        theme_classic() +
        theme(axis.title.x=element_blank()) +
        ylab("Number Individuals") +
        ggtitle(paste0("Nonmatching Ancestry\n(",nrow(ancestry_decisions), " of ", nrow(ancestry), " individuals)")) +
        theme(plot.title = element_text(hjust = 0.5)) 


    ggsave(plot_grid(p_ancestry1, p_ancestry2, rel_heights = c(2.5, 2), ncol = 1, nrow = 2, align = "v", axis = "lr"), filename = paste0(outdir, "/ancestry_summary.png"), width = 2 + length(unique(ancestry_combined_long$Ancestry)))


} else {
    ggsave(p_ancestry1, filename = paste0(outdir, "/sex_summary.png"))
}

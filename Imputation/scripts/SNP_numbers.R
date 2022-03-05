.libPaths("/usr/local/lib/R/site-library")
library(data.table)
library(ggplot2)


##### Set directories #####
args <- commandArgs(trailingOnly=TRUE)

dir <- args[1]
outdir <- args[2]



##### Make a list of the ancestry files #####
vcf_list <- paste0(dir,"/vcf_merged_by_ancestries/",list.files(paste0(dir,"/vcf_merged_by_ancestries/"), pattern = "_imputed_hg38.vcf.gz$"))
names(vcf_list) <- paste0(gsub(paste0(dir,"/vcf_merged_by_ancestries/"), "", gsub("_imputed_hg38.vcf.gz$","", vcf_list)),"\nImputed\nSNPs")

vcf_list[["All\nImputed\nSNPs"]] <- paste0(dir,"/vcf_all_merged/imputed_hg38_info_filled.vcf.gz")
vcf_list[["All\nImputed\nSNPs\nFiltered"]] <- paste0(dir,"/vcf_all_merged/imputed_hg38_R2_0.3_MAF0.05.vcf.gz")
vcf_list[["All\nImputed\nSNPs\nFiltered\nGenes"]] <- paste0(dir,"/vcf_4_demultiplex/imputed_hg38_R2_0.3_MAF0.05_exons_sorted.vcf")



vcf_length_list <- lapply(vcf_list, function(x){
    if(grepl(".gz$",x)){
        read.table(pipe(paste0("gunzip -c ", x, " | grep -v '#' | wc -l")))[[1]]
    } else {
        read.table(pipe(paste0("grep -v '#' ", x, " | wc -l")))[[1]]
    }
})


vcf_length <- unlist(vcf_length_list)


dt <- data.table(Group = names(vcf_length_list), N_SNPs = vcf_length)
dt$Group <- factor(dt$Group, levels = names(vcf_length_list))


plot <- ggplot(dt, aes(Group, log10(N_SNPs))) +
    geom_bar(stat = "identity") +
    theme_classic() +
    ylab("log(Number SNPs)") +
    theme(axis.title.x=element_blank()) +
    geom_text(aes(label = prettyNum(N_SNPs,big.mark=",",scientific=FALSE), y = 0.1, angle = 90, hjust = 0), color = "white")

ggsave(plot, filename = paste0(outdir, "/Number_SNPs.png"), width = 0.5 + 0.6 * length(dt$Group), height = 4)



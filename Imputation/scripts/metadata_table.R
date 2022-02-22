library(tidyverse)


##### Read in argument variables #####
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
psam <- as.character(arguments[1,])
outfile <- arguments[2,]



##### Read in the data #####
psam_df <- read_delim(psam, delim = "\t", col_types = cols("#FID" = "c", "IID" = "c", "Ancestry" = "c"))

meta <- data.frame("donor_id" = paste0(psam_df$`#FID`,"_",psam_df$IID), "sex" = gsub(1,"M",gsub(2,"F",psam_df$SEX)), psam_df[,6:ncol(psam_df)])
colnames(meta) <- gsub("^Ancestry$", "ethnicity_super_population_code", colnames(meta))

meta$imputation_server <- "SangerImputationServer"

write_delim(meta, outfile, delim = "\t")
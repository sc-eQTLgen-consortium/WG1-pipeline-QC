# Script information ------------------------------------------------------

# title: Remove individuals by heterozygosity rate
# author: Jose Alquicira Hernandez
# date: 2019/03/12
# description: Removing individuals who deviate ±3 SD fromthe samples' heterozygosity 
# rate mean


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("magrittr")

# Manage arguments --------------------------------------------------------

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]


# Read data ---------------------------------------------------------------

het <- read.table(input, head = TRUE)


# Remove individuals ------------------------------------------------------

het$HET_RATE <- (het$"N.NM." - het$"O.HOM.") / het$"N.NM."
het_fail <- subset(het, (het$HET_RATE < mean(het$HET_RATE) - 3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE) +3 * sd(het$HET_RATE)));
het_fail$HET_DST <- (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)


# Write results -----------------------------------------------------------

write.table(het_fail, output, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

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
output_fail <- args[2]
output_pass <- args[3]
output_pass_inds <- args[4]


# Read data ---------------------------------------------------------------

het <- read_delim(input, delim = "\t")


if (nrow(het) >= 3){
	# Remove individuals ------------------------------------------------------

	het$HET_RATE <- (het$"N_SITES" - het$"O(HOM)") / het$"N_SITES"
	het_fail <- subset(het, (het$HET_RATE < mean(het$HET_RATE) - 3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE) +3 * sd(het$HET_RATE)));
	het_pass <- subset(het, (het$HET_RATE > mean(het$HET_RATE) - 3*sd(het$HET_RATE)) & (het$HET_RATE < mean(het$HET_RATE) +3 * sd(het$HET_RATE)));
	het_fail$HET_DST <- (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)
	het_pass$HET_DST <- (het_pass$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)


	# Write results -----------------------------------------------------------

	write.table(het_fail, output_fail, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	write.table(het_pass, output_pass, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	write.table(het_pass$INDV, output_pass_inds, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

} else {
	write.table(het$INDV, output_pass_inds, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	system(paste0("touch ", output_fail))
	system(paste0("touch ", output_pass))
}
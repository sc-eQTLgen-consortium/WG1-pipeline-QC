# Script information ------------------------------------------------------

# title: Remove individuals by heterozygosity rate
# author: Jose Alquicira Hernandez
# date: 2019/03/12
# description: Removing individuals who deviate ±3 SD fromthe samples' heterozygosity 
# rate mean

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# Manage arguments --------------------------------------------------------

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, type="character", help="")
parser$add_argument("-t", "--threshold", required=TRUE, default=3, type="integer", help="")
parser$add_argument("-o", "--out", required=FALSE, type="character", help="")

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

print("Options in effect:")
paste0("  --input ", args$input)
paste0("  --out ", args$out)


# Import libraries --------------------------------------------------------

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(magrittr)))

# Read data ---------------------------------------------------------------

het <- read_delim(args$input, delim = "\t")


if (nrow(het) >= 3){
	# Remove individuals ------------------------------------------------------

	het$HET_RATE <- (het$"N_SITES" - het$"O(HOM)") / het$"N_SITES"
	het_fail <- subset(het, (het$HET_RATE < mean(het$HET_RATE) - args.threshold * sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE) + args.threshold * sd(het$HET_RATE)));
	het_pass <- subset(het, (het$HET_RATE > mean(het$HET_RATE) - args.threshold * sd(het$HET_RATE)) & (het$HET_RATE < mean(het$HET_RATE) + args.threshold * sd(het$HET_RATE)));
	het_fail$HET_DST <- (het_fail$HET_RATE - mean(het$HET_RATE)) / sd(het$HET_RATE)
	het_pass$HET_DST <- (het_pass$HET_RATE - mean(het$HET_RATE)) / sd(het$HET_RATE)


	# Write results -----------------------------------------------------------

	write.table(het_fail, paste0(args$out, "_failed.inds"), row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	write.table(het_pass, paste0(args$out, "_passed.inds"), row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	write.table(het_pass$INDV, paste0(args$out, "_passed.txt"), row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

} else {
	write.table(het$INDV, paste0(args$out, "_passed.txt"), row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
	system(paste0("touch ", paste0(args$out, "_failed.inds")))
	system(paste0("touch ", paste0(args$out, "_passed.inds")))
}
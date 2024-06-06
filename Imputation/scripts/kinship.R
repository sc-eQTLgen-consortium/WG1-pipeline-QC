#!/usr/bin/env Rscript
# Author:
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-k", "--king", required=TRUE, type="character", help="")
parser$add_argument("-ki", "--king_id", required=TRUE, type="character", help="")
parser$add_argument("-o", "--out", required=TRUE, type="character", help="")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

print("Options in effect:")
for (name in names(args)) {
	print(paste0("  --", name, " ", args[[name]]))
}
print("")

# check date is provided
if (any(unlist(lapply(c(args$king, args$king_id, args$out), is.na)))) {
  stop("required parameters (king, king_id & out) must be provided.")
}

kin <- read.delim(args$king,header=F)
kinIds <- read.delim(args$king_id)

kin <- kin * 2

rownames(kin) <- colnames(kin) <- kinIds[,1]

write.table(kin, args$out, sep="\t", quote=F, col.names=NA)

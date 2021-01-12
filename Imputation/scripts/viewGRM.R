library(tidyverse)

args <- commandArgs(TRUE)
in_file <- args[1]
out <- args[2]



##### Functions #####
omxReadGRMBin_updated <- function (prefix, AllN = FALSE, size = 4, returnList = FALSE) {sum_i = function(i) {
        return(sum(1:i))
    }
    BinFileName = paste(prefix, ".grm.bin", sep = "")
    NFileName = paste(prefix, ".grm.N.bin", sep = "")
    IDFileName = paste(prefix, ".grm.id", sep = "")
    id = read.table(IDFileName)
    n = dim(id)[1]
    BinFile = file(BinFileName, "rb")
    grm = readBin(BinFile, n = n * (n + 1)/2, what = numeric(0), 
        size = size)
    NFile = file(NFileName, "rb")
    if (AllN == T) {
        N = readBin(NFile, n = n * (n + 1)/2, what = numeric(0), 
            size = size)
    }
    else {
        N = readBin(NFile, n = 1, what = numeric(0), size = size)
    }
    i = sapply(1:n, sum_i)
    if (returnList) {
        return(list(diag = grm[i], off = grm[-i], id = id, N = N))
    }
    else {
        GRM <- matrix(0, nrow = n, ncol = n, dimnames = list(paste0(id$V1, "_", id$V2), paste0(id$V1, "_", id$V2)))
        GRM[!lower.tri(GRM, diag = T)] <- grm[-i]
        GRM <- GRM + t(GRM)
        diag(GRM) <- grm[i]
        return(GRM)
    }
}


##### Read in files #####
id <- read_delim(paste0(in_file, ".grm.id"), delim = "\t", col_names = c("FID","IID"))

matrix <- omxReadGRMBin_updated(in_file)
pairs <- data.frame(t(combn(paste0(id$FID, "_",id$IID), 2, simplify = TRUE)))
colnames(pairs) <- c("Individual1", "Individual2")
pairs$genetic_relatedness <- t(matrix)[lower.tri(t(matrix))]

pairs <- pairs[rev(order(pairs$genetic_relatedness)),]

write_delim(pairs, paste0(out,"paired_relatedness.tsv"), delim = "\t")
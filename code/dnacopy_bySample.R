
options(stringsAsFactors = FALSE)

library(DNAcopy)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

file.name     <- args[1]
output.folder <- args[2]
sample.name   <- gsub("^.+?_([^_]+).txt", "\\1", file.name)
#output.folder <- "results/dnacopy"
#plots.folder  <- "plots/dnacopy"

#X <- read.delim(file.name, check.names = FALSE)
X <- read.delim(file.name)

#array.info <- read.delim("data/Marazitafacialvariation_11779814_A.pfb")
#array.info = array.info[, -4]

#X <- left_join(X, array.info, by = c('Name' = 'SNP'))
#X <- left_join(X, array.info, by = c('SNP.Name' = 'Name'))

# remove probes labeled chr 0, these might be for quality control
X <- X[X$Chr != 0, ]

#lrr.column <- grep("Log.R.Ratio", colnames(X))
#sample.name <- gsub(".Log.R.Ratio", "", colnames(X)[lrr.column])

lrr.column <- grep("corrected.vals", colnames(X))
sample.name <- X[1, "Sample.Name"]

CNA.object <- CNA( 
    cbind(X[,lrr.column]), 
    X$Chr,
    X$MapInfo, 
    data.type = "logratio",
    sampleid = sample.name
)

## Clean Data
smoothed.CNA.object <- smooth.CNA(CNA.object)
segments <- segment(smoothed.CNA.object, min.width = 5, alpha = 0.1, verbose = 1)


write.table(segments$output, file = paste0(output.folder, "/segment_alpha0_1_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE) 



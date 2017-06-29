
options(stringsAsFactors = FALSE)

library(DNAcopy)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
file.name     <- args[1]
output.folder <- args[2]
sample.name   <- gsub("^.+?_([^_]+).txt", "\\1", file.name)

X <- read.delim(file.name)

# Assume we are using the GC corrected values generated from ArrayTV
lrr.column <- "corrected.vals"
if (! lrr.column %in% colnames(X))
{
	warning("The corrected.vals column was not found. For best results, perform GC correction before running DNAcopy.")
	lrr.column <- grep(".Log.R.Ratio$", colnames(X), value = TRUE)
}
stopifnot(all(c("Chr", "MapInfo", "Sample.Name", lrr.column) %in% colnames(X)))

sample.name <- X[1, "Sample.Name"]

# remove probes labeled chr 0, these might be for quality control
X <- filter(X, Chr != 0)

CNA.object <- CNA( 
    cbind(X[[lrr.column]]), 
    X[["Chr"]],
    X[["MapInfo"]], 
    data.type = "logratio",
    sampleid = sample.name
)

## Clean Data
smoothed.CNA.object <- smooth.CNA(CNA.object)
segments <- segment(smoothed.CNA.object, min.width = 5, alpha = 0.1, verbose = 1)

write.table(segments$output, file = paste0(output.folder, "/segment_alpha0_1_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE) 


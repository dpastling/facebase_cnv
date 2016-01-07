
options(stringsAsFactors = FALSE)

library(DNAcopy)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

file.name     <- args[1]
output.folder <- args[2]
sample.name   <- gsub("^.+?_([^_]+).txt", "\\1", file.name)
#output.folder <- "results/dnacopy"
plots.folder  <- "plots/dnacopy"

#X <- read.delim(file.name, check.names = FALSE)
X <- read.delim(file.name)

array.info <- read.delim("data/Human_Omni25exome.pfb")
array.info = array.info[, -4]

X <- left_join(X, array.info, by = c('Name' = 'SNP'))

# remove probes labeled chr 0, these might be for quality control
X <- X[X$Chr != 0, ]

lrr.column <- grep("Log.R.Ratio", colnames(X))
sample.name <- gsub(".Log.R.Ratio", "", colnames(X)[lrr.column])

CNA.object <- CNA( 
    cbind(X[,lrr.column]), 
    X$Chr,
    X$Position, 
    data.type = "logratio",
    sampleid = sample.name
)

## Clean Data
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = 1)



## Visualize Data

custom_png <- function(file)
{
	png(file, width = 1800, height = 1500, res = 300, bg = "white")
}

### Plot all chromosomes on one scatter plot
custom_png(file = paste0(plots.folder, "/genome_wide_", sample.name, ".png"))
plot(segment.smoothed.CNA.object, plot.type = "w")
dev.off()

### Individual faceted plots for each chromosome 
custom_png(file = paste0(plots.folder, "/by_chr_", sample.name, ".png"))
plot(segment.smoothed.CNA.object, plot.type = "s")
dev.off()

### Don't understand the point of this plot, maybe means for each chromosome?
custom_png(file = paste0(plots.folder, "/segment_means_", sample.name, ".png"))
plot(segment.smoothed.CNA.object, plot.type = "p")
dev.off()


# Remove segments that are not at least 3 SDs apart
sdundo.CNA.object <- segment(
	smoothed.CNA.object,
	undo.splits = "sdundo",
	undo.SD = 3,
	verbose = 1
)

segment.smoothed.p <- segments.p(segment.smoothed.CNA.object)
sdundo.p <- segments.p(sdundo.CNA.object)

write.table(segment.smoothed.CNA.object$output, file = paste0(output.folder, "/segment_smoothed_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE) 
write.table(sdundo.CNA.object$output, file = paste0(output.folder, "/sdundo_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(segment.smoothed.p, file = paste0(output.folder, "/segment_pvalue_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sdundo.p, file = paste0(output.folder, "/sdundo_pvalue_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

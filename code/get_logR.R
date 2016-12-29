
options(stringsAsFactors = FALSE)
#library(plyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 1)


sample.name      <- args[1]
sample.name.tidy <- gsub("@", ".", sample.name)
output.folder    <- "results/logR"

# for penncnv and VICE the first three columns are the chrom, pos1, pos2, like in a bed file
# Penncnv doesn't have headers
# the first column of dnacopy is the sample name.

penncnv <- read.delim(paste0("results/penncnv/sample_", sample.name, ".tabcnv"), header = FALSE)
colnames(penncnv) <- c("chr", "pos1", "pos2", "state", "file", "start_probe", "end_probe", "score", "n")

vice <- read.delim(paste0("results/vanilla_ice/VICE_", sample.name, ".txt"))
vice <- mutate(vice, seqnames = gsub("chr", "", seqnames))

dnacopy <- read.delim(paste0("results/dnacopy_alpha0_1/segment_alpha0_1_", sample.name, ".txt"))
# Note this filtering needs to be done elsewhere
#dnacopy <- filter(dnacopy, num.mark >= 5)
corrected.vals <- read.delim(paste0("data/corrected_vals/corrected_vals_", sample.name, ".txt"))

get_logR <- function(chr, pos1, pos2)
{
	corrected.vals %>% 
	filter(Chr == chr & MapInfo >= pos1 & MapInfo <= pos2) %>% 
	summarise(
		mean.logR = mean(Log.R.Ratio, na.rm = TRUE), 
		mean.norm.logR = mean(corrected.vals, na.rm = TRUE)
#		probes = n()
	) %>%
	as.data.frame()
}

# There's probably a more elegant way to do this, but just am brute forcing it

foo <- NULL
for (i in 1:nrow(vice))
{
	foo <- rbind(foo, get_logR(vice[i, "seqnames"], vice[i, "start"], vice[i, "end"]))
}
vice <- cbind(vice, foo)

foo <- NULL
for (i in 1:nrow(penncnv))
{
	foo <- rbind(foo, get_logR(penncnv[i, "chr"], penncnv[i, "pos1"], penncnv[i, "pos2"]))
}
penncnv <- cbind(penncnv, foo)

foo <- NULL
for (i in 1:nrow(dnacopy))
{
	foo <- rbind(foo, get_logR(dnacopy[i, "chrom"], dnacopy[i, "loc.start"], dnacopy[i, "loc.end"]))
}
dnacopy <- cbind(dnacopy, foo)


write.table(penncnv, file = paste0(output.folder, "/penncnv/", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dnacopy, file = paste0(output.folder, "/dnacopy/", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(vice, file = paste0(output.folder, "/vice/", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)




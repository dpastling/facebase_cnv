
options(stringsAsFactors = FALSE)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
sample.name <- args[1]
bed.file    <- args[2]
result.file <- args[3]
sample.name.tidy <- gsub("@", ".", sample.name)

get_logR <- function(chr, pos1, pos2)
{
	corrected.vals %>% 
	filter(Chr == chr & MapInfo >= pos1 & MapInfo <= pos2) %>% 
	summarise(
		mean.logR = mean(Log.R.Ratio, na.rm = TRUE), 
		mean.norm.logR = mean(corrected.vals, na.rm = TRUE), 
		probes = n()
	) %>%
	as.data.frame()
}

bed <- read.delim(bed.file, header = FALSE)
corrected.vals <- read.delim(paste0("data/corrected_vals/corrected_vals_", sample.name, ".txt"))
log.R <- NULL
for (i in 1:nrow(bed))
{
        log.R <- rbind(log.R, get_logR(bed[i, 1], bed[i, 2], bed[i, 3]))
}
bed <- cbind(bed, log.R)
write.table(bed, file = result.file, sep = "\t", quote = FALSE, row.names = FALSE)


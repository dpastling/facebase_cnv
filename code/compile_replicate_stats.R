
options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

results.folder <- "results/vanilla_ice_qc"

files <- list.files(results.folder, pattern = "((del)|(dup)).bed", full.names = TRUE)

df <- lapply(files, function(x) {
	X <- read.delim(x, header = FALSE)
	X <- X[,c(1,2,3,8,11,12)]
})
files <- gsub(results.folder, "", files)
files <- gsub("/", "", files)
files <- gsub(".bed", "", files)
names(df) <- files
df <- bind_rows(df, .id = "sample")
colnames(df) <- c("sample", "chr", "pos1", "pos2", "stat", "logR", "probes")
df <- mutate(df, length = pos2 - pos1)

df.stats <- df %>% 
	group_by(sample) %>%
	summarise(
		cnvs = n(), 
		mean.logR = mean(logR), 
		max.logR = max(logR), 
		min.logR = min(logR), 
		mean.len = mean(length), 
		mean.probes = mean(probes),
		mean.stat = mean(stat),
		min.stat = min(stat),
		max.stat = max(stat)
	)


write.table(df, file = "combined_replicates.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df.stats <- mutate(df.stats, hapmap.id = gsub("^(\\d+)[_\\.].+?$", "\\1", sample))
df.stats <- mutate(df.stats, alg = gsub("^.+?_([^_]+?)_((del)|(dup))", "\\1", sample))
df.stats <- mutate(df.stats, type = gsub("^.+?_([^_]+?)_((del)|(dup))", "\\2", sample))
df.stats <- mutate(df.stats, combo = ifelse(grepl("^\\d+_", sample), "two", ""))
df.stats <- mutate(df.stats, combo = ifelse(grepl("^\\d+\\.", sample), "single", combo))
df.stats <- mutate(df.stats, combo = ifelse(grepl("_all_", sample), "all", combo))
df.stats <- arrange(df.stats, alg, hapmap.id, type, combo)

write.table(df.stats, file = "replicate_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)






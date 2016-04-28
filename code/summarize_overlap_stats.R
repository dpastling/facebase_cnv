
options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

df <- read.delim("stats_20160416.txt", header = FALSE)
colnames(df) <- c("file", "chr", "pos1", "pos2", "n", "stat")

df <- mutate(df, file = gsub("results/vanilla_ice_qc/", "", file))
df <- mutate(df, file = gsub("_overlap.txt", "", file))
df <- mutate(df, file = gsub(".bed", "_bed", file))
df <- filter(df, pos2 - pos1 > 0)

df <- df %>% 
    group_by(file) %>% 
	summarise(
		number.of.regions = length(pos1),
		mean.log.length = mean(log(pos2 - pos1, 10)),
		mean.length = mean(pos2 - pos1),
		mean.probes = mean(n),
		mean.score  = mean(stat, na.rm = TRUE)
		)

df <- separate(df, file, c("patient", "id", "alg.1", "alg.2"), sep = "[_\\.]")

write.table(df, file = "full_comparison.txt", sep = "\t", quote = FALSE, row.names = FALSE)




options(stringsAsFactors = FALSE)
library(dplyr)

annotation.file <- "data/Marker_Info_Files/HumanOmni25Exome-8v1_A.csv"

features <- read.csv(annotation.file, header = TRUE, skip = 7)
features <- features %>%
                rename(Position = MapInfo) %>%
                filter(! Chr %in% c("0", "MT", "X", "XY", "Y")) %>%
                filter(! is.na(Position)) %>%
                mutate(Intensity.Only = as.integer(! grepl("\\[", SNP))) %>%
                filter(Chr == "16") %>% 
                arrange(as.numeric(Chr), Position) %>%
                select(Chr, Position, Name, Intensity.Only)

#features <- features[1:10000, ]

write.table(features, file = "SNP_info.csv", sep = ",", row.names = FALSE, quote = FALSE)

extdir <- "data/Spritz_release_genotype_files"
data.files <- list.files(extdir, pattern = "FinalReport", full.names = TRUE)
samples.of.interest <- c(
    "06985",
    "18502",
    "18505",
    "18508",
    "18858",
    "18861"
)
sample.names <- gsub("^.+?_([^_]+)@.+?$", "\\1", data.files)
data.files <- data.files[sample.names %in% samples.of.interest]
#data.files <- data.files[1:5]

for (i in 1:length(data.files))
{
	X <- read.csv(data.files[i], header = TRUE, check.names = FALSE, skip = 10)
	X <- X[X[, "SNP Name"] %in% features[, "Name"], ]
	write.table(X, file = paste0("sample_", i, ".txt"), sep = ",", quote = FALSE, row.names = FALSE)
}



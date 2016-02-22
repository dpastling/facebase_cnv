#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(plyr)
library(dplyr)

file <- "data/Spritz_release_genotype_files/Spritz_Omni2-5PlusExome_release_FinalReport_5566@1070895973.csv"
annotation.data <- "data//Marker_Info_Files/HumanOmni25Exome-8v1_A.csv"
annotation <- read.csv(annotation.data, header = TRUE, skip = 7)
annotation <- select(annotation, Name, Chr, MapInfo)

#X <- read.csv(file, header = TRUE, skip = 10)
#X <- left_join(X, annotation, by = c("SNP.Name" = "Name"))

calculate_MAD <- function(x) abs(x - median(x, na.rm = TRUE))

calculate_acf_lag10 <- function(x) {
    x.acf <- acf(x, lag.max = 10, na.action = na.exclude, plot = FALSE)
    i <- which(x.acf$lag == 10)
    return(x.acf$acf[i])
}

files <- list.files(path = "data/Spritz_release_genotype_files", pattern = ".csv$", full.names = TRUE)

samples.of.interest <- c(
    "06985",
    "18502",
    "18505",
    "18508",
    "18858",
    "18861"
)
sample.names <- gsub("^.+?_([^_]+)@.+?$", "\\1", files)
files <- files[sample.names %in% samples.of.interest]

mad.scores <- adply(files, 1, function(x) {
        X <- read.csv(x, header = TRUE, skip = 10)
        X <- left_join(X, annotation, by = c("SNP.Name" = "Name"))
        X <- X %>% 
            tbl_df() %>%
            filter(! Chr %in% c("0", "MT", "X", "XY", "Y")) %>% 
            arrange(as.numeric(Chr), MapInfo) %>%
            mutate(mad = calculate_MAD(Log.R.Ratio))
        sample.name <- gsub("^.+?/([^/]+).csv$", "\\1", x)
        sample.name <- gsub("^.+?FinalReport_(.+?)$", "\\1", sample.name)
        png(file = paste0("plots/mad_plot_", sample.name, ".png"), width = 2000, height = 500)
        plot(X$mad, ylab = "MAD Score")
        dev.off()
        result <- X %>% 
            summarise(mad = median(mad, na.rm = TRUE), acf = calculate_acf_lag10(Log.R.Ratio)) %>%
            mutate(sample = sample.name)
        return(as.data.frame(result))
        })
mad.scores <- mad.scores[, -1]

write.table(mad.scores, file = "qc_scores_controls.txt", sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)




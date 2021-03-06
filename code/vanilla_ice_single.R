
options(stringsAsFactors = FALSE)

library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ArrayTV)
library(VanillaICE)
library(doSNOW)
#library(doParallel)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 3)

# for testing
# file.name <- "data/Spritz_release_genotype_files/Spritz_Omni2-5PlusExome_release_FinalReport_1110@0123841633.csv"
# output.folder <- "results/vanilla_ice"

file.name     <- args[1]
output.folder <- args[2]
features.data <- args[3]
sample.name   <- gsub("^.+?/([^/]+)$", "\\1", file.name)
sample.name   <- gsub("FinalReport_", "", sample.name)
sample.name   <- gsub(".csv$", "", sample.name)
sample.name   <- gsub(".txt$", "", sample.name)

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)

features <- read.csv(features.data, header = TRUE, skip = 7)

# VanillaICE requires a field that identifies a probe as a SNP. All SNPs seem to
# be labeled like [A/T], while control probes are called something like "DNP (Bgnd)"
# we will use the presence of an opening bracket as evidence for a SNP
# Anyway, it turns out that all mapped probes have an IntensityOnly value of 0
features <- dplyr::mutate(features, IntensityOnly = as.integer(! grepl("\\[", SNP)))
features <- dplyr::select(features, Name, Chr, MapInfo, IntensityOnly)

X <- read.csv(file.name, header = TRUE, skip = 10)
X <- inner_join(X, features, by = c("SNP.Name" = "Name"))
X <- X %>% 
    dplyr::filter(! Chr %in% c("0", "MT", "XY")) %>% 
    dplyr::filter(! is.na(Log.R.Ratio)) %>%
    dplyr::arrange(as.numeric(Chr), MapInfo)
X <- as.data.frame(X)

max.window <- c(100, 10e3, 1e6)
increms <- c(20, 2000, 200e3)

# For some reason, ArrayTV requires the data to be in a matrix. 
# the object needs to be a one rolumn matrix for it to work.# alternativly you can pass along objects of type ‘BeadStudioSet’, ‘BafLrrSet’, or ‘BafLrrSetList’
# but is much easier just to select the single column
temp   <- dplyr::select(X, Log.R.Ratio, MapInfo)
temp   <- as.matrix(temp)
tvList <- gcCorrect(object = temp[, "Log.R.Ratio", drop = FALSE],
	chr = paste0("chr", X$Chr),
	starts = X$MapInfo,
	increms = increms,
	maxwins = max.window,
	build = 'hg19'
	)

X <- dplyr::mutate(X, corrected.vals = tvList$correctedVals)

# save corrected values
write.table(X, file = paste0(output.folder, "/corrected_vals_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# plot before and after Log.R.Ratios


######################################
# VanillaICE
######################################
    
# think about the Intensity Only flag
fgr            <- GRanges(
                      paste0("chr", X$Chr), 
                      IRanges(X$MapInfo, width = 1), 
                      isSnp = X$IntensityOnly == 0
                      )
fgr            <- SnpGRanges(fgr)
names(fgr)     <- X[["SNP.Name"]]
sl             <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlevels(fgr) <- sl[sl %in% seqlevels(fgr)]
seqinfo(fgr)   <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(fgr), ]
fgr            <- sort(fgr)


temp <- X[, c("corrected.vals", "B.Allele.Freq")]
temp <- as.matrix(temp)
snp_expt <- SnpArrayExperiment(
                cn = temp[, "corrected.vals", drop = FALSE], 
                baf = temp[, "B.Allele.Freq", drop = FALSE], 
                rowRanges = fgr
            )
param    <- EmissionParam()
fit     <- hmm2(snp_expt, param)

filter_param <- FilterParam(numberFeatures = 5)
result <- cnvSegs(fit, filter_param)
result <- as.data.frame(result)
write.table(result, file = paste0(output.folder, "/VICE_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

stopCluster(cl)


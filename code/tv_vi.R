#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(plyr)
library(dplyr)
library(ArrayTV)
library(VanillaICE)
library(BSgenome.Hsapiens.UCSC.hg19)
library(doSNOW)
#library(doParallel)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

file.name     <- args[1]
output.folder <- args[2]
sample.name   <- gsub("^.+?_([^_]+).csv", "\\1", file.name)

cl <- makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

#cl <- makeCluster(8)
#registerDoParallel(cl)

annotation.data <- "data//Marker_Info_Files/HumanOmni25Exome-8v1_A.csv"
annotation <- read.csv(annotation.data, header = TRUE, skip = 7)

# VanillaICE requires a field that identifies a probe as a SNP. All SNPs seem to
# be labeled like [A/T], while control probes are called something like "DNP (Bgnd)"
# we will use the presence of an opening bracket as evidence for a SNP
# Anyway, it turns out that all mapped probes have an IntensityOnly value of 0
annotation <- mutate(annotation, IntensityOnly = as.integer(! grepl("\\[", SNP)))
annotation <- select(annotation, Name, Chr, MapInfo, IntensityOnly)

X <- read.csv(file.name, header = TRUE, skip = 10)
X <- inner_join(X, annotation, by = c("SNP.Name" = "Name"))
X <- X %>% 
    filter(! Chr %in% c("0", "MT", "X", "XY", "Y")) %>% 
    filter(! is.na(Log.R.Ratio)) %>%
    arrange(as.numeric(Chr), MapInfo)
#X <- as.data.frame(X)

max.window <- c(100, 10e3, 1e6)
increms <- c(20, 2000, 200e3)

# For some reason, ArrayTV requires the data to be in a matrix. 
# the object needs to be a one rolumn matrix for it to work.# alternativly you can pass along objects of type ‘BeadStudioSet’, ‘BafLrrSet’, or ‘BafLrrSetList’
# but is much easier just to select the single column
temp   <- select(X, Log.R.Ratio, MapInfo)
temp   <- as.matrix(temp)
tvList <- gcCorrect(object = temp[, "Log.R.Ratio", drop = FALSE],
	chr = paste0("chr", X$Chr),
	starts = X$MapInfo,
	increms = increms,
	maxwins = max.window,
	build = 'hg19'
	)

X <- mutate(X, corrected.vals = tvList$correctedVals)

# save corrected values
#write.table(X, file = paste0(output.folder, "/corrected_vals_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# plot before and after Log.R.Ratios


######################################
#
# VanillaICE
#
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

## Option 1

#temp <- X[, c("corrected.vals", "B.Allele.Freq")]
temp <- X[, c("Log.R.Ratio", "B.Allele.Freq")]
temp <- as.matrix(temp)

# Test the randomness
#set.seed(12345)
set.seed(54321)

#snp_expt <- SnpArrayExperiment(cn = temp[, "corrected.vals", drop = F], baf = temp[, "B.Allele.Freq", drop = F], rowRanges = fgr)
snp_expt <- SnpArrayExperiment(
                cn = temp[, "Log.R.Ratio", drop = FALSE], 
                baf = temp[, "B.Allele.Freq", drop = FALSE], 
                rowRanges = fgr
            )
param    <- EmissionParam()
fit1     <- hmm2(snp_expt, param)

#save(fit1, file = "hmm_fit1.RData")

stopCluster(cl)

filter_param <- FilterParam(numberFeatures = 5, probability = 0.95)
result <- cnvSegs(fit1, filter_param)
result <- as.data.frame(result)
#write.table(result, file = paste0(output.folder, "/VICE_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result, file = paste0(output.folder, "/VICE_noTV_test3_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)



## Option 2
if (FALSE)
{
parsed_dir <- "vanilla_ice"
sample_id  <- "5566@1070895973"

r <- as.integer(round(1000 * X[, "Log.R.Ratio"], 0))
b <- as.integer(round(1000 * X[, "B.Allele.Freq"], 0))
gt <- paste0(X[, "Allele1...AB"], X[, "Allele2...AB"])
if(any(! gt %in% c("AA", "AB", "BB"))){
    gt.missing <- which(! gt %in% c("AA", "AB", "BB"))
    gt[gt.missing] <- NA
}
gt <- as.numeric(factor(gt), levels = c("AA", "AB", "BB"))
gt <- as.integer(round(1000 * gt, 0))

lrr_file <- file.path(parsed_dir, paste0(sample_id, “_lrr.rds”)
baf_file <- file.path(parsed_dir, paste0(sample_id, “_baf.rds”)
gt_file <- file.path(parsed_dir, paste0(sample_id, “_baf.rds”)
saveRDS(r, file = lrr_file)
saveRDS(b, file = baf_file)
saveRDS(gt, file = gt_file)

# Now propogate the above information to a views object to facilitate easy access to this data:

views <- ArrayViews(rowRanges = fgr, sourcePaths = files, parsedPath = parsedDir)
lrrFile(views) <- file.path(parsed_dir, basename(fileName(views, "lrr")))
views@bafFiles <- file.path(parsed_dir, basename(fileName(views, "baf")))
views@gtFiles <- file.path(parsed_dir, basename(fileName(views, "gt")))
colnames(views) <- gsub(".csv", "", colnames(views))
#show(views)

# Save the views object somewhere.

## PARALLELIZATION:
# After all the above work, you can easily parallelize the HMM over an index for the sample. For example:

param <- EmissionParam()
## In the following, j could be an index passed to Rscript in an array job
snp_exp <- SnpExperiment(views[, 1])
fit2 <- hmm2(snp_exp, param)
saveRDS(file, file = "hmm_fit2.rds")

save(fit2, file = "hmm_fit2.RData")

# select_columns <- match(c("SNP.Name", "Allele1...AB", "Allele2...AB", "Log.R.Ratio", "B.Allele.Freq"), names(X))
# index_genome   <- match(names(fgr), X[["SNP.Name"]])
# scan_params <- CopyNumScanParams(
# 	index_genome = index_genome, 
# 	select = select_columns, 
# 	cnvar = "Log.R.Ratio",
# 	bafvar = "B.Allele.Freq",
# 	gtvar = c("Allele1...AB", "Allele2...AB")
# 	)
}





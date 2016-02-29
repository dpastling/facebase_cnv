
options(stringsAsFactors = FALSE)

library(VanillaICE)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(foreach)
library(doSNOW)

# use this to supress the parallelization of VanillaICE
# registerDoSEQ()

cl <- makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

output.folder <- "results/vanilla_ice"
annotation.file <- "data//Marker_Info_Files/HumanOmni25Exome-8v1_A.csv"
files <- list.files(
	path = "data/Spritz_release_genotype_files", 
	pattern = "FinalReport", 
	recursive = TRUE, 
	full.names = TRUE
	)


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


features <- read.csv(annotation.file, skip = 7)
features <- features[! features$Chr %in% c("0", "MT", "X", "XY", "Y"), ]
features <- features[! is.na(features$MapInfo), ]
features[, "IntensityOnly"] <- as.integer(! grepl("\\[", features$SNP))
fgr <- GRanges(
	paste0("chr", features$Chr), 
	IRanges(features$MapInfo, width = 1), 
	isSnp = features$IntensityOnly == 0
	)
fgr <- SnpGRanges(fgr)
names(fgr) <- features[["Name"]]
sl <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlevels(fgr) <- sl[sl %in% seqlevels(fgr)]
seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(fgr),]
fgr <- sort(fgr)

parsedDir <- tempdir()
views <- ArrayViews(rowRanges = fgr, sourcePaths = files, parsedPath = parsedDir)
lrrFile(views) <- file.path(parsedDir, basename(fileName(views, "lrr")))
views@bafFiles <- file.path(parsedDir, basename(fileName(views, "baf")))
views@gtFiles <- file.path(parsedDir, basename(fileName(views, "gt")))
colnames(views) <- gsub(".csv", "", colnames(views))

dat <- fread(files[1])
select_columns <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB", "Log R Ratio", "B Allele Freq"), names(dat))
index_genome <- match(names(fgr), dat[["SNP Name"]])


scan_params <- CopyNumScanParams(
	index_genome = index_genome, 
	select = select_columns, 
	cnvar = "Log R Ratio",
	bafvar = "B Allele Freq",
	gtvar = c("Allele1 - AB", "Allele2 - AB")
	)
parseSourceFile(views, scan_params)
snp_exp <- SnpExperiment(views)

param <- EmissionParam()
fit <- hmm2(snp_exp, param)

filter_param <- FilterParam(numberFeatures = 5, probability = 0.95)

for (i in 1:length(fit))
{
	result <- cnvSegs(fit[[i]], filter_param)
	result <- as.data.frame(result)
	sample.name <- names(fit)[i]
	sample.name <- gsub("^.+?_([^_]+).csv$", "\\1", sample.name)
	write.table(result, file = paste0(output.folder, "/VICE_view_test2_", sample.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

stopCluster(cl)


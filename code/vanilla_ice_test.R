
options(stringsAsFactors = FALSE)

library(VanillaICE)
#library(BSgenome.Hsapiens.UCSC.hg18)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(foreach)
#library(doSNOW)
library(data.table)

registerDoSEQ()

#cl <- makeCluster(12, type = "SOCK")
#registerDoSNOW(cl)

#extdir <- system.file("extdata", package = "VanillaICE", mustWork = TRUE)
#data.files  <- list.files(extdir, pattern = "FinalReport", full.names = TRUE)

extdir <- "./"
annotation.file <- list.files(extdir, pattern = "SNP_info", full.names = TRUE)
data.files <- list.files(extdir, pattern = "sample_", full.names = TRUE)

features <- read.csv(annotation.file)

#file.index <- 2

##########################################################
#
# Option 1: Load first sample and process individually
#
##########################################################

fit_opt1 <- list()
for(file.index in 1:length(data.files))
{
#X <- read.csv(data.files[file.index], header = TRUE, skip = 9)
X <- read.csv(data.files[file.index], header = TRUE)
X <- inner_join(X, features, by = c("SNP.Name" = "Name"))
X <- X %>% arrange(as.numeric(Chr), Position)

fgr <- GRanges(
    paste0("chr", X$Chr), 
    IRanges(X$Position, width = 1), 
    isSnp = X$Intensity.Only == 0
    )
fgr            <- SnpGRanges(fgr)
names(fgr)     <- X[["SNP.Name"]]
sl             <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlevels(fgr) <- sl[sl %in% seqlevels(fgr)]
seqinfo(fgr)   <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(fgr), ]
fgr            <- sort(fgr)

temp <- X[, c("Log.R.Ratio", "B.Allele.Freq")]
temp <- as.matrix(temp)

snp_expt <- SnpArrayExperiment(
                cn = temp[, "Log.R.Ratio", drop = FALSE], 
                baf = temp[, "B.Allele.Freq", drop = FALSE], 
                rowRanges = fgr
            )
param    <- EmissionParam()
fit1     <- hmm2(snp_expt, param)
fit_opt1[[file.index]] <- fit1
}

##########################################################
#
# Option 2: Load samples together
#
##########################################################

parsedDir <- tempdir()
views <- ArrayViews(rowRanges = fgr, sourcePaths = data.files, parsedPath = parsedDir)
lrrFile(views) <- file.path(parsedDir, basename(fileName(views, "lrr")))
views@bafFiles <- file.path(parsedDir, basename(fileName(views, "baf")))
views@gtFiles <- file.path(parsedDir, basename(fileName(views, "gt")))
colnames(views) <- gsub(".txt", "", colnames(views))

dat <- fread(data.files[1])
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
fit2 <- hmm2(snp_exp, param)


##########################################################
#
# Compare results
#
##########################################################

#filter_param <- FilterParam(numberFeatures = 5, probability = 0.95)
filter_param <- FilterParam()

for ( file.index in 1:length(fit_opt1))
{
    result1 <- cnvSegs(fit_opt1[[file.index]], filter_param)
    result1 <- as.data.frame(result1)
    result1 <- result1[, 1:8]
    
    result2 <- cnvSegs(fit2[[file.index]], filter_param)
    result2 <- as.data.frame(result2)
    
    #stopCluster(cl)
    
    results.match <- identical(result1, result2)

    cat("testing file ", file.index, "\n")    
    if (results.match)
    {
        cat("The results are the same!\n")
    } else {
        cat("The results do not match\n")
        print(result1)
        print(result2)
    }
}





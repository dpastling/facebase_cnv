
library(oligoClasses)
library(VanillaICE)
library(crlmm)
library(SNPchip)
library(IRanges)
library(foreach)
library(ff)

#library(snow)
#library(doSNOW)
#cl <- makeCluster(2, type = "SOCK")
#registerDoSNOW(cl)


load("cache/cnSet_final.RData")

se <- as(cnSet, "SnpArrayExperiment")

library(ArrayTV)
i <- seq_len(ncol(se))
increms <- c(10, 1000, 100e3)
wins <- c(100, 10e3, 1e6)
res <- gcCorrect(
	lrr(se), 
	increms = increms, 
	maxwins = wins, 
	returnOnlyTV = FALSE, 
	verbose = TRUE, 
	build = "hg19", 
	chr = chromosome(se), 
	starts = start(se)
)

se2 <- se
assays(se2)[["cn"]] <- res$correctedVals


res <- hmm2(se2)


save(res, file = "cache/hmm_fit.RData")



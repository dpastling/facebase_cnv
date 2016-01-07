options(stringsAsFactors = FALSE)

library(crlmm)
#library(ff)
library(humanomni25quadv1bCrlmm)
library(foreach)
library(snow)
library(doSNOW)
library(lattice)
library(parallel)

#options(ffcaching = "ffeachflush")
#outdir <- paste("crlmm/", getRversion(), "/illumina_vignette", sep = "")
#ldPath(outdir)
#dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# We will also store cached computations in the directory outdir.
# We declare that crlmm should process 150,000 markers at a time and/or 500 samples at a time (when possible) to reduce the memory footprint. As our example dataset in this vignette contains fewer than 500 samples, all samples will be processed simultaneously.
ocProbesets(150e3)
ocSamples(3)

datadir <- "data/Array/Illumina HumanOmni2.5 BeadChip"
samplesheet = read.csv(file.path(datadir, "SampleSheet.csv"), skip = 8, header = TRUE, as.is = TRUE)
arrayNames <- file.path(datadir, paste(samplesheet[, "SentrixBarcode_A"], samplesheet[, "Sample_ID"], sep = "/"))

arrayInfo <- list(barcode = "SentrixBarcode_A", position = "SentrixPosition_A")

all(file.exists(paste(arrayNames, "_Grn.idat", sep="")))
all(file.exists(paste(arrayNames, "_Red.idat", sep="")))

## almost works = humanomniexpress12v1b

## get complete list with validCdfNames()
# cdfName <- "humanomniexpexome8v1p1b"
# cdfName <- "humanomni258v1a"
# cdfName <- "humanomni258v1p1b"
cdfName <- "humanomni25quadv1b"
# cdfName <- "humanomniexpress12v1b"  ## This works, but is the wrong chip
batch <- rep("1", nrow(samplesheet))

cnSet <- genotype.Illumina(
	sampleSheet = samplesheet,
	arrayNames = arrayNames,
	arrayInfoColNames = arrayInfo,
	cdfName = cdfName,
	batch = batch,
	call.method = 'krlmm'
)

save(cnSet, file = "cache/cnSet.RData")

## Quality Control


invisible(open(cnSet$SNR))
snr <- cnSet$SNR[]
close(cnSet$SNR)

pdf(file = "plots/crlmm_histogram.pdf")
print(histogram(~snr,
	panel=function(...){
		panel.histogram(...)
		},
	breaks = 25, 
	xlim = range(snr), 
	xlab = "SNR")
)
dev.off()

## Copy Number

crlmmCopynumber(cnSet)

tmp <- totalCopynumber(cnSet, i=seq_len(nrow(cnSet)))

save(cnSet, file = "cache/cnSet_final.RData")
save(tmp, file = "cache/tmp.RData")


## Downstream Analysis
#library(oligoClasses)
library(VanillaICE)
#library(SNPchip)
#library(IRanges)
library(ArrayTV)
cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)
ocSamples(2)

oligoList <- BafLrrSetList(cnSet)


### Wave correction
i <- seq_len(ncol(oligoList))
# Process 20 samples at a time
ocSamples(20)
increms <- c(10, 1000, 100e3)
wins <- c(100, 10e3, 1e6)

gc.result <- gcCorrect(oligoList, increms = increms, maxwins = wins, verbose = TRUE)

assays(oligoList)[["cn"]] <- gc.result$correctedVals

hmm.result <- hmm2(oligoList)

# tvScores <- gcCorrect(oligoList, increms = increms, maxwins = wins, verbose = TRUE)
# #tvScores <- gcCorrect(oligoList, increms = increms, maxwins = wins, returnOnlyTV = TRUE, verbose = TRUE)
# order(tvScores[[1]][, 2], decreasing = TRUE)
#
# # Next, we select a small window of 10 bp and a larger window of 10,000 bp and pre-compute the gc composition
# gc.matrix <- computeGC(oligoList, c(10, 10e3), c(10,10e3))
#
# oligoList2 <- gcCorrect(oligoList, increms = c(10, 10e3), maxwins = c(10, 10e3), providedGC = gc.matrix)


## HMM
res <- hmm2(oligoList, p.hom = 0.1, nupdates = 5, TAUP = 1e10)


gr <- unlist(res)
chr <- paste("chr", chromosome(oligoList), sep="")
brList <- oligoList[chr %in% chromosome(gr)]
brList <- brList[, match(sampleNames(gr)[1], sampleNames(brList))]
se <- as(brList, "SummarizedExperiment")
df <- dataFrame(gr, se, maxgap=500e3)



stopCluster(cl)






colors <- c("red", "orange", "white", "white", "lightblue", "blue")[state(gr)]
figs <- latticeFigs(gr, df, colors=colors)









bar2 <- NULL
for (i in 1:length(bar))
{
    bar2 <- rbind(bar2, data.frame(bar[[i]][,], chr = i))
}



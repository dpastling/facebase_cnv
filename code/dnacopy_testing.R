
options(stringsAsFactor = FALSE)
library(dplyr)

options(stringsAsFactors = FALSE)

library(DNAcopy)

load("results/CIDR_18505.RData")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
param.index <- args[1]

output.folder <- "results/dnacopy_optimize"

# test default settings
# param.index <- 165

######################## 
# Generate parameters
######################## 
alpha       <- c(0.001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2)
p.method    <- c("hybrid", "perm")
nperm       <- c(10000, 50000, 100000)
eta         <- c(0.001, 0.025, 0.05, 0.075, 0.1)
undo.splits <- c("none", "sdundo")
parameters  <- expand.grid(alpha, eta, nperm, undo.splits, p.method, stringsAsFactors = FALSE)
colnames(parameters) <- c("alpha", "eta", "nperm", "undo.splits", "p.method")


######################## 
# Find segments
######################## 

segment_with_params <- function(
	CNA.object,
	alpha = parameters[param.index, "alpha"], 
	p.method = parameters[param.index, "p.method"], 
	nperm = parameters[param.index, "nperm"], 
	eta = parameters[param.index, "eta"], 
	undo.splits = parameters[param.index, "undo.splits"]
) {
	cna.segments <- segment(
		CNA.object, alpha = alpha, p.method = p.method, 
		nperm = nperm, eta = eta, undo.splits = undo.splits
	)
	return(cna.segments$output)
}

segments1 <- segment_with_params(sample1)
segments2 <- segment_with_params(sample2)
segments3 <- segment_with_params(sample3)


bed_intersect <- function(bed1, bed2, inverse = FALSE)
{
        matches <- NULL
        for (i in 1:nrow(bed1))
        {
                chromosomes.match <- bed1[i, 1] == bed2[, 1]
                pos1.overlap <- bed1[i, 2] >= bed2[, 2] & bed1[i, 2] <= bed2[, 3]
                pos2.overlap <- bed1[i, 3] >= bed2[, 2] & bed1[i, 3] <= bed2[, 3]
                full.overlap <- bed1[i, 2] <= bed2[, 2] & bed1[i, 3] >= bed2[, 3]
				logR.matches <- sign(bed1[i, 5]) == sign(bed2[, 5])

                result <- any(chromosomes.match & (pos1.overlap | pos2.overlap | full.overlap) & logR.matches)
                matches <- c(matches, result)
        }
        if (inverse == TRUE)
        {
                matches <- ! matches
        }
        return(bed1[matches, ])
}

# For testing
# segments1 <- read.delim("results/dnacopy_alpha0_1/segment_alpha0_1_CIDR_18505.1054736194.txt")
# segments2 <- read.delim("results/dnacopy_alpha0_1/segment_alpha0_1_CIDR_18505.1054736220.txt")
# segments3 <- read.delim("results/dnacopy_alpha0_1/segment_alpha0_1_CIDR_18505.1064692236.txt")

filter_dnacopy <- function(df)
{
	df <- df %>%
		  select(-ID) %>%
		  filter(seg.mean >= 0.1 | seg.mean <= -0.2) %>%
		  filter(! chrom %in% c("MT", "X", "Y", "XY")) %>%
		  filter(num.mark >= 5) %>%
		  filter(loc.end - loc.start < 1e7)  
	return(df)
}

segments1 <- filter_dnacopy(segments1)
segments2 <- filter_dnacopy(segments2)
segments3 <- filter_dnacopy(segments3)


############################################
# Find Overlaps Among DNAcopy Replicates
############################################
dna.1_v_2 <- bed_intersect(segments1, segments2)
dna.all   <- bed_intersect(dna.1_v_2, segments3)

dna.inverse.1 <- bed_intersect(segments1, segments2, inverse = TRUE)
dna.inverse.1 <- bed_intersect(dna.inverse.1, segments3, inverse = TRUE)

dna.inverse.2 <- bed_intersect(segments2, segments1, inverse = TRUE)
dna.inverse.2 <- bed_intersect(dna.inverse.2, segments3, inverse = TRUE)

dna.inverse.3 <- bed_intersect(segments3, segments1, inverse = TRUE)
dna.inverse.3 <- bed_intersect(dna.inverse.3, segments2, inverse = TRUE)

FP.dna <- bind_rows(dna.inverse.1, dna.inverse.2, dna.inverse.3)


############################################
# Find Overlaps with HapMap Gold Standard
############################################
gold.standard <- read.delim("data/validation/estd20_NA18505_5plus_probes.bed", header = FALSE)
colnames(gold.standard) <- c("chr", "pos1", "pos2", "Name", "CNV", "n.markers")
gold.standard <- mutate(gold.standard, CNV = ifelse(grepl("gain", CNV), +1, -1))
gold.standard <- select(gold.standard, chr, pos1, pos2, n.markers, CNV)
dna.gold <- bed_intersect(dna.all, gold.standard)

TP.gold <- bed_intersect(dna.all, gold.standard)
FP.gold <- bed_intersect(dna.all, gold.standard, inverse = TRUE)


############################################
# Find Overlaps with PennCNV
############################################
penn1 <- read.delim("results/penncnv/sample_CIDR_18505@1054736194.tabcnv", header = FALSE)
penn2 <- read.delim("results/penncnv/sample_CIDR_18505@1054736220.tabcnv", header = FALSE)
penn3 <- read.delim("results/penncnv/sample_CIDR_18505@1064692236.tabcnv", header = FALSE)

penn1 <- penn1[, c(1,2,3,9,4)]
penn2 <- penn2[, c(1,2,3,9,4)]
penn3 <- penn3[, c(1,2,3,9,4)]

colnames(penn1) <- c("chr", "pos1", "pos2", "n.markers", "CNV")
colnames(penn2) <- c("chr", "pos1", "pos2", "n.markers", "CNV")
colnames(penn3) <- c("chr", "pos1", "pos2", "n.markers", "CNV")

# transform CNV call from 0, 1, 3, 4 into positive and negative numbers
penn1 <- mutate(penn1, CNV = (CNV / 2) - 1)
penn2 <- mutate(penn2, CNV = (CNV / 2) - 1)
penn3 <- mutate(penn3, CNV = (CNV / 2) - 1)

penn.1_v_2 <- bed_intersect(penn1, penn2)
penn.all   <- bed_intersect(penn.1_v_2, penn3)

TP.penn <- bed_intersect(dna.all, penn.all)
FP.penn <- bed_intersect(dna.all, penn.all, inverse = TRUE)


############################################
# Evaluate performance of each sample
# In practice we want each sample to have
# a high signal to noise ratio
############################################
TP.penn.1 <- bed_intersect(segments1, penn.all)
FP.penn.1 <- bed_intersect(segments1, penn.all, inverse = TRUE)
TP.penn.2 <- bed_intersect(segments2, penn.all)
FP.penn.2 <- bed_intersect(segments2, penn.all, inverse = TRUE)
TP.penn.3 <- bed_intersect(segments3, penn.all)
FP.penn.3 <- bed_intersect(segments3, penn.all, inverse = TRUE)


TP.gold.1 <- bed_intersect(segments1, gold.standard)
FP.gold.1 <- bed_intersect(segments1, gold.standard, inverse = TRUE)
TP.gold.2 <- bed_intersect(segments2, gold.standard)
FP.gold.2 <- bed_intersect(segments2, gold.standard, inverse = TRUE)
TP.gold.3 <- bed_intersect(segments3, gold.standard)
FP.gold.3 <- bed_intersect(segments3, gold.standard, inverse = TRUE)

TP.dna.1 <- bed_intersect(segments1, dna.all)
FP.dna.1 <- bed_intersect(segments1, dna.all, inverse = TRUE)
TP.dna.2 <- bed_intersect(segments2, dna.all)
FP.dna.2 <- bed_intersect(segments2, dna.all, inverse = TRUE)
TP.dna.3 <- bed_intersect(segments3, dna.all)
FP.dna.3 <- bed_intersect(segments3, dna.all, inverse = TRUE)



############################################
# Report Results
############################################
size.1 <- segments1[,3] - segments1[,2]
size.2 <- segments2[,3] - segments2[,2]
size.3 <- segments3[,3] - segments3[,2]
mean.size <- mean(c(size.1, size.2, size.3))

result <- data.frame(
	parameters[param.index, ], 
	mean.size = mean.size,
	TP.replicates = nrow(dna.all),
	FP.replicates = nrow(FP.dna),
	TP.gold = nrow(TP.gold),
	FP.gold = nrow(FP.gold),
	TP.penn = nrow(TP.penn),
	FP.penn = nrow(FP.penn),
	TP.penn.1 = nrow(TP.penn.1),
	FP.penn.1 = nrow(FP.penn.1),
	TP.penn.2 = nrow(TP.penn.2),
	FP.penn.2 = nrow(FP.penn.2),
	TP.penn.3 = nrow(TP.penn.3),
	FP.penn.3 = nrow(FP.penn.3),
	TP.gold.1 = nrow(TP.gold.1),
	FP.gold.1 = nrow(FP.gold.1),
	TP.gold.2 = nrow(TP.gold.2),
	FP.gold.2 = nrow(FP.gold.2),
	TP.gold.3 = nrow(TP.gold.3),
	FP.gold.3 = nrow(FP.gold.3),
	TP.dna.1 = nrow(TP.dna.1),
	FP.dna.1 = nrow(FP.dna.1),
	TP.dna.2 = nrow(TP.dna.2),
	FP.dna.2 = nrow(FP.dna.2),
	TP.dna.3 = nrow(TP.dna.3),
	FP.dna.3 = nrow(FP.dna.3)
)

cat("done!\n")

write.table(result, file = paste0(output.folder, "/parameter_", param.index, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)







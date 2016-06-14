options(stringsAsFactors = FALSE)

library(DNAcopy)

load("results/CIDR_18505.RData")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
param.index <- args[1]


######################## 
# Generate parameters
######################## 

alpha       <- c(0.001, 0.005, 0.01, 0.05, 0.075, 0.1, 0.15, 0.175, 0.2)
p.method    <- c("hybrid", "perm")
nperm       <- c(10000, 50000, 100000)
eta         <- c(0.001, 0.01, 0.2, 0.05, 0.1)
undo.splits <- c("none", "sdundo")
parameters  <- expand.grid(alpha, p.method, nperm, eta, undo.splits, stringsAsFactors = FALSE)
colnames(parameters) <- c("alpha", "p.method", "nperm", "eta", "undo.splits")


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


######################## 
# Intersect Bedss
######################## 

bed_intersect <- function(bed1, bed2, inverse = FALSE)
{
	matches <- NULL
	for (i in 1:nrow(bed1))
	{
		chromosomes.match <- bed1[i, 1] == bed2[, 1] 
		pos1.overlap <- bed1[i, 2] >= bed2[, 2] & bed1[i, 2] <= bed2[, 3]
		pos2.overlap <- bed1[i, 3] >= bed2[, 2] & bed1[i, 3] <= bed2[, 3]
		full.overlap <- bed1[i, 2] <= bed2[, 2] & bed1[i, 3] >= bed2[, 3]

		result <- any(chromosomes.match & (pos1.overlap | pos2.overlap | full.overlap))
		matches <- c(matches, result)
	}
	if (inverse == TRUE)
	{
		matches <- ! matches
	}
	return(bed1[matches, ])
}

# filter by coverage


# intersect 1 and 2
# intersect 1 and 3
# intersect 1_2 versus 1_3

# reverse intersect 1 and 2
# reverse intersect r_1_2 and 3
# reverse intersect 2 and 1
# reverse intersect r_2_1 and 3
# reverse intersect 3 and 1
# reverse intersect r_3_1 and 2

# combine r_1_2_3, r_2_1_3, and r_3_1_2

# count intersect
# count reverse intersect
# save parameters, numbers




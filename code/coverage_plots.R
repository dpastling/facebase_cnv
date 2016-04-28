
options(stringsAsFactors = FALSE)

library(dplyr)

df <- read.delim("data/corrected_vals/corrected_vals_CIDR_18505@1054736194.txt")

df <- select(df, Chr, MapInfo, Log.R.Ratio, corrected.vals, SNP.Name)

plot_segment <- function(chr, x1, x2, x.min, x.max)
{
	temp <- filter(df, Chr == chr, MapInfo > x.min, MapInfo < x.max)
	with(temp, plot(MapInfo, corrected.vals, ylim = c(-1, 1), xlab = "Genome Coordinate", ylab = "GC Corrected LogR", main = paste0("chr", chr, ":", x1, "-", x2)))
	mean.logR <- filter(df, Chr == chr, MapInfo >= x1, MapInfo <= x2)
	mean.logR <- mean(mean.logR[["corrected.vals"]])
	segments(x1, mean.logR, x2, mean.logR, col = "red", lwd = 2)
}


pdf(file = "coverage_1_95134903_95154785.pdf")
plot_segment("1", 95134903, 95154785, 95000000, 95400000)
dev.off()

pdf(file = "coverage_1_95134903_95154785_zoom.pdf")
plot_segment("1", 95134903, 95154785, 94500000, 95780000)
dev.off()

pdf(file = "coverage_1_196823300_196901753.pdf")
plot_segment("1", 196823300, 196901753, 196500000, 197200000)
dev.off()

pdf(file = "coverage_5_101169234_101209795.pdf")
plot_segment("5", 101169234, 101209795, 101050000, 101300000)
dev.off()

pdf(file = "coverage_8_13615303_13649956.pdf")
plot_segment("8", 13615303, 13649956, 13400000, 13865303)
dev.off()

pdf(file = "coverage_1_99577011_103473443.pdf")
plot_segment("1", 99577011, 103473443, 9.5e+07, 1.08e+08)
dev.off()

pdf(file = "coverage_3_75913757_90311605.pdf")
plot_segment("3", 75913757, 90311605, 70913757, 95311605)
dev.off()

pdf(file = "coverage_10_135242873_135379710.pdf")
plot_segment("10", 135242873, 135379710, 135e6, 135.6e6)
dev.off()


## Big Deletion
pdf(file = "coverage_2_89978069_90256530.pdf")
plot_segment("2", 89978069, 90256530, 80e6, 120e6)
dev.off()

# Normal Deletion
pdf(file = "coverage_3_62665222_62694321.pdf")
plot_segment("3", 62665222, 62694321, 62.1e6, 63e6)
dev.off()

# duplications
pdf(file = "coverage_10_135242873_135369532.pdf")
plot_segment("10", 135242873, 135369532, 134e6, 137e6)
dev.off()

pdf(file = "coverage_4_168029135_168080180.pdf")
plot_segment("4", 168029135, 168080180, 167.5e6, 168.5e6)
dev.off()





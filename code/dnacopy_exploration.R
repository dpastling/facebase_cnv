
options(stringsAsFactors = FALSE)

library(dplyr)

df <- read.delim("results/logR/dnacopy/CIDR_18505@1054736194.txt")

#files <- list.files(path = "results/logR/dnacopy", pattern = ".txt", full.names = TRUE)
#df <- lapply(files, read.delim)
#df <- bind_rows(df)

df <- filter(df, num.mark >= 5, ! chrom %in% c("MT", "X", "Y", "XY")) 
df <- mutate(df, length = loc.end - loc.start)
point.size <- (df[["length"]] / max(df[["length"]])) * 5

pdf("plots/exploratory/dnacopy_logR.pdf")
with(df, plot(log10(pval), mean.norm.logR, ylim = c(-1,1), xlim = c(-100, 0)))
dev.off()

pdf("plots/exploratory/dnacopy_hist_logR.pdf")
temp <- filter(df, mean.norm.logR >= -1, mean.norm.logR <= 1)
hist(temp[["mean.norm.logR"]], breaks = 10, xlim = c(-1, 1), xlab = "GC Corrected LogR", main = "Sample CIDR_18505@1054736194")
dev.off()

#pdf("plots/exploratory/dnacopy_hist_logR.pdf")
#temp <- filter(df, mean.norm.logR >= -1, mean.norm.logR <= 1)
#hist(temp[["mean.norm.logR"]], xlab = "GC Corrected LogR", main = "Combined HapMap Samples")
#dev.off()

pdf("plots/exploratory/dnacopy_len_vs_logR.pdf")
with(df, plot(log10(length), mean.norm.logR, ylim = c(-1, 1)))
dev.off()

############################
# intersect with replicates
############################

df_1_vs_2 <- read.delim("results/vanilla_ice_qc/18505_1054736194_1054736220_dna.bed", header = FALSE)
df_1_vs_3 <- read.delim("results/vanilla_ice_qc/18505_1054736194_1064692236_dna.bed", header = FALSE)

df_rep <- mutate(df, cnv.id = paste0(chrom, ":", loc.start, "-", loc.end))
df_1_vs_2 <- paste0(df_1_vs_2[, 1], ":", df_1_vs_2[, 2], "-", df_1_vs_2[, 3])
df_1_vs_3 <- paste0(df_1_vs_3[, 1], ":", df_1_vs_3[, 2], "-", df_1_vs_3[, 3])

#df_rep <- df_rep %>% filter(cnv.id %in% df_1_vs_2 & cnv.id %>% df_1_vs_3)
df_rep = df_rep[df_rep[["cnv.id"]] %in% df_1_vs_2 & df_rep[["cnv.id"]] %in% df_1_vs_3, ]


pdf("plots/exploratory/dnacopy_logR_replicates.pdf")
with(df, plot(log10(pval), mean.norm.logR, ylim = c(-1,1), xlim = c(-100, 0)))
with(df_rep, points(log10(pval), mean.norm.logR, col = "red"))
dev.off()


penn <- read.delim("results/vanilla_ice_qc/18505.1054736194_dna_penn_overlap.txt", header = FALSE)
df_penn <- mutate(df, cnv.id = paste0(chrom, ":", loc.start, "-", loc.end))
penn <- paste0(penn[, 1], ":", penn[, 2], "-", penn[, 3])
df_penn <- df_penn[df_penn[["cnv.id"]] %in% penn, ]
df_vice <- df_vice[df_penn[["length"]] < 100e6, ]

pdf("plots/exploratory/dnacopy_logR_penn.pdf")
with(df, plot(log10(pval), mean.norm.logR, ylim = c(-1,1), xlim = c(-100, 0)))
with(df_penn, points(log10(pval), mean.norm.logR, col = "red"))
dev.off()

pdf("plots/exploratory/dnacopy_hist_logR_penn.pdf")
temp <- filter(df_penn, mean.norm.logR >= -1, mean.norm.logR <= 1)
hist(temp[["mean.norm.logR"]], breaks = 10, xlim = c(-1, 1), xlab = "GC Corrected LogR", main = "PennCNV Overlap\nSample CIDR_18505@1054736194")
dev.off()



vice <- read.delim("results/vanilla_ice_qc/18505.1054736194_dna_vice_overlap.txt", header = FALSE)
df_vice <- mutate(df, cnv.id = paste0(chrom, ":", loc.start, "-", loc.end))
vice <- paste0(vice[, 1], ":", vice[, 2], "-", vice[, 3])
df_vice <- df_vice[df_vice[["cnv.id"]] %in% vice, ]
df_vice <- df_vice[df_vice[["length"]] < 100e6, ]

pdf("plots/exploratory/dnacopy_logR_vice.pdf")
with(df, plot(log10(pval), mean.norm.logR, ylim = c(-1,1), xlim = c(-100, 0)))
with(df_vice, points(log10(pval), mean.norm.logR, col = "red", pch = 16))
dev.off()

pdf("plots/exploratory/dnacopy_hist_logR_vice.pdf")
temp <- filter(df_vice, mean.norm.logR >= -1, mean.norm.logR <= 1)
hist(temp[["mean.norm.logR"]], breaks = 10, xlim = c(-1, 1), xlab = "GC Corrected LogR", main = "VanillaICE Overlap\nSample CIDR_18505@1054736194")
dev.off()





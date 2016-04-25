options(stringsAsFactors = FALSE)

library(dplyr)

df <- read.csv("HumanOmni25Exome-8v1_A.csv", skip = 7)

# the sequence data takes up a lot of memory
df <- select(df, -ends_with("Seq"))


temp <- df %>% filter(Chr == 1, MapInfo > 231e6, MapInfo < 233e6) %>% select(Chr, MapInfo)
pdf(file = "plots/map_locations.pdf", h = 4, w = 12)
plot(temp[["MapInfo"]], jitter(rep(1, nrow(temp)), amount = 0.2), xlab = "Genome Coordinate", ylab = "", yaxt = "n", ylim = c(0.5,1.5))
dev.off()

#d <- density(temp$MapInfo)
pdf("plots/denisty.pdf")
hist(temp$MapInfo, breaks = 100, xlab = "Genome Coordinate (in bp)", main = "Number of Probes per 10 kb")
dev.off()



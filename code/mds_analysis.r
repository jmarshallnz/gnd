#' MDS analyses
library(vegan)
library(RColorBrewer)
source("code/read_metadata.R")


pdf("figures/mds.pdf", width=10, height=8)

abund_samp = read.csv("no_error_abundance_with_ctrl.csv", row.names=1)
abund_samp = t(abund_samp)
abund_per_sample = abund_samp / rowSums(abund_samp)
abund_dist = as.matrix(dist(abund_per_sample))
meta <- read_metadata(animal_cols = rownames(abund_samp))

# table for things
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 22, 21), size=c(2, 2, 1, 1))
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 23, 24), size=c(1, 1, 1, 1))
cols <- read.csv("data/colours.csv", stringsAsFactors = FALSE) %>% arrange(animal_tag)

meta = meta %>% mutate(animal = as.numeric(animal)) %>% left_join(cols, by=c('animal' = 'animal_tag')) %>% left_join(shapes, by=c('source' = 'treatment'))

set.seed(5)


mds <- monoMDS(abund_dist,k=2)

x <- mds$points[,1]
y <- mds$points[,2]

plot(x, y, main="", type="n", xlab="Coordinate 1", ylab="Coordinate 2", xlim=c(min(x), max(x)+0.5))
points(x,y,
       col="black", pch=meta$shape, bg=meta$hex, cex=meta$size)
legend("topright", legend=c(cols$animal_tag, "ctrl"), pch=c(rep(19, 20), 1), col=c(cols$hex, "black"), cex=0.8, pt.cex=c(rep(0.8, 20), 2), bty="n")
legend("topleft", legend=shapes$treatment, pch=shapes$shape, col="grey30", cex=0.8, bty="n")
#legend(2, 2, legend=c("chicken", "cow", "sheep", "water"), fill=1:4)

wch = is.na(meta$animal)
points(x[wch], y[wch], cex=2, col="grey30")

dev.off()

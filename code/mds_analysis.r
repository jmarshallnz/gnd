#' MDS analyses
library(MASS)
library(RColorBrewer)
source("code/read_metadata.R")

abund_samp = read.csv("no_error_abundance_with_ctrl.csv", row.names=1)
abund_samp = t(abund_samp)
abund_per_sample = abund_samp / rowSums(abund_samp)
abund_dist = as.matrix(dist(abund_per_sample))
samp.meta = read_metadata(animal_cols = rownames(abund_samp))

set.seed(5)


mds <- monoMDS(abund_dist,k=2)

x <- mds$points[,1]
y <- mds$points[,2]

pdf("figures/mds.pdf", width=10, height=8)
cols = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
plot(x, y, main="", type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(min(x), max(x)+0.5))
points(x,y,col=cols[as.numeric(as.factor(samp.meta$animal))], pch = 16+as.numeric(as.factor(samp.meta$source)))
legend("topright", legend=levels(as.factor(samp.meta$animal)), pch=19, col=cols, cex=0.8)
legend("topleft", legend=levels(as.factor(samp.meta$source)), pch=16+1:4, col=cols[24], cex=0.8)
#legend(2, 2, legend=c("chicken", "cow", "sheep", "water"), fill=1:4)

wch = samp.meta$animal == "ctrl"
points(x[wch], y[wch], cex=2)

dev.off()

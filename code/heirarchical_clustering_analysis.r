#' code for doing heirarchical clustering of abundance data

library(readxl)
library(dplyr)
source("code/read_metadata.R")

abund_samp = read.csv("no_error_abundance_with_ctrl.csv", row.names=1)
abund_samp = t(abund_samp)
abund_per_sample = abund_samp / rowSums(abund_samp)
abund_dist = as.matrix(dist(abund_per_sample))


animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", rownames(abund_per_sample))
rep    <- sub("X([0-9]+)_.*", "\\1", rownames(abund_per_sample))
rownames(abund_dist) = colnames(abund_dist) = paste0(animal, "_", rep)

# add labels
rownames(abund_dist) = colnames(abund_dist) = read_metadata(animal_cols = rownames(abund_per_sample))$label

# do some clustering
samp.hc = hclust(as.dist(abund_dist))
pdf("sample_tree1.pdf", width=12, height=8)
plot(samp.hc, main="Sample clustering (method 1)", xlab="", sub="", cex=0.7)
dev.off()

samp.hc2 = hclust(as.dist(abund_dist), method="average")
pdf("sample_tree2.pdf", width=12, height=8)
plot(samp.hc2, main="Sample clustering (method 2)", xlab="", sub="", cex=0.7)
dev.off()

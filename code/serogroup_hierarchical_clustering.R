#' Serogroup clustering
#' 
#' This is just hierarchical clustering on serogroups, ignoring anything else.

library(cluster)
source("code/read_fasta.R")
source("code/read_abundance.R")

fa15 = read_fasta()

#' Now compute the distances between these things
fa15.dist = as.matrix(daisy(fa15 %>% select(starts_with('X')) %>% as.matrix %>% as.data.frame))
rownames(fa15.dist) = fa15$serogroup
colnames(fa15.dist) = fa15$serogroup

#' Reduce it down to those we're removing (do we want 80 or 84?)
y = rowSums(read_abundance(removed=c(97,98,12,'ctrl')))

fa15.dist = fa15.dist[names(y), names(y)]

#' Now do the clustering
fast.hc = hclust(as.dist(fa15.dist))
pdf("figures/gnd_clustering.pdf", width=12, height=8)
plot(fast.hc, main="Serogroup clustering", xlab="", sub="", cex=0.2)
dev.off()

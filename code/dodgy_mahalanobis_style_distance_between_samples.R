#' ---
#' title: GND data analysis
#' author: Jonathan Marshall
#' ---
#' 
#' Read the data in

source("code/read_metadata.R")
fa15.dist = as.matrix(read.csv("temp/sero_dist15.csv", row.names=1))
fa15.abund = read.csv("temp/sero_abundance.csv", row.names=1)

#' now, how can we add this up so we get differences between isolates?
#'
#' We'll need to use some sort of Mahalanobis style measure.
#' 
#' We use (n_i - n_j)' d (n_i - n_j)
#' 
#' where n_i, n_j are the number (proportion?) of reads of each serogroup across samples i and j.
#' (i.e. n_i, n_j are the vector of proportions), and d is a similarity (not distance!) matrix between
#' each of the serogroups (i.e. d = 1 on the diagonal, d < 1 off the diagonals).
#' 
#' We probably want proportions for the differences I guess, so we are normalised for the number of reads?
#'
#' First we want to grab the numbers...

samp_dist = matrix(NA, ncol(fa15.abund), ncol(fa15.abund))
for (i in 1:ncol(fa15.abund)) {
  n_i = fa15.abund[,i]/sum(fa15.abund[,i])
  for (j in 1:ncol(fa15d)) {
    n_j = fa15.abund[,j]/sum(fa15.abund[,j])
    d_ij = n_i - n_j
    samp_dist[i,j] = t(d_ij) %*% (1-fa15.dist) %*% d_ij
  }
}
samp_dist[1:5,1:5]

# add labels
rownames(samp_dist) = colnames(samp_dist) = read_metadata(animal_cols = colnames(fa15.abund))$label

# do some clustering
samp.hc = hclust(as.dist(samp_dist))
pdf("figures/sample_tree1_dodgy_jm_method.pdf", width=12, height=8)
plot(samp.hc, main="Sample clustering (Dodgy JM method 1)", xlab="", sub="", cex=0.7)
dev.off()

samp.hc2 = hclust(as.dist(samp_dist), method="average")
pdf("figures/sample_tree2_dodgy_jm_method.pdf", width=12, height=8)
plot(samp.hc2, main="Sample clustering (Dodgy JM method 2)", xlab="", sub="", cex=0.7)
dev.off()

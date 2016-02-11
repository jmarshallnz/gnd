#' ---
#' title: GND data analysis
#' author: Jonathan Marshall
#' ---
#' 
#' Read the data in

library(seqinr)
library(dplyr)
library(cluster)

fa15 = read.fasta("solexaQA15thou_nucleotideGE10.fa")
fa15meta = read.table("solexaQA15thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# now label the fa15 by serotroup rather than MD5, and merge to give id/sequence dataframe
fa_names <- substring(names(fa15),1,32)

#' Hmm, the following rows are missing stuff
rows <- which(!(fa15meta$md5 %in% fa_names))

#' reshape the fa data into a data frame
fa15 = data.frame(do.call(rbind, fa15))
fa15$md5 = substring(rownames(fa15), 1, 32)
rownames(fa15) <- NULL

#' join our id table
fa15 = fa15 %>% left_join(fa15meta %>% select(md5, serogroup))

#' Now compute the distances between these things
fa15.dist = as.matrix(daisy(fa15[,1:284]))
rownames(fa15.dist) = fa15$serogroup
colnames(fa15.dist) = fa15$serogroup

#' Do some clustering for fun
fa.hc1 <- hclust(as.dist(fa15.dist*284))
plot(fa.hc1, cex=0.2)

fa.hc2 <- hclust(as.dist(fa15.dist*284), method='average')
plot(fa.hc2, cex=0.2)

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

samp_dist = matrix(NA, ncol(fa15d), ncol(fa15d))
for (i in 1:ncol(fa15d)) {
  n_i = fa15d[,i]/sum(fa15d[,i])
  for (j in 1:ncol(fa15d)) {
    n_j = fa15d[,j]/sum(fa15d[,j])
    d_ij = n_i - n_j
    samp_dist[i,j] = t(d_ij) %*% (1-fa15.dist) %*% d_ij
  }
}
samp_dist[1:5,1:5]
animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", colnames(fa15d))
rep    <- sub("X([0-9]+)_.*", "\\1", colnames(fa15d))
rownames(samp_dist) = colnames(samp_dist) = paste0(animal, "_", rep)

#' Read in metadata
library(readxl)
meta = read_excel("gnd_seqs_metadata.xlsx", sheet=2)

samp.meta = data.frame(sample=substring(colnames(fa15d),2), stringsAsFactors = FALSE)

samp.meta = samp.meta %>% left_join(meta, by=c('sample'='Library Name')) %>% select(sample, treatment=Treatment, source=`Description [Optional] `)
samp.meta$treatment = factor(samp.meta$treatment)
levels(samp.meta$treatment) = c("bf", "ct")
samp.meta$source = factor(samp.meta$source)
levels(samp.meta$source) = c("fec", "pob", "por", "pre")
samp.meta = samp.meta %>% mutate(animal = sub(".*_([0-9ctrl]+)_.*", "\\1", sample),
                                 label = paste0(animal, "_", treatment, "_", source))

# relabel the ctrls
wch <- samp.meta$animal == "ctrl"
ctrl_labs = sub("([0-9ctrl]+)_.*", "\\1", samp.meta$sample)
samp.meta$label[wch] <- paste0("ctrl_", ctrl_labs[wch])

# add labels
rownames(samp_dist) = colnames(samp_dist) = samp.meta$label

# do some clustering
samp.hc = hclust(as.dist(samp_dist))
plot(samp.hc)

samp.hc2 = hclust(as.dist(samp_dist), method="average")
plot(samp.hc2)

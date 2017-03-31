#' Create distance matrix based on 284 bases, and create abundance matrix from metadata

library(cluster)
source("code/read_fasta.R")

fa15 = read_fasta()

#' Now compute the distances between these things
fa15.dist = as.matrix(round(daisy(fa15 %>% select(starts_with('X')) %>% as.matrix %>% as.data.frame)*284))
rownames(fa15.dist) = fa15$serogroup
colnames(fa15.dist) = fa15$serogroup
write.csv(fa15.dist, "temp/sero_dist15.csv", row.names=TRUE)

#' Create abundance matrix
fa15meta = read.table("data/gnd2/solexaQA15thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
fa15d = fa15 %>% select(serogroup) %>% left_join(fa15meta) %>% select(serogroup, starts_with('X'))
rownames(fa15d) = fa15d$serogroup
fa15d = fa15d %>% select(starts_with('X')) %>% as.matrix
write.csv(fa15d, "temp/sero_abundance.csv", row.names=TRUE)

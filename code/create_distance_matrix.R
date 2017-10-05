#' Create distance matrix based on 284 bases, and create abundance matrix from metadata

library(cluster)
library(tibble)
library(tidyr)

source("code/read_fasta.R")

qa <- 15

fa <- read_fasta(qa=qa)

#' Now compute the distances between these things
fa.dist = as.matrix(round(daisy(fa %>% select(starts_with('X')) %>% as.matrix %>% as.data.frame)*284))
rownames(fa.dist) = fa$serogroup
colnames(fa.dist) = fa$serogroup

# save it for fast loading, and as a .csv
save(fa.dist, file=paste0("temp/sero_dist",qa,".Rda"))
write.csv(fa.dist, paste0("temp/sero_dist",qa,".csv"), row.names=TRUE)

# reshape this into something we can use in dplyr etc.
fa.dist.df <- fa.dist %>% as.data.frame %>% tibble::rownames_to_column("gST") %>%
  gather(gST2, Distance, -gST)

save(fa.dist.df, file=paste0("temp/sero_dist",qa,"_df.Rda"))

#' Create abundance matrix
file_meta <- paste0("data/gnd2/solexaQA",qa,"thou_minTotalGE10.txt")
fa_meta = read.table(file_meta, header=TRUE, sep="\t", stringsAsFactors = FALSE)
fa_d = fa %>% select(serogroup) %>% left_join(fa_meta) %>% select(serogroup, starts_with('X'), starts_with('Ctrl'))
rownames(fa_d) = fa_d$serogroup
fa_d = fa_d %>% select(starts_with('X'), starts_with('Ctrl')) %>% as.matrix
write.csv(fa_d, paste0("temp/sero_abundance",qa,".csv"), row.names=TRUE)

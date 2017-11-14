#' Create distance matrix based on 284 bases, and create abundance matrix from metadata

library(cluster)
library(tibble)
library(tidyr)

source("code/read_fasta.R")

farm <- 1
qa <- 3
minTotal <- 10

fa <- read_fasta(farm=farm, qa=qa, minTotal=minTotal)

#' Chop down fa...
fa$sequence <- apply(fa[,-c(1:2)], 1, paste0, collapse='')

#' Close distances
#fa.close.df = calc_dist(fa$serogroup, fa$sequence, 4) # use 4 cores
#save(fa.close.df, file=paste0("temp/sero_close",qa,"_", minTotal, "_df.Rda"))

#' Now compute the distances between these things
#fa.dist = as.matrix(round(daisy(fa %>% select(starts_with('X')) %>% as.matrix %>% as.data.frame)*284))
#rownames(fa.dist) = fa$serogroup
#colnames(fa.dist) = fa$serogroup

# save it for fast loading, and as a .csv
#save(fa.dist, file=paste0("temp/sero_dist",qa,"_", minTotal, ".Rda"))
#write.csv(fa.dist, paste0("temp/sero_dist",qa,"_", minTotal, ".csv"), row.names=TRUE)

# reshape this into something we can use in dplyr etc.
#fa.dist.df <- fa.dist %>% as.data.frame %>% tibble::rownames_to_column("gST") %>%
#  gather(gST2, Distance, -gST)

#save(fa.dist.df, file=paste0("temp/sero_dist",qa,"_", minTotal, "_df.Rda"))

#' Create abundance matrix
base_dir <- paste0("data/gnd2/farm", farm)
file_meta <- paste0("solexaQA",qa,"thou_minTotalGE", minTotal, ".txt")
fa_meta = read.table(file.path(base_dir, file_meta), header=TRUE, sep="\t", stringsAsFactors = FALSE)
fa_d = fa %>% select(serogroup) %>% left_join(fa_meta) %>% select(serogroup, starts_with('X'), starts_with('Ctrl'))

out_dir <- paste0("temp/farm", farm)
rownames(fa_d) = fa_d$serogroup
fa_dm = fa_d %>% select(starts_with('X'), starts_with('Ctrl')) %>% as.matrix
out_file = paste0("sero_abundance",qa,"_", minTotal, ".csv")
write.csv(fa_dm, file.path(out_dir, out_file), row.names=TRUE)

#' Find likely parents. We need only bother with those that are high abundance
parent_min_abund <- 20
parents <- fa_d %>% gather(Library, Count, -serogroup) %>% group_by(serogroup) %>%
  summarize(Count=sum(Count)) %>%
  filter(Count >= parent_min_abund) %>%
  left_join(fa %>% select(serogroup, sequence))

#' Compute the distance between each sequence and the parents
Rcpp::sourceCpp('src/close_dist.cpp')
dist.close <- list()
system.time(for (i in 1:nrow(fa)) {
  dist = dist_to_parents(fa$sequence[i], parents$sequence)
  # technically we don't need *all* the distances.
  # We need only those within, say, 5, and the minimum one.
  wch <- which(dist <= 5)
  if (length(wch) == 0) {
    wch <- which.min(dist)
  }
  # OK, now construct the distance dataframe
  dist.close[[i]] <- data.frame(gST = fa$serogroup[i], gST2 = parents$serogroup[wch], Distance = dist[wch], stringsAsFactors = FALSE)
  if (i %% 1000 == 0) {
    cat("Up to i=", i, "of", nrow(fa), "\n")
  }
})
fa.close.df <- do.call(rbind, dist.close)
out_file <- paste0("sero_close",qa,"_", minTotal, "_df.Rda")
save(fa.close.df, file=file.path(out_dir, out_file))


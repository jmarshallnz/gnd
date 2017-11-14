# test mapping of distances
library(Rcpp)

library(cluster)
library(tibble)
library(tidyr)

source("code/read_fasta.R")

qa <- 3
minTotal <- 0

fa <- read_fasta(qa=qa, minTotal=minTotal)

#' Chop down fa...
fa$sequence <- apply(fa[,-c(1:2)], 1, paste0, collapse='')

dist_to_parent <- function(names, sequences, parent_names, parent_sequences, maxdist=5) {
  # iterate over the sequences. Find the closest among the parents.
  for (i in 1:length(sequences)) {
    dist = dist_to_parents(sequences[i], parent_sequences)
  }
}


calc_dist <- function(names, dna, no_cores = 1) {
  # Formula for converting from linear x through a distance matrix to i,j column, row
  # i <- trunc(n+1/2 - sqrt(n^2 - n - 2*x + 0.5))
  # j <- x - (i-1)*n + (i-1)*i/2 + i

  no_cores <- min(parallel::detectCores(), no_cores)

  # internal function that does the work in chunks
  calc_dist_internal <- function(chain, dna, i_start, i_end) {
    Rcpp::sourceCpp('src/close_dist.cpp')
    data.frame(calcDist(dna, i_start[chain], i_end[chain], progress=(chain==1)))
  }

  if (no_cores > 1) {
    # break up the linear distance array, then convert to row/columns
    # Alternatively we could just use clusterApplyLB() though that wouldn't have any feedback
    x <- round(seq(1, n*(n-1)/2, length.out=no_cores+1))
    i <- trunc(n+1/2 - sqrt(n^2 - n - 2*x + 0.5))
    i_start <- i[-length(i)]
    i_end <- i[-1]-1; i_end[length(i_end)] <- n-1

    ## run the parallel version
    cluster <- parallel::makeCluster(no_cores, outfile="")
    out <- parallel::clusterApply(cluster, 1:no_cores, calc_dist_internal, dna=dna, i_start=i_start, i_end=i_end)
    parallel::stopCluster(cluster)
    out <- do.call(rbind, out)
  } else {
    out <- calc_dist_internal(1, dna, i_start = 1, i_end = length(dna)-1)
  }
  # Extend to duplicate results so we have faster downstream processing
  out2 <- data.frame(gST=out$gST2, gST2=out$gST, Distance=out$Distance)

  # And replace all gST/gST2 with the names
  final <- rbind(out, out2)
  final$gST <- names[final$gST+1]
  final$gST2 <- names[final$gST2+1]

  # Done!
  final
}

calc_single_dist <- function(seq) {
  dna1 <- fa[fa$serogroup == seq,3:286]
  dist <- numeric(nrow(fa))
  for (i in 1:nrow(fa)) {
    if (i %% 1000 == 0)
      cat("Up to ", i, " of ", nrow(fa), "\n")
    dna2 = fa[i,3:286]
    dist[i] = sum(dna1 != dna2)
  }
  data.frame(gST2=fa$serogroup, Distance=dist)
}


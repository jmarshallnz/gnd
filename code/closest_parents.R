Rcpp::sourceCpp('src/close_dist.cpp')

closest_parents <- function(id, sequence, parent_id, parent_sequence, min_dist = 5) {
  dist.close <- list()
  system.time(for (i in 1:length(id)) {
    dist = dist_to_parents(sequence[i], parent_sequence)
    # technically we don't need *all* the distances.
    # We need only those within, say, 5, and the minimum one.
    wch <- which(dist <= min_dist)
    if (length(wch) == 0) {
      wch <- which.min(dist)
    }
    # OK, now construct the distance dataframe
    dist.close[[i]] <- data.frame(id = id[i], id2 = parent_id[wch], Distance = dist[wch], stringsAsFactors = FALSE)
    if (i %% 1000 == 0) {
      cat("Up to i=", i, "of", length(id), "\n")
    }
  })
  do.call(rbind, dist.close)
}

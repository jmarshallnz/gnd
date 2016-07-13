#' Code for reading in abundance information, potentially ignoring several animals
read_abundance <- function(file="sero_abundance.csv", removed=c(97,98,120,"ctrl")) {
  y = read.csv(file, row.names=1)
  
  # filter out the animals we don't need
  if (length(removed) > 0) {
    animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", names(y))
    y = y[,!(animal %in% as.character(removed))]
  }
  
  #' now filter those out who have less than 10 for consistency (previous
  #' GSTs were also filtered by this criteria)
  GSTs <- rowSums(y) < 10
  y[!GSTs,]
}

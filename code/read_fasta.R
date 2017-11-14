library(seqinr)
library(dplyr)
library(digest)

read_fasta <- function(farm=1, qa=15, minTotal=10) {

  base_dir <- "data/gnd2"
  farm_dir <- paste0("farm", farm)
  file_meta <- paste0("solexaQA",qa,"thou_minTotalGE", minTotal, ".txt")
  # Read in the metadata (it has the abundances) and remove the non-functional groups
  fa_meta = read.table(file.path(base_dir, farm_dir, file_meta), header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(functional != 'no')

  # read in the fasta file
  file_fa <- paste0("solexaQA",qa,"thou_nucleotideGE", minTotal, ".fa")
  fa = read.fasta(file.path(base_dir, farm_dir, file_fa))

  # The fasta file does not include the O serogroups for some reason, so add them in as well
  fa_O = read.fasta(file.path(base_dir, "gnd_DB_09012017.fas"))
  # remap to use MD5
  seq <- unlist(lapply(fa_O, function(x) { toupper(paste(x,collapse='')) }))
  md5 <- unlist(lapply(seq, function(x) { digest(x, serialize=FALSE)}))
  names(fa_O) <- md5
  fa_O <- fa_O[!duplicated(md5)]

  # FOR ADRIAN:
#  known_serogroups <- data.frame(serogroup=names(md5), md5=md5, sequence=seq)
#  write.csv(known_serogroups, "known_serogroups_MD5.csv", row.names=FALSE)
  # FOR ADRIAN

  # combine together
  fa <- c(fa, fa_O)

  # pull out the MD5s so we can map to serogroup in the metadata file
  fa_md5 <- substring(names(fa),1,32)

  # check we have sequences for all the gSTs

  # we should have no sequences that aren't in the metadata file
  if (sum(!fa_meta$md5 %in% fa_md5) != 0) {
    print(which(!fa_meta$md5 %in% fa_md5))
    stop("Error: We have MD5's in the minTotal file that we don't have sequences for!")
  }

  #' HMM, SEEMS TO BE AN ISSUE
  wch <- which(lengths(fa) != 284)
  if (length(wch) > 0) {
    warning(paste("There seems to be", length(wch), "sequences in the fasta file that are the wrong length"))
    fa <- fa[-wch]
  }
  #' reshape the fa data into a data frame
  fa = data.frame(do.call(rbind, fa), stringsAsFactors = FALSE)
  fa$md5 = substring(rownames(fa), 1, 32)
  rownames(fa) <- NULL

  #' join our id table. NOTE use of inner_join here: We only want stuff we have both abundance and distance
  #' information
  final <- fa_meta %>% select(md5, serogroup) %>% inner_join(fa)

  final
}

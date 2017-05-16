library(seqinr)
library(dplyr)

read_fasta <- function() {

  # Read in the metadata (it has the abundances) and remove the non-functional groups
  fa_meta = read.table("data/gnd2/solexaQA15thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(functional != 'no')

  # read in the fasta file
  fa = read.fasta("data/gnd2/solexaQA15thou_nucleotideGE10.fa")

  # The fasta file does not include the O serogroups for some reason, so add them in as well
  fa_O = read.fasta("data/gnd2/gnd_DB_09012017.fas")
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

  #' reshape the fa data into a data frame
  fa = data.frame(do.call(rbind, fa), stringsAsFactors = FALSE)
  fa$md5 = substring(rownames(fa), 1, 32)
  rownames(fa) <- NULL

  #' join our id table
  final <- fa_meta %>% select(md5, serogroup) %>% left_join(fa)

  final
}

library(seqinr)
library(dplyr)

read_fasta <- function() {
  fa15 = read.fasta("data/gnd2/solexaQA15thou_nucleotideGE10.fa")
  fa15meta = read.table("data/gnd2/solexaQA15thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  #' Hmm, some of the names in the fasta file are O's...
  o_rows <- which(substring(names(fa15), 1, 1) == "O")
  
  # now label the fa15 by serotroup rather than MD5, and merge to give id/sequence dataframe
  fa_names <- substring(names(fa15),1,32)
  
  #' Hmm, the following rows are missing stuff
  rows <- which(!(fa15meta$md5 %in% fa_names) & !fa15meta$serogroup %in% fa_names)
  
  #' reshape the fa data into a data frame
  fa15 = data.frame(do.call(rbind, fa15), stringsAsFactors = FALSE)
  fa15$md5 = substring(rownames(fa15), 1, 32)
  rownames(fa15) <- NULL
  
  #' join our id table
  fa15_m = fa15 %>% left_join(fa15meta %>% select(md5, serogroup)) %>% filter(!is.na(serogroup))
  fa15_o = fa15 %>% left_join(fa15meta %>% select(md5, serogroup), by=c('md5'='serogroup')) %>% filter(!is.na(md5.y)) %>% rename(serogroup=md5, md5=md5.y)
  fa15 = union(fa15_m, fa15_o)
  
  fa15
}

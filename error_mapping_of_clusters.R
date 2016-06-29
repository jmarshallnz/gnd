
# Preliminary stuff...

library(seqinr)
library(dplyr)
library(cluster)
library(readxl)

# The read_fasta command may need some alteration to make sure you have the .fa and .txt (metadata)
# and pull out the various names etc.

read_fasta <- function() {
  #fa15 = read.fasta("solexaQA18thou_nucleotideGE10.fa")
  #fa15meta = read.table("solexaQA18thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
  fa15 = read.fasta("15thou_403seqs.fas")
  fa15meta = read.table("solexaQA15thou_minTotalGE10.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  #' Hmm, some of the names in the fasta file are O's...
  o_rows <- which(substring(names(fa15), 1, 1) == "O")
  
  # now label the fa15 by serotroup rather than MD5, and merge to give id/sequence dataframe
  fa_names <- substring(names(fa15),1,32)
  
  #' Hmm, the following rows are missing stuff
  rows <- which(!(fa15meta$md5 %in% fa_names) & !fa15meta$serogroup %in% fa_names)
  
  #' reshape the fa data into a data frame
  fa15 = data.frame(do.call(rbind, fa15))
  fa15$md5 = substring(rownames(fa15), 1, 32)
  rownames(fa15) <- NULL
  
  #' join our id table
  fa15_m = fa15 %>% left_join(fa15meta %>% select(md5, serogroup)) %>% filter(!is.na(serogroup))
  fa15_o = fa15 %>% left_join(fa15meta %>% select(md5, serogroup), by=c('md5'='serogroup')) %>% filter(!is.na(md5.y)) %>% rename(serogroup=md5, md5=md5.y)
  fa15 = union(fa15_m, fa15_o)
  
  fa15
}

# The error mapping file is just 2 columns (id, parent). Columns are named 'sample' and 'probable_source' as that is what
# I use below (i.e. just change it if you like below)
mapping = read.csv("error_mapping.csv")

fa15 <- read_fasta()

#' compute difference maps and higlight them
map_source <- mapping %>% left_join(fa15 %>% select(-md5), by=c('sample' = 'serogroup'))
map_dest   <- mapping %>% left_join(fa15 %>% select(-md5), by=c('probable_source' = 'serogroup'))

# now check the difference between them and highlight it...
# TODO: CHANGE THIS BELOW TO DO THE MAPPING TO COLOURS BASED ON DIFFERENT DIFFERENCES (A->G etc.)
diff = map_source[,-(1:2)] != map_dest[,-(1:2)]
map_diff = cbind(map_source[,1:2], diff)

map_diff = map_diff %>% arrange(probable_source)
map_matrix = as.matrix(map_diff[,-(1:2)])
# add in an alternating thingee
map_matrix = map_matrix + 2*(as.numeric(map_diff[,2]) %% 2)
colnames(map_matrix) = 1:284
rownames(map_matrix) = map_diff[,1]
pdf(paste0("error_maps_by_base_pairs.pdf"), width=10, height=20)
par(mar=c(4,6,2,6))
image(1:ncol(map_matrix), 1:nrow(map_matrix), t(map_matrix), col=c("white", "black", "grey80", "black"), xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, 1:nrow(map_matrix), rownames(map_matrix), las=2, cex.axis=0.6)
axis(4, 1:nrow(map_matrix), map_diff[,2], las=2, cex.axis=0.6)
dev.off()

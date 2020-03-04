library(tidyverse)

source('code/read_fasta.R')

#abundance_file <- "data/patrick_costs/solexaQA1thou_minTotalGE10.txt"
#fasta_file <- "data/patrick_costs/solexaQA1thou_top50000NucleotideGE0.fa"

abundance_file <- "data/patrick_costs/pear_minTotalGE10.txt"
fasta_file <- "data/patrick_costs/pear_nucleotideGE0_manually_hacked.fa"

# read in the metadata file
abund <- read.table(abundance_file, sep="\t", header=TRUE, fill=TRUE)

md5s <- abund %>% filter(functional == "yes") %>% pull(md5)

# MORE MANUAL HACKERY
ignore <- c(9115, 9141, 9221, 9241, 9247, 9259, 9271, 9292,
            9353, 9396, 9464, 9488, 9556, 9569, 9629, 9649,
            9658, 9693, 9760, 9771, 9811, 9829, 9838, 9925, 9935, 9940)

md5s <- md5s[-ignore]
# MORE MANUAL HACKERY

fasta <- read_fasta_files(md5s, fasta_file, read_length = 484)

# yay, assemble down to what we need
fasta$sequence <- apply(fasta %>% select(-md5), 1, paste0, collapse='')

# find most abundant (> 100 reads)
minimum_parent_count <- 100
parents <- abund %>% select(md5, count = overall.counts) %>% filter(count >= minimum_parent_count) %>%
  left_join(fasta %>% select(md5, sequence))

# find closest most abundant one for all others
source("code/closest_parents.R") #< builds c++

fa.close.df <- closest_parents(fasta$md5, fasta$sequence, parents$md5, parents$sequence, min_dist = 5)

min_dist <- fa.close.df %>% group_by(id) %>% summarise(MinDist = min(Distance)) %>%
  left_join(abund %>% select(id = md5, count = overall.counts))

min_dist %>% filter(MinDist >= 10) %>% arrange(desc(count))

#save(fa.close.df, file=file.path(out_dir, out_file))

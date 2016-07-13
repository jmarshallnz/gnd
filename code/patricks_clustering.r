#' figuring out patrick's clustering
library(dplyr)

pat = read.table("renamedClusters.txt", header=TRUE)
nrow(pat)

abund = read.csv("sero_abundance.csv", row.names = 1)

abund = data.frame(serogroup=rownames(abund), abundance=rowSums(abund))

abund = abund %>% left_join(pat)

clust = abund %>% group_by(ClusterNum) %>% summarise(total=sum(abundance), max=max(abundance), wch=serogroup[which.max(abundance)])

clust %>% left_join(abund)

cent_clust <- abund %>% left_join(clust) %>% filter(mostAbundant == 'y') %>% dplyr::select(-clustValue, -seqLength)
all(cent_clust$serogroup == cent_clust$wch)

other_clust <- abund %>% left_join(clust) %>% filter(mostAbundant == 'n') %>% dplyr::select(-clustValue, -seqLength)

other_clust = other_clust %>% left_join(abund %>% filter(mostAbundant == 'y'), by=c('ClusterNum'='ClusterNum'))
all(other_clust$wch == other_clust$serogroup.y)

all_clust <- abund %>% left_join(clust) %>% dplyr::select(-clustValue, -seqLength)

# read in patricks stuff
pat_abund = read.table("Patrick_Proportions_15thou.txt", header=TRUE, row.names=1)
read_counts = colSums(read.csv("sero_abundance.csv", row.names=1))
patrick = sweep(patrick, 2, read_counts/colSums(patrick), '*')
y.patrick = rowSums(patrick)

cent_clust = all_clust %>% filter(serogroup == wch)

abund_pat = data.frame(serogroup=names(y.patrick), abundance = y.patrick)

cent_clust %>% left_join(abund_pat, by=c('serogroup' = 'serogroup'))

# seems to be some differences...
setdiff(abund_pat$serogroup, cent_clust$serogroup)
setdiff(cent_clust$serogroup, abund_pat$serogroup)

# create new abundance thingee for Patrick's stuff, using CDHIT assignment across all.

abundance = read.csv("sero_abundance.csv", row.names = 1)
abundance$serogroup = rownames(abundance)

all_clust$wch = droplevels(all_clust$wch)
new_abund = abundance %>% left_join(all_clust, by = c("serogroup" = "serogroup"))

new_abund2 = apply(new_abund[,1:96], 2, function(x) { tapply(x, new_abund$wch, sum) })
new_abund2 = data.frame(serogroup = levels(new_abund$wch), new_abund2)
write.csv(new_abund2, "patrick_cdhit_abundance.csv", row.names=FALSE)


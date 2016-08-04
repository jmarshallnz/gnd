#' figuring out patrick's clustering
library(dplyr)

source("code/read_abundance.R")

pat = read.table("temp/tableFromCluster348.txt", header=FALSE)
names(pat) <- c("serogroup", "seqLength", "ClusterNum", "repSequence")

abund = read_abundance()

abund = data.frame(serogroup=rownames(abund), abundance=rowSums(abund))

abund = abund %>% left_join(pat)

clust = abund %>% group_by(ClusterNum) %>% summarise(total=sum(abundance), max=max(abundance), wch=serogroup[which.max(abundance)])

all_clust <- abund %>% left_join(clust) %>% dplyr::select(-seqLength)

# create new abundance thingee for Patrick's stuff, using CDHIT assignment across all.

abundance = read.csv("sero_abundance.csv", row.names = 1)
abundance$serogroup = rownames(abundance)

new_abund = abundance %>% left_join(all_clust, by = c("serogroup" = "serogroup"))
new_abund$wch = factor(new_abund$wch)

new_abund2 = apply(new_abund[,1:96], 2, function(x) { tapply(x, new_abund$wch, sum) })
new_abund2 = data.frame(serogroup = levels(new_abund$wch), new_abund2)
write.csv(new_abund2, "temp/patrick_cdhit_abundance.csv", row.names=FALSE)


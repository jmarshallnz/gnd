#' compressed version of the matrix: Ignore any more than distance d away
library(dplyr)
source("code/read_abundance.R")

#' The sample. Can optionally include controls here (e.g. for MDS plot in the paper)
#appendix <- ""
#y = rowSums(read_abundance())
fa15d <- read_abundance()
fa15.dist = as.matrix(read.csv("temp/sero_dist15.csv", row.names=1))
#d = d[names(y), names(y)]



d <- 2
#dist.comp <- vector("list", nrow(fa15.dist))
dist.comp <- apply(fa15.dist, 1, function(x, d) { lapply(1:d, function(y) { which(x == y) }) }, d=2)

#' Find all pairs that are distance 1 from each other
wch_1 <- which(fa15.dist == 1)
col <- (wch_1-1) %/% nrow(fa15.dist) + 1
row <- (wch_1-1) %% nrow(fa15.dist) + 1

#' filter out the duplicates
dedupe <- unique(cbind(pmin(col,row), pmax(col,row)))
row <- dedupe[,1]
col <- dedupe[,2]

#' now pull out the abundances
abund <- rowSums(fa15d)
d <- data.frame(node1 = rownames(fa15.dist)[row],
                node2 = rownames(fa15.dist)[col],
                abund1 = abund[row],
                abund2 = abund[col], stringsAsFactors = FALSE)
d <- d %>% mutate(total = abund1 + abund2, ratio = pmin(abund1, abund2)/total)

# see how many are potentially errors, and swap n1,n2 based on abundance
dd <- d %>% filter(ratio < 0.05) %>%
  mutate(one_bigger = pmax(abund1, abund2) == abund1,
         n1 = ifelse(one_bigger, node1, node2),
         n2 = ifelse(one_bigger, node2, node1))

# find the fakers
fakers <- dd %>% select(n2,n1) %>% unique %>% arrange(n2)

# some of these are duplicates
duped <- fakers$n2[duplicated(fakers$n2)]

# link with their abundance
fa15d$st <- rownames(fa15d)
abund <- fa15d %>% left_join(fakers, by=c('st'='n2'))


# Distribute abundances across the duplicated rows
duped <- abund$st[duplicated(abund$st)]
for (i in seq_along(duped)) {
  r <- which(abund$st == duped[i])
  abund[r,1:92] <- abund[r,1:92] / length(r)
}

# Fill in n1 for those that are empty
abund <- abund %>% mutate(n1 = ifelse(is.na(n1), st, n1))

# Group by n1 and sum everything up
final <- abund %>% select(-st) %>% group_by(n1) %>% summarise_all(sum)

final$Totals <- rowSums(final[,-1])

final <- final %>% rename(ST = n1)

write.csv(final, "unfaked.csv", row.names=FALSE)

summary(final$Totals)

fakers[fakers$n2 %in% duped,]

# check none of the fakers are in the parents
which(fakers$n2 %in% parent$n1)

# add in the non-fakers


nrow(fakers)  
nrow(parent)

pdf("GND ratio of abundances.pdf", width=8, height=11)
par(mfrow=c(3,1))
hist(log10(d$ratio), xaxt = 'n', main = "Ratios from GNDs 1 SNP apart", xlab="Ratio")
axis(1, at=log10(10^c(-4:0)), labels = formatC(10^c(-4:0), width=5))

wch_2 <- which(fa15.dist == 2)
col <- (wch_2-1) %/% nrow(fa15.dist) + 1
row <- (wch_2-1) %% nrow(fa15.dist) + 1

#' now pull out the abundances
abund <- rowSums(fa15d)
d <- data.frame(node1 = rownames(fa15.dist)[row],
                node2 = rownames(fa15.dist)[col],
                abund1 = abund[row],
                abund2 = abund[col])
d <- d %>% mutate(total = abund1 + abund2, ratio = pmin(abund1, abund2)/total)

hist(log10(d$ratio), xaxt = 'n', main = "Ratios from GNDs 2 SNP apart", xlab="Ratio")
axis(1, at=log10(10^c(-4:0)), labels = formatC(10^c(-4:0), width=5))

hist(log10(d$ratio), xaxt = 'n', main = "Ratios from GNDs 2 SNP apart (zoomed in)", xlab="Ratio", ylim=c(0,350))
axis(1, at=log10(10^c(-4:0)), labels = formatC(10^c(-4:0), width=5))
dev.off()


# check which control gSTs are potential fakers, given the known ones
ctrl_gst <- rownames(fa15d)[which(fa15d[,91] > 0)]
know_gst <- rownames(fa15d)[which(fa15d[,91] > 110)]

error_gst <- setdiff(ctrl_gst, know_gst)

# find distance between error_gst and know_gst
diff <- apply(fa15.dist[error_gst, know_gst], 1, min)
where <- know_gst[apply(fa15.dist[error_gst, know_gst], 1, which.min)]
d1 <- data.frame(faker=error_gst, count=fa15d[error_gst,91], distance=diff, from=where,control="Ctrl01R1_S91")

ctrl_gst <- rownames(fa15d)[which(fa15d[,92] > 0)]
know_gst <- rownames(fa15d)[which(fa15d[,92] > 110)]
error_gst <- setdiff(ctrl_gst, know_gst)

# find distance between error_gst and know_gst
diff <- apply(fa15.dist[error_gst, know_gst], 1, min)
where <- know_gst[apply(fa15.dist[error_gst, know_gst], 1, which.min)]
d2 <- data.frame(faker=error_gst, count=fa15d[error_gst,92], distance=diff, from=where,control="Ctrl02R1_S92")

d <- rbind(d1, d2)
write.csv(d, "control_errors.csv", row.names=FALSE)

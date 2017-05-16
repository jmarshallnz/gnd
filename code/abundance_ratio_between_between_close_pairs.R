#' compressed version of the matrix: Ignore any more than distance d away

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
dd <- d %>% filter(ratio < 0.01) %>%
  mutate(one_bigger = pmax(abund1, abund2) == abund1,
         n1 = ifelse(one_bigger, node1, node2),
         n2 = ifelse(one_bigger, node2, node1))

fakers <- dd %>% select(n2) %>% unique
parent <- dd %>% select(n1) %>% unique

# check none of the fakers are in the parents
which(fakers$n2 %in% parent$n1)

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

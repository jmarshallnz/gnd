# script for comparing original data between methods
source("code/read_abundance.R")

y.all = read_abundance()
y.jm = rowSums(read.csv("no_error_abundance.csv", row.names=1)[,names(y.all)])
y.pb = rowSums(read.csv("temp/patrick_cdhit_abundance.csv", row.names=1)[,names(y.all)])
y.all = rowSums(y.all)

d = as.matrix(round(read.csv("sero_dist15.csv", row.names=1) * 284))
d = d[names(y.all), names(y.all)]
diag(d) <- Inf

# right, now create some summaries of the differences between them

# differences between the two

# These are what Patrick is missing:
extras_jm = setdiff(names(y.jm), names(y.pb))

# These are what Patrick has that jm is missing
extras_pb = setdiff(names(y.pb), names(y.jm))


y.all[extras_pb]

jm_closest_other = names(y.jm)[apply(d[extras_jm,names(y.jm)], 1, which.min)]
cbind(y.all[extras_jm], jm_closest_other, y.all[jm_closest_other])

# for each observation, find the closest other serogroup (among highest if that makes sense)
closest_abundant <- function(x, abund) {
  possibles = which(x == min(x))
  max_abund = which.max(abund[possibles])
  possibles[max_abund]
}

close_abund = colnames(d)[apply(d, 1, closest_abundant, y.all)]

mapping = data.frame(serogroup = names(y.all),
                     abundance = y.all,
                     closest = close_abund,
                     close_abund = y.all[close_abund],
                     close_dist = d[cbind(names(y.all), close_abund)], jm = names(y.all) %in% names(y.jm), pb = names(y.all) %in% names(y.pb))
rownames(mapping) = NULL
write.csv(mapping, "temp/inconsistent_methods.csv", row.names=FALSE)

#' now compare abundances for the samples.

y.all = read_abundance()
y.jm = read.csv("no_error_abundance.csv", row.names=1)[,names(y.all)]
y.pb = read.csv("patrick_cdhit_abundance.csv", row.names=1)[,names(y.all)]

# make y.jm and y.pb into the joint size we need
serogroups = union(rownames(y.pb), rownames(y.jm))

# now we need length(serogroups) colours...

prop.jm = matrix(0, length(serogroups), ncol(y.all))
rownames(prop.jm) = serogroups
colnames(prop.jm) = colnames(y.jm)
prop.jm[rownames(y.jm),] = sweep(as.matrix(y.jm), 2, colSums(y.jm), "/")
barplot(prop.jm)

prop.pb = matrix(0, length(serogroups), ncol(y.all))
rownames(prop.pb) = serogroups
colnames(prop.pb) = colnames(y.all)
prop.pb[rownames(y.pb),colnames(y.pb)] = sweep(as.matrix(y.pb), 2, colSums(y.pb), "/")
barplot(prop.pb)

# now figure out all the animal names (yay, this again...)

#' Read in metadata
source("code/read_metadata.R")

labels <- read_metadata(animal_cols = colnames(y.all))$label

# ok, now plot them somehow...
colnames(prop.pb) = labels
barplot(prop.pb, las=2)

colnames(prop.jm) = labels
o = order(colnames(prop.jm))
barplot(prop.jm[,o], las=2)
barplot(prop.pb[,o], las=2)

# choose some colours
library(RColorBrewer)
cols = brewer.pal(8, "Set1")

# now extend this to way more, by picking 32 shades of each of these 8...
# 1. Pick a base colour
# 2. Calculate the hue
# 3. Use hue = baseHue + ((240/pieces)) * piece % 240

cols = hsv(seq(0,0.8,length.out=length(serogroups)), s = 0.7, v = 0.9)

set.seed(4)
cols = sample(cols)

png("figures/barplot_jm.png", width=800, height=500)
par(mar=c(5,4,2,2))
barplot(prop.jm[,o], las=2, col=cols, border=NA, xaxs="i", cex.names = 0.8)
dev.off()

png("figures/barplot_pb.png", width=800, height=500)
par(mar=c(5,4,2,2))
barplot(prop.pb[,o], las=2, col=cols, border=NA, xaxs="i", cex.names = 0.8)
dev.off()

library(animation)
saveGIF({
  for (i in 1:30) {
    cols = sample(cols)
    par(mar=c(5,4,2,2))
    barplot(prop.pb[,o], las=2, col=cols, border=NA, xaxs="i", cex.names = 0.8)
  }}, movie.name = "test.gif", interval=runif(30, 0.01, 1), nmax=30, ani.width=600, ani.height=400)

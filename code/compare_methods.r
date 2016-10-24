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

source("code/read_abundance.R")

y.all = read_abundance()
y.jm = read.csv("no_error_abundance.csv", row.names=1)[,names(y.all)]
y.pb = read.csv("temp/patrick_cdhit_abundance.csv", row.names=1)[,names(y.all)]

# make y.jm and y.pb into the joint size we need
serogroups = union(rownames(y.pb), rownames(y.jm))

# work out proportions
prop.jm = matrix(0, length(serogroups), ncol(y.all))
rownames(prop.jm) = serogroups
colnames(prop.jm) = colnames(y.jm)
prop.jm[rownames(y.jm),] = sweep(as.matrix(y.jm), 2, colSums(y.jm), "/")

prop.pb = matrix(0, length(serogroups), ncol(y.all))
rownames(prop.pb) = serogroups
colnames(prop.pb) = colnames(y.all)
prop.pb[rownames(y.pb),colnames(y.pb)] = sweep(as.matrix(y.pb), 2, colSums(y.pb), "/")

# figure out which ones we want to keep
threshold = 0.01 # only those that make up at least 5% of one of the samples
wch.jm = names(which(apply(prop.jm, 1, max) >= threshold))
wch.pb = names(which(apply(prop.pb, 1, max) >= threshold))
wch = union(wch.jm, wch.pb)

prop.jm = prop.jm[wch,]
prop.pb = prop.pb[wch,]
prop.jm = rbind(prop.jm, Other = 1-colSums(prop.jm))
prop.pb = rbind(prop.pb, Other = 1-colSums(prop.pb))

serogroups = wch
# now figure out all the animal names (yay, this again...)

#' Read in metadata
source("code/read_metadata.R")

meta <- read_metadata(animal_cols = colnames(y.all))
colnames(prop.pb) = meta$label
colnames(prop.jm) = meta$label
o = order(as.numeric(meta$animal), meta$source)

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

png("figures/barcode_jm_pb.png", width=1000, height=800)
layout(matrix(1:3,ncol=3), width = c(10,10,1.5),height = c(1,1,0.3))
par(oma=c(6,5,1,1))
par(mar=c(0,0,0,0.5))
# try a heatmap instead
trans <- function(x) {
  x
}
image(x = seq_len(nrow(prop.jm)), y = seq_len(ncol(prop.jm)), trans(prop.jm[,o]), col=grey(seq(1, 0, length.out=21)), axes=FALSE, xlab="", ylab="")
axis(1, at=seq_len(nrow(prop.jm)), labels=row.names(prop.jm), las=2, cex.axis=0.6)
axis(2, at=seq_len(ncol(prop.jm)), labels=colnames(prop.jm)[o], las=2, cex.axis=0.8)
par(mar=c(0,0.5,0,0))
image(x = seq_len(nrow(prop.pb)), y = seq_len(ncol(prop.pb)), trans(prop.pb[,o]), col=grey(seq(1, 0, length.out=21)), axes=FALSE, xlab="", ylab="")
axis(1, at=seq_len(nrow(prop.pb)), labels=row.names(prop.pb), las=2, cex.axis=0.6)
#axis(2, at=seq_len(ncol(prop.pb)), labels=colnames(prop.pb)[o], las=2)

legend_image <- as.raster(matrix(grey(seq(0, 1, length.out=21)), ncol=1))

par(mar=c(0,0.5,2,0.5))
plot(c(0,2),c(-5,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Proportion')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0.5, 0, 1, 1)
dev.off()

# Now work out the difference between the two
diff <- prop.jm - prop.pb

# grab a palette
pal <- rev(brewer.pal(9, "BrBG"))

# work out the breaks
up_bks = quantile(diff[diff > 0], seq(1,9,by=2)/9)
lo_bks = quantile(diff[diff < 0], seq(8,0,by=-2)/9)
bks = round(pmax(up_bks, abs(lo_bks)), 5)
bks = seq(0.1, 0.9, length.out=4)

png("figures/barcode_jm_v_pb.png", width=500, height=800)
par(mfrow=c(1,1), mar=c(6,6,2,2))
vals <- cut(diff, breaks = c(-bks, bks))
prop.df <- matrix(as.numeric(vals), nrow(diff))
image(x = seq_len(nrow(prop.jm)), y = seq_len(ncol(prop.jm)), trans(prop.df[,o]), col=pal, axes=FALSE, xlab="", ylab="")
axis(1, at=seq_len(nrow(prop.jm)), labels=row.names(prop.jm), las=2, cex.axis=0.6)
axis(2, at=seq_len(ncol(prop.jm)), labels=colnames(prop.jm)[o], las=2, cex.axis=0.8)
legend('topleft', legend=levels(vals), fill = pal, cex=0.6)
dev.off()

library(animation)
saveGIF({
  for (i in 1:30) {
    cols = sample(cols)
    par(mar=c(5,4,2,2))
    barplot(prop.pb[,o], las=2, col=cols, border=NA, xaxs="i", cex.names = 0.8)
  }}, movie.name = "test.gif", interval=runif(30, 0.01, 1), nmax=30, ani.width=600, ani.height=400)

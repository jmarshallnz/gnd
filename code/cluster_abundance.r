library(readxl)
library(dplyr)

abund_samp = read.csv("no_error_abundance_with_ctrl.csv", row.names=1)
abund_samp = t(abund_samp)
abund_per_sample = abund_samp / rowSums(abund_samp)
abund_dist = as.matrix(dist(abund_per_sample))

abund.hc1 = hclust(abund_dist, method='complete')
plot(abund.hc1)

animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", rownames(abund_per_sample))
rep    <- sub("X([0-9]+)_.*", "\\1", rownames(abund_per_sample))
rownames(abund_dist) = colnames(abund_dist) = paste0(animal, "_", rep)

#' Read in metadata
meta = read_excel("gnd_seqs_metadata.xlsx", sheet=2)

samp.meta = data.frame(sample=substring(rownames(abund_per_sample),2), stringsAsFactors = FALSE)

samp.meta = samp.meta %>% left_join(meta, by=c('sample'='Library Name')) %>% dplyr::select(sample, treatment=Treatment, source=`Description [Optional] `)
samp.meta$treatment = factor(samp.meta$treatment)
levels(samp.meta$treatment) = c("bf", "ct")
samp.meta$source = factor(samp.meta$source)
levels(samp.meta$source) = c("fec", "pob", "por", "pre")
samp.meta = samp.meta %>% mutate(animal = sub(".*_([0-9ctrl]+)_.*", "\\1", sample),
                                 label = paste0(animal, "_", treatment, "_", source))

# relabel the ctrls
wch <- samp.meta$animal == "ctrl"
ctrl_labs = sub("([0-9ctrl]+)_.*", "\\1", samp.meta$sample)
samp.meta$label[wch] <- paste0("ctrl_", ctrl_labs[wch])

# add labels
rownames(abund_dist) = colnames(abund_dist) = samp.meta$label

# do some clustering
samp.hc = hclust(as.dist(abund_dist))
pdf("sample_tree1.pdf", width=12, height=8)
plot(samp.hc, main="Sample clustering (method 1)", xlab="", sub="", cex=0.7)
dev.off()

samp.hc2 = hclust(as.dist(abund_dist), method="average")
pdf("sample_tree2.pdf", width=12, height=8)
plot(samp.hc2, main="Sample clustering (method 2)", xlab="", sub="", cex=0.7)
dev.off()

#' MDS below here by the looks
library(MASS)
set.seed(5)
mds <- monoMDS(abund_dist,k=2)

x <- mds$points[,1]
y <- mds$points[,2]
#z <- mds$points[,3]

library(RColorBrewer)
#plot(x, y, main="MDS plot for MLST", type="n")
#text(x,y,labels=mlst2$ST, col=as.numeric(mlst2$Group))
#pdf("mds.pdf", width=5, height=4)
#ar(pin=c(5,4), omi=rep(0,4))
pdf("mds.pdf", width=10, height=8)
cols = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
plot(x, y, main="", type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(min(x), max(x)+0.5))
points(x,y,col=cols[as.numeric(as.factor(samp.meta$animal))], pch = 16+as.numeric(as.factor(samp.meta$source)))
legend("topright", legend=levels(as.factor(samp.meta$animal)), pch=19, col=cols, cex=0.8)
legend("topleft", legend=levels(as.factor(samp.meta$source)), pch=16+1:4, col=cols[24], cex=0.8)
#legend(2, 2, legend=c("chicken", "cow", "sheep", "water"), fill=1:4)

wch = samp.meta$animal == "ctrl"
points(x[wch], y[wch], cex=2)

dev.off()

wch = wch | rowSums(abund_samp) < 200
perm_diss_red = abund_dist[!wch, !wch]
meta_red = samp.meta[!wch,]

mds.red <- monoMDS(perm_diss_red,k=2)

x <- mds.red$points[,1]
y <- mds.red$points[,2]

#pdf("mds.pdf", width=5, height=4)
#ar(pin=c(5,4), omi=rep(0,4))
#pdf("mds.pdf", width=10, height=8)
cols = c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))
plot(x, y, main="", type="n", xlab="", ylab="", xaxt="n", yaxt="n")
points(x,y,col=cols[as.factor(meta_red$animal)], pch = 16+as.numeric(as.factor(meta_red$source)))
legend("topright", legend=levels(as.factor(meta_red$animal)), pch=19, col=cols, cex=0.8)
legend("topleft", legend=levels(as.factor(meta_red$source)), pch=16+1:4, col=cols[24], cex=0.8)

# TODO: BARPLOT FOR PATRICKS AND MINE, WITH HOPEFULLY NOT HORRID COLOURS!
barplot(t(abund_samp/rowSums(abund_samp)))

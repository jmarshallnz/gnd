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

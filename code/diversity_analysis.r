# some diversity stuff
library(vegan)
source("code/read_abundance.R")

abund_samp = read.csv("no_error_abundance.csv", row.names=1)
abund_samp = t(abund_samp)
abund_all = read_abundance()
abund_all = t(abund_all)

# compute diversity based on no-error data and on all data
div.simpson = vegan::diversity(abund_samp, "simpson")
div.shannon = vegan::diversity(abund_samp, "shannon")
div.simpson.all = vegan::diversity(abund_all, "simpson")
div.shannon.all = vegan::diversity(abund_all, "shannon")

pdf("figures/simpson.pdf", width=6, height=6)
plot(div.simpson, div.simpson.all, xlab="Simpson diversity (corrected)", ylab="Simpson diversity (all)")
abline(0,1)
dev.off()

pdf("figures/shannon.pdf", width=6, height=6)
plot(div.shannon, div.shannon.all, xlab="Shannon diversity (corrected)", ylab="Shannon diversity (all)")
abline(0,1)
dev.off()

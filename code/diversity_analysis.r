# some diversity stuff
library(vegan)

source("code/read_abundance.R")
source("code/read_metadata.R")

abund_samp = read.csv("no_error_abundance.csv", row.names=1)
abund_samp = t(abund_samp)
abund_all = read_abundance()
abund_all = t(abund_all)

# compute diversity based on no-error data and on all data
div.simpson = vegan::diversity(abund_samp, "simpson")
div.shannon = vegan::diversity(abund_samp, "shannon")
div.simpson.all = vegan::diversity(abund_all, "simpson")
div.shannon.all = vegan::diversity(abund_all, "shannon")

meta <- read_metadata(animal_cols = rownames(abund_samp))

# table for things
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 22, 21), size=c(2, 2, 1, 1))
cols <- read.csv("data/colours.csv", stringsAsFactors = FALSE)

meta = meta %>% mutate(animal = as.numeric(animal)) %>% left_join(cols, by=c('animal' = 'animal_tag')) %>% left_join(shapes, by=c('source' = 'treatment'))

pdf("figures/simpson.pdf", width=6, height=6)
plot(div.simpson, div.simpson.all, xlab="Simpson diversity (corrected)", ylab="Simpson diversity (all)", col="black", pch=meta$shape, bg=meta$hex, cex=meta$size)
abline(0,1)
dev.off()

pdf("figures/shannon.pdf", width=6, height=6)
plot(div.shannon, div.shannon.all, xlab="Shannon diversity (corrected)", ylab="Shannon diversity (all)", col="black", pch=meta$shape, bg=meta$hex, cex=meta$size)
abline(0,1)
dev.off()

pdf("figures/shannon_facet.pdf", width=6, height=7)
par(mfrow=c(5,4), mar=rep(0.5, 4), oma=c(4,4,1,1))
meta$shannon = div.shannon
meta$shannon.all = div.shannon.all
animals <- unique(meta$animal)
for (i in seq_along(animals)) {
  an_vals <- meta %>% filter(animal == animals[i])
  plot(NULL, xlim=c(0,3), ylim=c(0,3), axes = FALSE)
  abline(0,1)
  points(an_vals$shannon, an_vals$shannon.all, col="black", pch=an_vals$shape, bg=an_vals$hex, cex=an_vals$size)
  if (i %% 4 == 1) {
    axis(2, at=0:3, cex.axis=0.8, las=1, col.axis = "grey50", col = "grey50", tck = -.02, labels=FALSE)
    mtext(0:3, side = 2, at=0:3, line=0.5, las=1, cex = 0.5, col = "grey50")
  }
  if (i >= 16) {
    axis(1, at=0:3, cex.axis=0.8, las=1, col.axis = "grey50", col = "grey50", tck = -.02, labels=FALSE)
    mtext(0:3, side = 1, at=0:3, line=0.5, las=1, cex = 0.5, col = "grey50")
  }
  box(col = "grey50")
}
mtext("Shannon diversity (reduced data)", side = 1, outer=TRUE, at = 0.5, line=2, col = "grey50")
mtext("Shannon diversity (all data)", side = 2, outer=TRUE, at = 0.5, line=2, col = "grey50")
dev.off()

# some diversity stuff
library(vegan)

source("code/read_abundance.R")
source("code/read_metadata.R")

abund_samp = read.csv("no_error_abundance.csv", row.names=1)
abund_samp = t(abund_samp)
abund_all = read_abundance()
abund_all = t(abund_all)

# bootstrap method for standard error of diversities
boot_diversity <- function(abund, samples = 100) 
{
  npop  <- nrow(abund) # number of populations
  nspec <- ncol(abund) # number of species
  tspec <- apply(abund, 1, sum) # total number of species

  div  <- vegan::diversity(abund)
  boot <- matrix(NA, npop, samples)
  rownames(boot) <- rownames(abund)
  for (j in 1:samples) {
    # generate bootstrap sample of our population
    boot_abund <- t(apply(abund, 1, function(x) { rmultinom(n = 1, size = sum(x), prob = x) }))
    boot[,j] <- vegan::diversity(boot_abund)
  }
  return(list(diversity = div, bootstrap = boot))
}

# compute shannon diversity index
div.samp <- boot_diversity(abund_samp, 100)
div.samp$lci <- apply(div.samp$bootstrap, 1, quantile, 0.025)
div.samp$uci <- apply(div.samp$bootstrap, 1, quantile, 0.975)

meta <- read_metadata(animal_cols = rownames(abund_samp))

# table for things
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 22, 21), size=c(2, 2, 1, 1))
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 23, 24), size=c(1, 1, 1, 1))

cols <- read.csv("data/colours.csv", stringsAsFactors = FALSE)

meta = meta %>% mutate(animal = as.numeric(animal)) %>% left_join(cols, by=c('animal' = 'animal_tag')) %>% left_join(shapes, by=c('source' = 'treatment'))

# totals...
totals <- apply(abund_samp, 2, function(x) { tapply(x, meta$source, sum) })
totals.div <- boot_diversity(totals)
totals.div$lci <- apply(totals.div$bootstrap, 1, quantile, 0.025)
totals.div$uci <- apply(totals.div$bootstrap, 1, quantile, 0.975)

# faceted version
pdf("figures/shannon_facet.pdf", width=6, height=7)
par(mfrow=c(5,4), mar=rep(0.5, 4), oma=c(4,4,1,1))
meta$shannon = div.samp$diversity
meta$shannon.lc = div.samp$lci
meta$shannon.uc = div.samp$uci
#meta$shannon.all = div.all$bootstrap
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

pdf("figures/mds_and_shannon.pdf", width=8, height=10)
par(mfrow=c(2,1), mar=c(4,4,2,1))

abund_samp = read.csv("no_error_abundance_with_ctrl.csv", row.names=1)
abund_samp = t(abund_samp)
abund_per_sample = abund_samp / rowSums(abund_samp)
abund_dist = as.matrix(dist(abund_per_sample))
meta <- read_metadata(animal_cols = rownames(abund_samp))

# table for things
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 22, 21), size=c(2, 2, 1, 1))
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 23, 24), size=c(1, 1, 1, 1))
cols <- read.csv("data/colours.csv", stringsAsFactors = FALSE) %>% arrange(animal_tag)

meta = meta %>% mutate(animal = as.numeric(animal)) %>% left_join(cols, by=c('animal' = 'animal_tag')) %>% left_join(shapes, by=c('source' = 'treatment'))

set.seed(5)


mds <- monoMDS(abund_dist,k=2)

x <- mds$points[,1]
y <- mds$points[,2]

plot(x, y, main="", type="n", xaxt='n', xlab="Coordinate 1", ylab="Coordinate 2", xlim=c(min(x), max(x)+0.5), bty='n')
points(x,y,col="black", pch=meta$shape, bg=meta$hex, cex=meta$size)
axis(1, at=seq(-2.5,2.5,by=0.5))
legend("topright", legend=c(cols$animal_tag, "ctrl"), pch=c(rep(19, 20), 1), col=c(cols$hex, "black"), cex=0.8, pt.cex=c(rep(0.8, 20), 2), bty="n")
legend("topleft", legend=shapes$treatment, pch=shapes$shape, col="grey30", cex=0.8, bty="n")
#legend(2, 2, legend=c("chicken", "cow", "sheep", "water"), fill=1:4)

wch = is.na(meta$animal)
points(x[wch], y[wch], cex=2, col="grey30")

#par(xpd = NA)
mtext('a.', side=2, at=2, cex=2, font = 2, las=1)

#dev.off()

axis_col <- 'black'
#pdf("figures/shannon_animal.pdf", width=9, height=5)

#par(mfrow=c(1,1), mar=c(4,4,2,2))


abund_samp = read.csv("no_error_abundance.csv", row.names=1)
abund_samp = t(abund_samp)
abund_all = read_abundance()
abund_all = t(abund_all)

# bootstrap method for standard error of diversities
boot_diversity <- function(abund, samples = 100) 
{
  npop  <- nrow(abund) # number of populations
  nspec <- ncol(abund) # number of species
  tspec <- apply(abund, 1, sum) # total number of species
  
  div  <- vegan::diversity(abund)
  boot <- matrix(NA, npop, samples)
  rownames(boot) <- rownames(abund)
  for (j in 1:samples) {
    # generate bootstrap sample of our population
    boot_abund <- t(apply(abund, 1, function(x) { rmultinom(n = 1, size = sum(x), prob = x) }))
    boot[,j] <- vegan::diversity(boot_abund)
  }
  return(list(diversity = div, bootstrap = boot))
}

# compute shannon diversity index
div.samp <- boot_diversity(abund_samp, 100)
div.samp$lci <- apply(div.samp$bootstrap, 1, quantile, 0.025)
div.samp$uci <- apply(div.samp$bootstrap, 1, quantile, 0.975)

meta <- read_metadata(animal_cols = rownames(abund_samp))

# table for things
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 22, 21), size=c(2, 2, 1, 1))
shapes <- data.frame(treatment=c("pob", "por", "fec", "pre"), shape=c(22, 21, 23, 24), size=c(1, 1, 1, 1))

cols <- read.csv("data/colours.csv", stringsAsFactors = FALSE)

meta = meta %>% mutate(animal = as.numeric(animal)) %>% left_join(cols, by=c('animal' = 'animal_tag')) %>% left_join(shapes, by=c('source' = 'treatment'))

# totals...
totals <- apply(abund_samp, 2, function(x) { tapply(x, meta$source, sum) })
totals.div <- boot_diversity(totals)
totals.div$lci <- apply(totals.div$bootstrap, 1, quantile, 0.025)
totals.div$uci <- apply(totals.div$bootstrap, 1, quantile, 0.975)

# faceted version
meta$shannon = div.samp$diversity
meta$shannon.lc = div.samp$lci
meta$shannon.uc = div.samp$uci

# animal version
meta = meta %>% mutate(xvar = as.numeric(as.factor(animal)) + as.numeric(as.factor(source))/4)

par(mar=c(4,4,2,0))
plot(NULL, ylim=c(0,3.2), xlim=c(1.25,23), axes = FALSE, xlab="", ylab="")
segments(meta$xvar, meta$shannon.lc, meta$xvar, meta$shannon.uc)
points(meta$xvar, meta$shannon, xlab="Animal", ylab="Shannon diversity", pch = meta$shape, bg= meta$hex, cex=meta$size)

# totals
totals.div = data.frame(treatment = names(totals.div$diversity), diversity = totals.div$diversity, lc = totals.div$lci, uc = totals.div$uci,
                        xvar = seq(22, 22.75, by=0.25))
totals.div = totals.div %>% left_join(shapes)
segments(totals.div$xvar, y0 = totals.div$lc, y1 = totals.div$uc)
points(totals.div$xvar, totals.div$diversity, xlab="Animal", ylab="Shannon diversity", pch = totals.div$shape, bg= "grey50", cex=totals.div$size)
axis(2, at=0:3, cex.axis=0.8, las=1, col.axis = axis_col, col = axis_col, tck = -.02, labels=FALSE)
mtext(0:3, side = 2, at=0:3, line=0.5, las=1, cex = 1, col = axis_col)
axis(1, at=c(1:20,22)+0.5, cex.axis=0.8, las=1, col.axis = axis_col, col = axis_col, tck = -.02, labels=FALSE)
mtext(levels(as.factor(meta$animal)), side = 1, at=1:20+0.5, line=0.5, las=1, cex = 1, col = axis_col)
mtext("Total", side = 1, at=22.5, line=0.5, las=1, cex = 1, col = axis_col)
abline(v=21.5, col=axis_col)

mtext("Shannon diversity", side = 2, at = 1.5, line=2, col = axis_col)
mtext("Animal", side = 1, at = 11, line=2.5, col = axis_col)

legend("bottomright", legend = totals.div$treatment, pch = totals.div$shape, pt.cex = totals.div$size, pt.bg = "grey50", bty = "n")

par(xpd = NA)
mtext('b.', side=2, at=3.5, cex=2, las=1, font=2)

dev.off()

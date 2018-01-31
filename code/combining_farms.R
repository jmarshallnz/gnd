library(dplyr)
library(tidyr)
library(readxl)

source("code/read_abundance.R")
source("code/read_fasta.R")

read_farm_data <- function(farm, qa, minTotal) {
  # Read in abundances for farms
  base_dir <- paste0("temp/farm", farm)
  fasta <- read_fasta(farm=farm, qa=qa, minTotal=minTotal) %>%
    select(md5, gST=serogroup)
  abund <- read_abundance(file=file.path(base_dir, paste0("sero_abundance", qa, "_", minTotal, ".csv"))) %>%
    tibble::rownames_to_column("gST") %>%
    left_join(fasta) %>%
    mutate(gST = ifelse(substring(gST, 1, 1) == "O", gST, md5)) %>%
    select(-md5) %>%
    gather(Library, Count, starts_with('X'), starts_with('Ctrl')) %>%
    extract(Library, into=c('Library','Junk'), 'X([0-9]+AGR)_(.*)') %>%
    filter(substring(Junk, 1, 1) != "B") %>%
    select(-Junk)
}

f1 <- read_farm_data(1, 3, 10)
f2 <- read_farm_data(2, 3, 10)

# read in the metadata
read_farm_meta <- function(farm) {
  read_excel("data/gnd2/meta/NZGL02276 metadata_021117.xlsx", sheet=farm) %>%
    rename(Library = `Library no.`, AnimalID=`Animal ID`, Age=`Age (d)`)
}
meta <- bind_rows(read_farm_meta(1), read_farm_meta(2))

# join everything up, removing animal id 69 and any library that is empty
f <- rbind(f1, f2) %>%
  spread(Library, Count, fill=0) %>%
  gather(Library, Count, -gST) %>%
  left_join(meta) %>%
  filter(AnimalID != 69, !is.na(Library)) %>%
  group_by(Library) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  ungroup

# Find the top 20 so we can do some consistent looking plots
top20 <- f %>% group_by(gST) %>% summarize(TotalProp = sum(Proportion)) %>%
  top_n(20, TotalProp) %>% pull(gST)

# Filter these out and assign an 'Other' group
extras <- f %>% filter(gST %in% top20) %>% summarize(Proportion=1-sum(Proportion)) %>%
  mutate(gST = "Other")

# embiggen everything so that all animals share these (otherwise they won't across
# farms)

top20_data <- f %>% filter(gST %in% top20)

library(ggplot2)
library(ggjoy)
library(gridExtra)
pdf("farms_top20.pdf", width=24, height=12)
g1 = ggplot(top20_data %>% filter(Farm == "GB"),
       aes(Age, y=factor(AnimalID), height=Proportion, fill=AnimalID)) +
  geom_ridgeline(size=0.1) +
  facet_wrap(~gST) +
  theme_bw() + ylab("Animal") + xlab("Age (days)") +
  guides(fill="none") +
  ggtitle("Farm 1")

g2 = ggplot(top20_data %>% filter(gST != "Other", Farm != "GB"),
       aes(Age, y=factor(AnimalID), height=Proportion, fill=AnimalID)) +
  geom_ridgeline(size=0.1) +
  facet_wrap(~gST) +
  theme_bw() + ylab("Animal") + xlab("Age (days)") +
  guides(fill="none") +
  ggtitle("Farm 2")
grid.arrange(g1, g2, ncol=2)
dev.off()


# Righto, do some MDS shenanigans now. First we need to filter out
# a bunch of types that have really low counts. We need to do this
# per library I guess and then find the ones that account for 99%
# of each library

topbunch <- f %>% group_by(Library) %>% arrange(desc(Proportion)) %>% mutate(CumProp=cumsum(Proportion)) %>%
  filter(CumProp < 0.90) %>% pull(gST) %>% unique

final_dat <- f %>% filter(gST %in% topbunch)

# OK, now compute the stuff we need for an MDS plot
final_dist <- final_dat %>%
  select(gST, Library, Proportion) %>%
  spread(Library, Proportion) %>% tibble::remove_rownames() %>% as.data.frame %>% tibble::column_to_rownames('gST') %>%
  t %>% dist %>% as.matrix

# do an MDS plot
library(vegan)
mds <- monoMDS(final_dist,k=2)
x <- mds$points[,1]
y <- mds$points[,2]

plot(x, y)
dat <- data.frame(Library=rownames(final_dist), x=x, y=y) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique)

pdf("mds_stressing_me_out.pdf", width=10, height=10)
dat <- data.frame(Library=rownames(final_dist), x=x, y=y) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique)
plot(y ~ x, col=factor(Farm), cex=sqrt(Age/50), data=dat)
dev.off()

# Try a PCA plot
final_pca <- final_dat %>%
  select(gST, Library, Proportion) %>%
  spread(Library, Proportion) %>% tibble::remove_rownames() %>% as.data.frame %>% tibble::column_to_rownames('gST') %>% t

# Try a logit transformation? Use minimum value I guess?
min_value <- min(final_pca[final_pca > 0]) / 10
final_pca[final_pca <= 0] <- min_value
f <- log(final_pca) - log(final_pca[,1])

pca <- prcomp(f[,-1])
dat <- data.frame(Library=rownames(final_pca), pca$x) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique) %>%
  mutate(Farm = factor(ifelse(Farm == "GB", 1, 2))) %>% mutate(Animal=factor(AnimalID))

png("pca_both_farms.png", width=960, height=640)
ggplot(dat, aes(x=PC1, y=PC2, col=Farm, size=Age)) + geom_point(alpha=0.5) +
  theme_bw(base_size = 20) + scale_color_hue(l=45)
dev.off()

png("pca_farm2.png", width=960, height=640)
#ggplot(dat %>% filter(Farm == 2, PC1 < 20), aes(x=PC1, y=PC2, col=Age)) + geom_point(alpha=0.5, col=scale_colour_hue(l=45)$palette(2)[2]) +
ggplot(dat %>% filter(Farm == 2, PC1 < 20), aes(x=PC1, y=PC2, col=Age)) + geom_point(alpha=0.5, size=6) +
  theme_bw(base_size = 20) # + scale_color_hue(l=45)
dev.off()

png("pca_farm2_animal.png", width=960, height=640)
#ggplot(dat %>% filter(Farm == 2, PC1 < 20), aes(x=PC1, y=PC2, col=Age)) + geom_point(alpha=0.5, col=scale_colour_hue(l=45)$palette(2)[2]) +
ggplot(dat %>% filter(Farm == 2, PC1 < 20), aes(x=PC1, y=PC2, col=Animal)) + geom_point(alpha=0.5, size=6) +
  theme_bw(base_size = 20)  + scale_color_discrete(l=45)
dev.off()

pdf("pca_both_farms_zoom2_by_animal.pdf", width=12, height=8)
ggplot(dat %>% filter(Farm == 2, PC1 < 20), aes(x=PC1, y=PC2, col=Animal, size=Age)) + geom_point(alpha=0.5) +
  theme_bw()
dev.off()

# shannon diversity...
div <- diversity(final_pca)
div_dat <- data.frame(Shannon=div, Library=names(div)) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique)

ggplot(div_dat, aes(x=Age, y=Shannon, col=factor(AnimalID))) + geom_line() + facet_wrap(~AnimalID) +
  theme_bw() + guides(col='none') + ylab("Shannon diversity")

final_div <- final_dat %>%
  select(gST, Library, Count) %>%
  spread(Library, Count) %>% tibble::remove_rownames() %>% as.data.frame %>% tibble::column_to_rownames('gST') %>% t

# hmm, the bootstrapping is doing something weird. Mostly, due to the huge numbers of reads, the SE
# are way too small to really be believable. But it is just a resampling, right?
boot_div <- boot_diversity(final_div)
boot_div$lci <- apply(boot_div$bootstrap, 1, quantile, 0.025)
boot_div$uci <- apply(boot_div$bootstrap, 1, quantile, 0.975)

div_dat <- data.frame(Shannon=boot_div$diversity, LI=boot_div$lci, UI=boot_div$uci, Library=names(boot_div$diversity)) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique)

pdf("shannon_by_animal.pdf", width=12, height=8)
ggplot(div_dat, aes(x=Age, y=Shannon, col=factor(AnimalID))) + 
  geom_line(size=1) + facet_wrap(~AnimalID) +
  theme_bw() + guides(col='none') + ylab("Shannon diversity")
dev.off()

# redo PCA only on farm 1
final_pca <- final_dat %>% filter(Farm != "GB") %>%
  select(gST, Library, Proportion) %>%
  spread(Library, Proportion) %>% tibble::remove_rownames() %>% as.data.frame %>% tibble::column_to_rownames('gST') %>% t

# Try a logit transformation? Use minimum value I guess?
min_value <- min(final_pca[final_pca > 0]) / 10
final_pca[final_pca <= 0] <- min_value
f <- log(final_pca) - log(final_pca[,1])

pca <- prcomp(f[,-1])
dat <- data.frame(Library=rownames(final_pca), pca$x) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique) %>%
  mutate(Farm = factor(ifelse(Farm == "GB", 1, 2)))

ggplot(dat, aes(x=PC1, y=PC2, col=factor(AnimalID), size=Age)) + geom_point(alpha=0.5) +
  theme_bw()

pdf("pca_on_logit_scale_seems_better.pdf", width=12, height=8)
plot(PC2 ~ PC1, data=dat, col=factor(Farm), cex=sqrt(Age/50), pch=as.numeric(factor(AnimalID)))
dev.off()

wch_odd_ones <- dat %>% filter(PC1 > 40) %>% select(-starts_with('PC')) %>% pull(Library)
barplot(f[wch_odd_ones[1],])

mds <- monoMDS(dist(f),k=2)
x <- mds$points[,1]
y <- mds$points[,2]
dat <- data.frame(Library=rownames(final_dist), x=x, y=y) %>%
  left_join(final_dat %>% select(Library, Farm, AnimalID, Age) %>% unique)
plot(y ~ x, col=factor(Farm), cex=sqrt(Age/50), data=dat)

pca_time <- final_dat %>%
  select(gST, Library, Proportion) %>%
  spread(Library, Proportion) %>% tibble::remove_rownames() %>% as.data.frame %>% tibble::column_to_rownames('gST') %>% t

# hmm, try rejigging so output is timeseries per gST. So model is
# n_{ijk} = gST[k] * Farm[i]:Animal[j] where n_{ijk} is a vector by sampling
time_outcome <- final_dat %>% select(gST, Proportion, Farm, AnimalID, Sampling) %>%
  spread(Sampling, Proportion)

pca_time <- prcomp(time_outcome %>% select(`02`:`09`))
plot(pca_time)

dat <- cbind(time_outcome, pca_time$x)
plot(PC2 ~ PC1, data=dat, col=factor(Farm))

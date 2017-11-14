library(dplyr)
library(tidyr)
library(readxl)

source("code/read_abundance.R")

farm <- 2
qa <- 3 # quality level
minTotal <- 10 # minimum abundance level

#' Read in abundances for farms
base_dir <- paste0("temp/farm", farm)
abund <- read_abundance(file=file.path(base_dir, paste0("sero_abundance", qa, "_", minTotal, ".csv"))) %>%
  tibble::rownames_to_column("gST") %>%
  gather(Library, Count, starts_with('X'), starts_with('Ctrl')) %>%
  extract(Library, into='Library', 'X([0-9]+AGR)')

#' Read in the metadata
farm <- read_excel("data/gnd2/meta/NZGL02276 metadata_021117.xlsx", sheet=farm) %>%
  rename(Library = `Library no.`, AnimalID=`Animal ID`, Age=`Age (d)`)

#' Find the largest 10 (or more if ties) from each library, excluding controls
l <- abund %>% filter(!is.na(Library)) %>% group_by(Library) %>%
  mutate(Proportion = Count/sum(Count))

#' Filter down to just the top 20 so we can do some consistent looking plots
top20 <- l %>% group_by(gST) %>% summarize(TotalProp = sum(Proportion)) %>%
  top_n(20, TotalProp) %>% pull(gST)

#' filter these out and assign an 'Other' group
extras <- l %>% filter(gST %in% top20) %>% summarize(Proportion=1-sum(Proportion)) %>%
  mutate(gST = "Other")

f <- l %>% filter(gST %in% top20) %>% bind_rows(extras) %>% left_join(farm)

# righto, let's see if we can plot this shiz
library(ggplot2)


ggplot(f %>% filter(AnimalID != 69), aes(Age, Proportion, fill=gST)) +
  geom_area() +
  facet_wrap(~AnimalID, nrow=2)

# hmm, that's kinda sucky, at least colour-wise.
ggplot(f %>% filter(AnimalID != 69), aes(Age, Proportion)) +
  geom_line() +
  facet_grid(AnimalID~gST)

pdf("farm2_top20.pdf", width=12, height=12)
ggplot(f %>% filter(gST != "Other"),
       aes(Age, y=factor(AnimalID), height=Proportion, fill=AnimalID)) +
  geom_ridgeline(size=0.1) +
  facet_wrap(~gST) +
  theme_bw() + ylab("Animal") + xlab("Age (days)") +
  guides(fill="none")
dev.off()

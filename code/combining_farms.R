library(dplyr)
library(tidyr)
library(readxl)

source("code/read_abundance.R")
source("code/read_fasta.R")

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
farm_meta <- read_excel("data/gnd2/meta/NZGL02276 metadata_021117.xlsx", sheet=farm) %>%
  rename(Library = `Library no.`, AnimalID=`Animal ID`, Age=`Age (d)`)

#' Find the largest 10 (or more if ties) from each library, excluding controls
l <- abund %>% filter(!is.na(Library)) %>% group_by(Library) %>%
  mutate(Proportion = Count/sum(Count))

f <- read_fasta(farm=farm, qa=qa, minTotal=minTotal)

l1 <- l %>% left_join(f %>% select(md5, gST=serogroup))

# YES, MD5's are different so need to use those primarily when doing this shizz

# So, simple idea is to get those that explain 99% of each type using cumsum or something?
l1 %>% ungroup %>%
  left_join(l2 %>% ungroup %>% select(gST, md5) %>% unique, by="md5") %>%
  select(gST.x, gST.y, md5) %>% filter(!is.na(gST.y))

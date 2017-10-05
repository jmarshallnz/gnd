library(dplyr)
library(tidyr)
library(cluster)

source("code/read_abundance.R")
source("code/read_fasta.R")

qa <- 3

#' Read in abundances
fa15d <- read_abundance(file=paste0("temp/sero_abundance", qa, ".csv"))

#' Read in sequences
fa15 <- read_fasta(qa = qa)

ctrl <- fa15d %>% select(Ctrl1 = Ctrl01R1_S91counts, Ctrl2 = Ctrl02R1_S92counts) %>%
  tibble::rownames_to_column("gST") %>% gather(Control, Count, Ctrl1:Ctrl2) %>%
  filter(Count > 0) %>% left_join(fa15, by=c('gST' = 'serogroup'))

ctrl1 <- ctrl %>% filter(Control == "Ctrl2")
genuine <- ctrl1 %>% filter(Count > 100) %>% pull(gST)

# compute distance matrix
DM <- ctrl1 %>% select(starts_with('X')) %>%
  as.matrix %>% as.data.frame %>% daisy %>% as.matrix
colnames(DM) <- rownames(DM) <- ctrl1$gST

dist <- DM %>% as.data.frame %>% select(one_of(genuine)) %>% tibble::rownames_to_column("gST") %>%
  gather(Parent, Distance, -gST) %>% group_by(gST) %>% summarize(Parent=Parent[which.min(Distance)], Distance = min(Distance)*284)

fakers <- ctrl1 %>% filter(Count < 300) %>% select(-Control) %>% left_join(dist) %>%
  left_join(ctrl1 %>% filter(Count > 300) %>% select(-Control, -Count, -md5), by=c("Parent"="gST")) %>%
  gather(Gene, Value, starts_with('X')) %>% extract(Gene, into=c('Gene', 'Which'), regex="X([0-9]+)\\.([xy])",convert=TRUE) %>%
  spread(Which, Value) %>% mutate(z = ifelse(x==y,"",paste0(y,x)), Gene=sprintf('%03i',Gene))

# draw a picture
fk <- fakers %>% select(-x,-y) %>% spread(Gene,z)

image(1:284,1:nrow(fk), t(fk[,-c(1:5)]!=""), col=c('white','black'))

# the errors don't appear to be what you'd expect if proportional to the parents.
# Errors from O121 and O153 are way overrepresented.
fk %>% filter(Distance == 1) %>% pull(Parent) %>% table

# Number of error gSTs and total count. It doesn't look proportional to parents really...
fk %>% filter(Distance == 1) %>% group_by(Parent) %>% summarize(Num=n(), Count=sum(Count)) %>%
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  mutate(ErrorRate = Count/ParentCount*100)

# Do the same thing, but by gST instead?
fk %>% filter(Distance == 1) %>% select(gST, Count, Parent) %>%
  left_join(ctrl1 %>% filter(Count > 300) %>% select(Parent=gST, ParentCount=Count)) %>%
  mutate(ErrorRate = Count/ParentCount*100) %>% filter(Count > 1) %>% arrange(desc(ErrorRate))

# hmm, this seems WAY too variable. What about counts across all sites?
all_counts <- fa15d %>% tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>% group_by(gST) %>% summarize(Count=sum(Count))

fk %>% filter(Distance == 1) %>% select(gST, Parent, Count) %>%
  left_join(all_counts %>% rename(Total=Count)) %>% group_by(Parent) %>%
  summarize(Num=n(), Count=sum(Count), TotalCount=sum(Total)) %>%
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  left_join(all_counts %>% rename(Parent=gST,AllParentCount=Count)) %>%
  mutate(ErrorRateCount = Count/ParentCount*100, ErrorRateNum=Num/ParentCount*100,
         AllErrorRateCount = Count/AllParentCount*100, AllErrorRateNum=Num/AllParentCount*100)

#' Have a look at this by GST
fk %>% filter(Distance == 1) %>% select(gST, Parent, Count) %>%
  left_join(all_counts %>% rename(Total=Count)) %>% 
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  left_join(all_counts %>% rename(Parent=gST,AllParentCount=Count)) %>%
  mutate(ErrorRateCount = Count/ParentCount*100,
         AllErrorRateCount = Count/AllParentCount*100) %>% filter(Count > 1) %>% arrange(desc(ErrorRateCount))

load(paste0("temp/sero_dist",qa,"_df.Rda"))
close <- fa.dist.df %>% filter(Distance == 1)

#' Hmm, need to find the distance from the ones we have here to all others, as they could have arisen from
#' a bunch of different parents
all_dist1 <- fk %>% select(gST, Parent, Count) %>% left_join(close) %>%
  left_join(all_counts %>% rename(gST2=gST, All=Count)) %>%
  group_by(gST, Parent) %>% summarize(AllDist1 = sum(All))

fk %>% filter(Distance == 1) %>% select(gST, Parent, Count) %>%
  left_join(all_counts %>% rename(Total=Count)) %>%
  left_join(all_dist1) %>%
  group_by(Parent) %>%
  summarize(Num=n(), Count=sum(Count), TotalCount=sum(Total), AllDist1Count=sum(AllDist1)) %>%
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  left_join(all_counts %>% rename(Parent=gST,AllParentCount=Count)) %>%
  mutate(ErrorRateCount = Count/ParentCount*100, ErrorRateNum=Num/ParentCount*100,
         AllErrorRateCount = Count/AllParentCount*100, AllErrorRateNum=Num/AllParentCount*100)

#' and for just the gSTs
fk %>% filter(Distance == 1) %>% select(gST, Parent, Count) %>%
  left_join(all_counts %>% rename(Total=Count)) %>%
  left_join(all_dist1) %>%
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  left_join(all_counts %>% rename(Parent=gST,AllParentCount=Count)) %>%
  mutate(ErrorRateCount = Count/ParentCount*100, AllErrorRateCount = Count/AllParentCount*100) %>%
  filter(Count > 1) %>% arrange(desc(ErrorRateCount))

#' Could it be possible that the 'fakers' we see in controls that seem at way too higher rates
#' proportionally to be errors from within the controls or even to be errors from outside
#' of controls, that they're genuine but have been labelled with the wrong library?
#' i.e. take a look at these ones and assess if this could be a possibility.
#' My guess is probably not as they seem to be at too high abundance for that to be
#' a plausible suggestion?
#' 
#' So why are they cropping up? What is the data generating process that leads to them being
#' there? If we can't describe it, then just filtering based on some abundance proportion
#' might be the best way to go?

#' What about the ones that are more than distance 1 away from their controls?

fk %>% filter(Distance > 1) %>% select(gST, Distance, Parent, Count) %>%
  left_join(all_counts %>% rename(Total=Count)) %>%
  left_join(ctrl1 %>% filter(Count > 100) %>% select(Parent=gST, ParentCount=Count)) %>%
  left_join(all_dist1)

#' These have reasonably clearly been slotted in here when they shouldn't have. Are these
#' errors in labelling of which library is which? Most are highly prevalent in total, and all except
#' one (highly prevalent) is an O serogroup. Also, the counts are low (1 or 2 max)

#' Should we repeat the analysis with a differing quality threshold?? We might expect fewer
#' fakers within the controls? YES, this is what happens.
#' 
#' Verdict: There is clear transfer/mislabelling of the library number. What is still unclear
#'          is whether it's one library dominates that transfers to others (or does this
#'          matter?) The error rate could be as large as 1/1000. e.g. id0015557 from Ctrl2.
#'          
#'          In addition, there seems to be evidence of both single SNP error AND library mislabelling
#'          e.g. in id0013775 from Ctrl2 seems like it's most likely a 1 SNP difference plus
#'          library mislabelling. The error rate might be very large again though? e.g. id0009599
#'          in Ctrl 2 suggests 1/100 if all of them are single SNP+library mislabel, though I guess
#'          we have evidence only that one is single SNP+library, so 1/2000. This seems to contradict
#'          the 1/1000 suggestion above, nonetheless, or we have quite a high single SNP error rate,
#'          or they're not independent (which seems unlikely?)
#'
#'          We note that this is 2 SNPs away from O55A though, so could arise if the dual-SNP error
#'          rate was around 1/3000 or so.

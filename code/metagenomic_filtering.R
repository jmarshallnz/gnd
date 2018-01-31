library(readxl)
library(dplyr)

# read in the data, shorten column names a bit
meta <- read_excel("data/metagenomics/metagenomics data.xlsx")

# shorten names a bit
meta <- meta %>% rename(QuerySequence=`Query (1st read and reverse complement of second read joined by string of N's)`,
                GNDmd5=`Hit name`,
                Percent=`Percent identity of match`,
                Length=`Match length`,
                QueryStart=`Query Start`,
                QueryEnd=`Query End`,
                GNDStart=`Hit Start`,
                GNDEnd=`Hit End`)

# find the best match in GND by using various criteria. First percent match
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Percent) %>%
  arrange(Percent)

# Now score (which seems ot be some combination of length and percent maybe?)
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  arrange(Score)

# Note that top_n(1, Foo) returns the top match (or matches, if identical) based on Score.
# A bunch of things match multiple different GNDs (poorly). So basically we don't care much about
# those. We don't care about which GND it matches either, should it match a bunch, so remove GNDmd5
# and kill the non-uniques
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  select(-GNDmd5) %>%
  unique %>%
  arrange(Score)

# We also don't care if there's identical matches at different query start/ends. All we care about is if they match
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  select(-GNDmd5, -QueryStart, -QueryEnd) %>%
  unique %>%
  arrange(Score)

# This still isn't unique, but it gives us a good idea. To figure out why, lets find the duplicates
dat <- meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  select(-GNDmd5, -QueryStart, -QueryEnd) %>%
  unique %>%
  arrange(Score)

dupes <- dat %>%
  select(QuerySequence) %>%
  duplicated %>%
  which

dat %>% filter(QuerySequence %in% dat$QuerySequence[dupes])

# Still not unique - can have same score but different percent by the looks. This is unique
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  select(Sample, QuerySequence, Score) %>%
  unique %>%
  arrange(Score)

# and can choose some arbitrary cutoff if we want. e.g. Scores <100 don't look like very good matches
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  filter(Score < 100) %>%
  select(-GNDmd5) %>%
  unique %>%
  arrange(desc(Score))

# while scores around 150 aren't too bad maybe?
meta %>%
  group_by(QuerySequence) %>%
  top_n(1, Score) %>%
  filter(Score < 150) %>%
  select(-GNDmd5) %>%
  unique %>%
  arrange(desc(Score))

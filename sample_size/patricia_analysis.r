patricia <- read.table('patricia.txt', header=TRUE)
head(patricia)

library(tidyr)
pat <- patricia %>% gather(Result, Count, Positive:Negative)

# expand it outwith random effect for calf etc
expand.grid(c("Positive", "Negative"), value=)

expand_rows <- function(m, column) {
  out <- m[rep(sequence(nrow(m)), m[[column]]), ]
  out[[column]] <- NULL
  out
}

p <- expand_rows(pat, "Count")

library(dplyr)
out <- p %>% group_by(Serogroup, When) %>% mutate(ID = 1:n()) %>% arrange(Serogroup, When, ID)
write.csv(out, "patricia_trans.csv", row.names=FALSE)

pat <- read.csv("patricia_trans.csv") %>% mutate(ID = factor(ID))

pat %>% group_by(Serogroup, When, Result) %>% summarize(Count=n())

library(lme4)

# Hmm, problem: We don't know which isolates correspond to which,
# so we'll be inaccurately summarising the either-or option, right? I guess we could estimate separately?

mod <- glmer(Result ~ When + (1 | ID), data = pat %>% filter(Serogroup=="O157"), family='binomial')
summary(mod)

library(sjstats)
icc(mod) # ICC is 0.999

delogit <- function(x) {
  1/(1+exp(-x))
}
delogit(predict(mod))

mod <- glmer(Result ~ When + (1 | ID), data = pat %>% filter(Serogroup=="O26"), family='binomial')
summary(mod)
icc(mod)

mod <- glm(Result ~ When, data = pat %>% filter(Serogroup=="O26"), family='binomial')

# Not much information from here really... Weird that the glmer doesn't fit well. I guess just 2 observations
# per random effect (most of them will be zero, with just one non-zero...) isn't enough to be any use...

# I guess we could do paired differences instead?
library(samplesize4surveys)

ss4dpH(100000, 0.1, 0.2, 0.05, DEFF = 2, conf = 0.95, power = 0.8, plot = TRUE)

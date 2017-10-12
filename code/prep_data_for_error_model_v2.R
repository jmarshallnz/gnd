#' Generate data for error model.
#' 
#' The model is:
#' 
#' $$
#' E_i \sym Poisson(\lambda_i)
#' $$
#' 
#' where $\lambda_i = q\left[p_i(1-p_i) N^1_i + p_i^2 N^2_i\right] + (1-q)\left[(1-p_i)^2 m_i + p_i(1-p_i)M^1_i + p^iM^2_i\right]
#' 
#' and $N^k_i$ is the number of true gsts within library k away from the i-th gst, m_i is the number
#' of reads in another library of gst i, and $M^k_i$ is the number of true gsts in other libraries
#' $k$ away from gst $i$.
#'
#' First step is to compute $m_i, N^k_i, M^k_i$.
#'
#' Read in abundance data and determine genuines and fakers
library(dplyr)
library(tidyr)
library(cluster)

source("code/read_abundance.R")
source("code/read_fasta.R")

qa <- 3 # quality level

#' Read in abundances
abund <- read_abundance(file=paste0( "temp/sero_abundance", qa, ".csv"))

#' Read in sequences
sequences <- read_fasta(qa = qa)

#' Read in the distances
load(paste0("temp/sero_dist",qa,"_df.Rda"))

ctrl <- abund %>% select(Ctrl1 = Ctrl01R1_S91counts, Ctrl2 = Ctrl02R1_S92counts) %>%
  tibble::rownames_to_column("gST") %>% gather(Control, Count, Ctrl1:Ctrl2)

#' Determine the genuines from the control columns
genuine1 <- ctrl %>% filter(Count > 100, Control=="Ctrl1")
genuine2 <- ctrl %>% filter(Count > 100, Control=="Ctrl2")

#' Determine the fakers from the control columns
faker1 <- ctrl %>% filter(Count <= 100, Control=="Ctrl1") %>% rename(E = Count)
faker2 <- ctrl %>% filter(Count <= 100, Control=="Ctrl2") %>% rename(E = Count)

#' Find out how close each E is to each N, and total up the ones at appropriate distance
genuine1_d1 <- fa.dist.df %>% filter(Distance == 1) %>% left_join(genuine1, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(N1=sum(Count, na.rm=TRUE))
genuine2_d1 <- fa.dist.df %>% filter(Distance == 1) %>% left_join(genuine2, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(N1=sum(Count, na.rm=TRUE))

genuine1_d2 <- fa.dist.df %>% filter(Distance == 2) %>% left_join(genuine1, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(N2=sum(Count, na.rm=TRUE))
genuine2_d2 <- fa.dist.df %>% filter(Distance == 2) %>% left_join(genuine2, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(N2=sum(Count, na.rm=TRUE))

faker1 %>% left_join(genuine1_d1) %>% left_join(genuine1_d2) %>% replace_na(list(N1 = 0, N2 = 0))
faker2 %>% left_join(genuine2_d1) %>% left_join(genuine2_d2) %>% replace_na(list(N1 = 0, N2 = 0))

# OK, now find out how close each E is to each M
# We have to be a bit careful here, as we want to only include non-fakers in m and M from the
# control columns. It's possible that mislabelling of the control columns has happened though,
# so we need to include the genuines there.
altlib1 <- abund %>% select(-Ctrl01R1_S91counts) %>% tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>%
  filter(Library != 'Ctrl02R1_S92counts' | gST %in% genuine2$gST) %>% group_by(gST) %>% summarize(Count=sum(Count))

altlib2 <- abund %>% select(-Ctrl02R1_S92counts) %>% tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>%
  filter(Library != 'Ctrl01R1_S91counts' | gST %in% genuine1$gST) %>% group_by(gST) %>% summarize(Count=sum(Count))

# OK, now find the ones distance 0, 1 and 2 away
altlib1_d1 <- fa.dist.df %>% filter(Distance == 1) %>% left_join(altlib1, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(M1=sum(Count, na.rm=TRUE))
altlib2_d1 <- fa.dist.df %>% filter(Distance == 1) %>% left_join(altlib2, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(M1=sum(Count, na.rm=TRUE))

altlib1_d2 <- fa.dist.df %>% filter(Distance == 2) %>% left_join(altlib1, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(M2=sum(Count, na.rm=TRUE))
altlib2_d2 <- fa.dist.df %>% filter(Distance == 2) %>% left_join(altlib2, by=c('gST2' = 'gST')) %>% group_by(gST) %>% summarize(M2=sum(Count, na.rm=TRUE))

# I think we're done now:

d1 <- faker1 %>%
  left_join(genuine1_d1) %>%
  left_join(genuine1_d2) %>%
  left_join(altlib1 %>% rename(m = Count)) %>%
  left_join(altlib1_d1) %>%
  left_join(altlib1_d2) %>% replace_na(list(N1 = 0, N2 = 0, m = 0, M1 = 0, M2 = 0))

d2 <- faker2 %>%
  left_join(genuine2_d1) %>%
  left_join(genuine2_d2) %>%
  left_join(altlib2 %>% rename(m = Count)) %>%
  left_join(altlib2_d1) %>%
  left_join(altlib2_d2) %>% replace_na(list(N1 = 0, N2 = 0, m = 0, M1 = 0, M2 = 0))

# Now we can fit our model...
dat <- rbind(d1, d2)

E <- dat$E
N1 <- dat$N1
N2 <- dat$N2
m <- dat$m
M1 <- dat$M1
M2 <- dat$M2

# It's a poisson model with beta priors on p and q

log_lik <- function(p, q) {

  if (p < 0 | p > 1)
    return(-Inf)

  if (q < 0 | q > 1)
    return(-Inf)

  # hmm, do we need (1-p) in here. If we do, should it have another multiplier
  # in front of it?
  # The idea of the (1-p) is that the rate of 1 error is p, so 2 errors will be
  # p^2, but 1 error then means (1-p)*p as you've NOT made 2 errors?
  # If the error rate is low, does this actually matter at all?
  lambda <- (1-q) * (p*(1-p)*N1 + p*p*N2) + q*((1-p)*(1-p)*m + p*(1-p)*M1 + p*p*M2)
  # now E ~ Poisson(lambda)
  sum(dpois(E, lambda, log=TRUE))
}

samples  <- 1000
thinning <- 1
burnin   <- 0

total <- burnin + samples*thinning

# prior values
p_a <- 1; p_b <- 100-1 # 1/100 on avg
p <- p_a / (p_a + p_b)
q_a <- 1; q_b <- 100-1 # 1/100 on avg
q <- q_a / (q_a + q_b)

# posterior storage
posterior <- matrix(0, samples, 2)
colnames(posterior) <- c("p", "q")

# accept reject
accept <- matrix(0,2,2)
colnames(accept) <- colnames(posterior)
rownames(accept) <- c("accept", "reject")

p_sd <- 0.01

# try P on a log scale with normal updating
log_p_mu <- log(0.01); log_p_sigma <- 2
log_p <- log_p_mu
prop_p_sd <- 0.2

log_q_mu <- log(0.01); log_q_sigma <- 2
log_q <- log_q_mu
prop_q_sd <- 1

p <- exp(log_p); q <- exp(log_q)
ll <- log_lik(p, q)

n <- 1
for (i in 1:total) {
  # update p. Stupid independence sampler from prior for now (DUMB!)
  #phat <- rbeta(1, p_a, p_b)
  
  # ok, more intelligent now - try MH. Proposal is I guess symmetric (we're
  # kinda relying on it not returning a -ve though...)
  #phat <- rnorm(1, p, p_sd)
  #lh_ratio = dbeta(phat, p_a, p_b, log=TRUE) - dbeta(p, p_a, p_b, log=TRUE)
  
  # ok, on a log scale now
  log_phat <- rnorm(1, log_p, prop_p_sd)

  # prior ratio (no need for proposal ratio though)
  lh_ratio = dnorm(log_phat, log_p_mu, log_p_sigma, log=TRUE) - dnorm(log_p, log_p_mu, log_p_sigma, log=TRUE)

  ll_hat <- log_lik(exp(log_phat), q)
  ll_ratio = ll_hat - ll

  # accept/reject
  log_alpha = ll_ratio + lh_ratio
  if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
    log_p <- log_phat
    p <- exp(log_p)
    ll <- ll_hat
    accept[1,1] <- accept[1,1]+1
  } else {
    accept[2,1] <- accept[2,1]+1
  }

  # update q
#  qhat <- rbeta(1, q_a, q_b)
#  lh_ratio = 0

  # prior ratio and proposal ratios cancel as we're doing a dumb independence sampler
  log_qhat <- rnorm(1, log_q, prop_q_sd)

  # prior ratio (no need for proposal ratio though)
  lh_ratio = dnorm(log_qhat, log_q_mu, log_q_sigma, log=TRUE) - dnorm(log_q, log_q_mu, log_q_sigma, log=TRUE)
  
  ll_hat <- log_lik(p, exp(log_qhat))
  ll_ratio = ll_hat - ll

  # accept/reject
  log_alpha = ll_ratio + lh_ratio
  if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
    log_q <- log_qhat
    q <- exp(log_q)
    ll <- ll_hat
    accept[1,2] <- accept[1,2]+1
  } else {
    accept[2,2] <- accept[2,2]+1
  }

  # sample
  if (i > burnin) {
    if ((i - burnin) %% thinning == 0) {
      posterior[n,] <- c(p,q)
      n <- n + 1
    }
  }
}

# OK, now assume the model is fitting (hah!) How well? What is the posterior fit?
lambda_fit <- function(x) {
  p <- x[1]; q <- x[2]
  (1-q) * (p*(1-p)*N1 + p*p*N2) + q*((1-p)*(1-p)*m + p*(1-p)*M1 + p*p*M2)
}

lambda <- rowMeans(apply(posterior[20:1000,], 1, lambda_fit))
out <- data.frame(gST = dat$gST, E=E, lambda=lambda, P=dpois(E, lambda))
out %>% top_n(-10, P) %>% arrange(P)

# OK, now cluster the data based on kmedoids or whatever
load("temp/sero_dist3.Rda")

# do a clustering
library(cluster)

fa15.dist <- as.matrix(fa15.dist)
# cluster the O-serogroups into 2...
wch <- substring(rownames(fa15.dist),1,1)=="O"
fa15.dist <- as.dist(fa15.dist)
clust2 <- pam(fa15.dist, 2)

clustO <- pam(as.dist(fa15.dist[wch,wch]),2)
df <- data.frame(gST=names(clustO$clustering), Cluster=clustO$clustering)

df <- data.frame(gST=names(clust2$clustering), Cluster=clust2$clustering)
df %>% filter(Cluster == 1)
df %>% filter(Cluster == 1, substring(gST,1,1) == "O") %>% arrange(gST)
write.csv(df, "temp/cluster2.csv", row.names=FALSE)

# OK, a better way to go maybe is to use hclust to cluster
# with single linkage, picking out any that are more than
# say 20 BP away from anything else?
hc <- hclust(as.dist(fa.dist), method="single")
dend <- as.dendrogram(hc)

# TODO: Run down the tree and grab everything that's not below h=20
cutted <- cut(dend, h=10)
# run down the
labels(cutted$lower[[1]])
plot(hc, labels=FALSE)
member_list <- lapply(cutted$lower, labels)
lengths(member_list)

weirdos <- unlist(member_list[-which.max(lengths(member_list))])

totals <- abund %>%
  tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>%
  group_by(gST) %>%
  summarize(Count = sum(Count), maxCount=max(Count))

weirdo_summary <- totals %>% filter(gST %in% weirdos)
fa.dist.df %>% filter(Distance != 0) %>% group_by(gST) %>% top_n(-1, Distance)

seq<-data.frame(sequence=seq, gST=sequences$serogroup, md5=sequences$md5)

weirdo_summary <- weirdo_summary %>% left_join(seq)
write.csv(weirdo_summary, "weirdos.csv", row.names = FALSE)
# NOW: refit model with 2 different error rates.

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

farm <- 1
qa <- 2 # quality level
minTotal <- 0 # minimum abundance level
max_dist_to_parent <- 2 # maximum distance to a parent

#' Read in abundances
abund <- read_abundance(file=paste0("temp/farm", farm, "/sero_abundance", qa, "_", minTotal, ".csv"))

#' Read in sequences
sequences <- read_fasta(qa = qa, minTotal=minTotal)

#' Read in the CLOSE distances
load(paste0("temp/sero_close",qa,"_", minTotal, "_df.Rda"))

ctrl <- abund %>% select(Ctrl1 = Ctrl01R1_S91counts, Ctrl2 = Ctrl02R1_S92counts) %>%
  tibble::rownames_to_column("gST") %>% gather(Control, Count, Ctrl1:Ctrl2)

#' Determine the genuines from the control columns
cut_off <- 100
genuine1 <- ctrl %>% filter(Count > cut_off, Control=="Ctrl1")
genuine2 <- ctrl %>% filter(Count > cut_off, Control=="Ctrl2")

#' Check they're the right length
if (nrow(genuine1) != 7) {
  stop("Don't have 7 genuines in control 1?")
}
if (nrow(genuine2) != 7) {
  stop("Don't have 7 genuines in control 2?")
}

#' Determine the fakers from the control columns
faker1 <- ctrl %>% filter(Count <= cut_off, Control=="Ctrl1") %>% rename(E = Count)
faker2 <- ctrl %>% filter(Count <= cut_off, Control=="Ctrl2") %>% rename(E = Count)

#' Find out how close each E is to each N, and total up the ones at appropriate distance
genuine1_d <- fa.close.df %>% filter(Distance > 0, Distance <= max_dist_to_parent) %>% left_join(genuine1, by=c('gST2' = 'gST')) %>% group_by(gST, Distance) %>%
  summarize(Count=sum(Count, na.rm=TRUE)) %>% rename(N=Distance) %>%
  spread(N, Count, fill=0, sep='')

genuine2_d <- fa.close.df %>% filter(Distance > 0, Distance <= max_dist_to_parent) %>% left_join(genuine2, by=c('gST2' = 'gST')) %>% group_by(gST, Distance) %>%
  summarize(Count=sum(Count, na.rm=TRUE)) %>% rename(N=Distance) %>%
  spread(N, Count, fill=0, sep='')

# OK, now find out how close each E is to each M
# We have to be a bit careful here, as we want to only include non-fakers in m and M from the
# control columns. It's possible that mislabelling of the control columns has happened though,
# so we need to include the genuines there.
altlib1 <- abund %>% select(-Ctrl01R1_S91counts) %>% tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>%
  filter(Library != 'Ctrl02R1_S92counts' | gST %in% genuine2$gST) %>% group_by(gST) %>% summarize(Count=sum(Count))

altlib2 <- abund %>% select(-Ctrl02R1_S92counts) %>% tibble::rownames_to_column("gST") %>% gather(Library, Count, -gST) %>%
  filter(Library != 'Ctrl01R1_S91counts' | gST %in% genuine1$gST) %>% group_by(gST) %>% summarize(Count=sum(Count))

# OK, now find the ones distance 1 and 2 away
altlib1_d <- fa.close.df %>% filter(Distance > 0, Distance <= max_dist_to_parent) %>% left_join(altlib1, by=c('gST2' = 'gST')) %>% group_by(gST, Distance) %>%
  summarize(Count=sum(Count, na.rm=TRUE)) %>% rename(M=Distance) %>%
  spread(M, Count, fill=0, sep='')
altlib2_d <- fa.close.df %>% filter(Distance > 0, Distance <= max_dist_to_parent) %>% left_join(altlib2, by=c('gST2' = 'gST')) %>% group_by(gST, Distance) %>%
  summarize(Count=sum(Count, na.rm=TRUE)) %>% rename(M=Distance) %>%
  spread(M, Count, fill=0, sep='')

# I think we're done now:
replacers <- c(paste0("N", 1:max_dist_to_parent), "m", paste0("M", 1:max_dist_to_parent))
replaced <- numeric(length(replacers))
names(replaced) <- replacers

d1 <- faker1 %>%
  left_join(genuine1_d) %>%
  left_join(altlib1 %>% rename(m = Count)) %>%
  left_join(altlib1_d) %>% replace_na(as.list(replaced))

d2 <- faker2 %>%
  left_join(genuine2_d) %>%
  left_join(altlib2 %>% rename(m = Count)) %>%
  left_join(altlib2_d) %>% replace_na(as.list(replaced))

# Before we model, split into two clusters to find out the weirdos...
# We need only cluster the parents though
parents <- fa.close.df %>% select(gST=gST2) %>% unique %>% left_join(sequences, by=c('gST'='serogroup'))

library(cluster)
parents.d <- daisy(parents %>% mutate_at(vars(starts_with("X")), function(x) as.factor(x)) %>% select(starts_with("X"), gST) %>% tibble::column_to_rownames("gST"))*284

cut_off <- 35
parents.hc <- as.dendrogram(hclust(parents.d, method="complete"))

par(mfrow=c(1,1))
plot(parents.hc, leaflab='none')

cutted <- cut(parents.hc, h=cut_off)
member_list <- lapply(cutted$lower, labels)
largest <- which.max(lengths(member_list))

# check the largest - should be no stragglers
plot(cutted$lower[[largest]], leaflab='none')

# now assign the weirdos in the parents
parent_grouping <- parents %>% mutate(ErrorRate = ifelse(gST %in% member_list[[largest]], 1, 2)) %>%
  select(gST, ErrorRate)

# OK, now assign all gSTs to the group of their closest parent
weirdo_mapping <- fa.close.df %>% group_by(gST) %>% top_n(1, -Distance) %>%
  left_join(parent_grouping %>% rename(gST2=gST)) %>% summarize(ErrorRate = min(ErrorRate), Conflict=length(unique(ErrorRate)) > 1)

# Check there's no conflict
weirdo_mapping %>% filter(Conflict)

dat <- rbind(d1, d2) %>% left_join(weirdo_mapping)

# check there aren't any really odd ones (e.g. ones we have errors for but have never seen...)
so_far_away_from_me <- dat %>% gather(Dist, Count, starts_with("N"), starts_with("M"), "m") %>%
  group_by(gST, Control, E) %>% summarise(Count=max(Count)) %>%
  filter(E > 0, Count == 0) %>%
  pull(gST)

dat <- rbind(d1, d2) %>% left_join(weirdo_mapping) %>% filter(!gST %in% so_far_away_from_me)

E <- dat$E
N <- dat %>% select(starts_with("N")) %>% as.matrix
m <- dat$m
M <- dat %>% select(starts_with("M", ignore.case=FALSE)) %>% as.matrix
WhichP <- dat$ErrorRate

# It's a poisson model with beta priors on p and q

log_lik <- function(p, q) {

  if (any(p < 0 | p > 1))
    return(-Inf)

  if (q < 0 | q > 1)
    return(-Inf)

  P = p[WhichP,]

  # hmm, do we need (1-p) in here. If we do, should it have another multiplier
  # in front of it?
  # The idea of the (1-p) is that the rate of 1 error is p, so 2 errors will be
  # p^2, but 1 error then means (1-p)*p as you've NOT made 2 errors?
  # If the error rate is low, does this actually matter at all?
  lambda <- (1-q) * (rowSums(P*N)) + q*(m + rowSums(P*M))
  # now E ~ Poisson(lambda)
  sum(dpois(E, lambda, log=TRUE))
}

log_lik <- function(p, q, r) {

  if (any(p < 0 | p > 1))
    return(-Inf)

  if (q < 0 | q > 1)
    return(-Inf)

  P = p[WhichP,]

  # hmm, do we need (1-p) in here. If we do, should it have another multiplier
  # in front of it?
  # The idea of the (1-p) is that the rate of 1 error is p, so 2 errors will be
  # p^2, but 1 error then means (1-p)*p as you've NOT made 2 errors?
  # If the error rate is low, does this actually matter at all?
  lambda <- (1-q) * (rowSums(P*N)) + q*(m + rowSums(P*M))

  p <- r / (r + lambda)

  # now E ~ NegBin(r, p)
  sum(dnbinom(E, size=r, prob = p, log=TRUE))
}

samples  <- 1000
thinning <- 10
burnin   <- 100

total <- burnin + samples*thinning

# try P on a log scale with normal updating
log_p_mu <- log(0.001); log_p_sigma <- 2
log_p <- matrix(log_p_mu, 2, max_dist_to_parent)
prop_p_sd <- 0.2

log_q_mu <- log(0.001); log_q_sigma <- 2
log_q <- log_q_mu
prop_q_sd <- 1

r_a <- 0.5; r_b <- 0.1; # gamma params for the dispersion
plot(function(x) { dgamma(x, r_a, rate=r_b) }, xlim=c(0,2))
r <- r_a / r_b # dispersion with a gamma prior
prop_r_sd <- 1;

p <- exp(log_p); q <- exp(log_q)

ll <- log_lik(p, q, r)

# posterior storage
posterior <- matrix(0, samples, 3 + length(p))
colnames(posterior) <- c(paste0('p', apply(expand.grid(1:nrow(p), 1:ncol(p)), 1, paste, collapse='')), "q", "r", "ll")

# accept reject
accept <- matrix(0,2,3) # We're sampling both P's at once
colnames(accept) <- c("p", "q", "r")
rownames(accept) <- c("accept", "reject")

n <- 1
for (i in 1:total) {
  # update p. Stupid independence sampler from prior for now (DUMB!)
  #phat <- rbeta(1, p_a, p_b)
  
  # ok, more intelligent now - try MH. Proposal is I guess symmetric (we're
  # kinda relying on it not returning a -ve though...)
  #phat <- rnorm(1, p, p_sd)
  #lh_ratio = dbeta(phat, p_a, p_b, log=TRUE) - dbeta(p, p_a, p_b, log=TRUE)
  
  # ok, on a log scale now
  
  # sample the p's one at a time
  for (j in 1:length(p)) {
    log_phat <- log_p
    log_phat[j] <- rnorm(1, log_p[j], prop_p_sd)

    # likelihood
    ll_hat <- log_lik(exp(log_phat), q, r)

    # prior ratio (no need for proposal ratio though)
    lh_ratio = dnorm(log_phat[j], log_p_mu, log_p_sigma, log=TRUE) - dnorm(log_p[j], log_p_mu, log_p_sigma, log=TRUE)
    ll_ratio = ll_hat - ll
    if (is.nan(ll_ratio))
      ll_ratio <- 0; # don't accept if we're nowhere in range still

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
  }
  if (0) {
  # NOTE: ATM We're sampling all p's both at once.
  log_phat <- matrix(rnorm(length(log_p), log_p, prop_p_sd), nrow(log_p), ncol(log_p))

  # prior ratio (no need for proposal ratio though)
  lh_ratio = sum(dnorm(log_phat, log_p_mu, log_p_sigma, log=TRUE)) - sum(dnorm(log_p, log_p_mu, log_p_sigma, log=TRUE))
  
  ll_hat <- log_lik(exp(log_phat), q, r)
  ll_ratio = ll_hat - ll
  if (is.nan(ll_ratio))
    ll_ratio <- 0; # don't accept if we're nowhere in range still

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
  }

  # update q
#  qhat <- rbeta(1, q_a, q_b)
#  lh_ratio = 0

  # prior ratio and proposal ratios cancel as we're doing a dumb independence sampler
  log_qhat <- rnorm(1, log_q, prop_q_sd)

  # prior ratio (no need for proposal ratio though)
  lh_ratio = dnorm(log_qhat, log_q_mu, log_q_sigma, log=TRUE) - dnorm(log_q, log_q_mu, log_q_sigma, log=TRUE)
  
  ll_hat <- log_lik(p, exp(log_qhat), r)
  ll_ratio = ll_hat - ll
  if (is.nan(ll_ratio))
    ll_ratio <- 0; # don't accept if we're nowhere in range still

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
  
  # sample r using a lognormal distribution
  r_hat <- rlnorm(1, log(r), prop_r_sd)

  # prior ratio
  lp_ratio = dgamma(r_hat, r_a, rate=r_b, log=TRUE) - dgamma(r, r_a, rate=r_b, log=TRUE)
  
  # proposal ratio
  lh_ratio = dlnorm(r, log(r_hat), prop_r_sd, log=TRUE) - dlnorm(r_hat, log(r), prop_r_sd, log=TRUE)

  ll_hat <- log_lik(p, exp(log_qhat), r_hat)
  ll_ratio = ll_hat - ll

  # TEST. Should be equivalent to
  if (is.nan(ll_ratio))
    ll_ratio <- 0; # don't accept if we're nowhere in range still
  
  # accept/reject
  log_alpha = ll_ratio + lp_ratio + lh_ratio
  if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
    r <- r_hat
    ll <- ll_hat
    accept[1,3] <- accept[1,3]+1
  } else {
    accept[2,3] <- accept[2,3]+1
  }
  

  # sample
  if (i > burnin) {
    if ((i - burnin) %% thinning == 0) {
      posterior[n,] <- c(p,q,r,ll)
      n <- n + 1
    }
  }
  # progress
  if (i %% 100 == 0) {
    cat("Up to iteration", i, "of", total, "\n")
  }
}

# check model traces
par(mfrow=c(ncol(posterior),1), mar=c(2,2,1,2))
apply(posterior[-c(1:50),], 2, function(x) { plot(x, type='l') })
apply(posterior[-c(1:50),], 2, mean)
# summary of p
matrix(apply(posterior[-c(1:50),1:length(p)], 2, mean), nrow(p), ncol(p))
matrix(apply(posterior[-c(1:50),1:length(p)], 2, sd), nrow(p), ncol(p))

# OK, now assume the model is fitting (hah!) How well? What is the posterior fit?
lambda_fit <- function(x) {
  p <- matrix(x[1:length(p)],nrow(p),ncol(p)); q <- x[length(p)+1]; r <- x[length(r)+2];
  P <- p[WhichP,]
  lambda <- (1-q) * (rowSums(P*N)) + q*(m + rowSums(P*M))
  lambda
}

dic <- function(theta, ll) {
  d_bar <- -2*mean(ll)
  d_hat <- -2*sum(dpois(E, lambda_fit(theta), log=TRUE))
  pD = d_bar - d_hat
  d_hat + 2 * pD
}

dic_nb <- function(posterior, burnin=50) {
  # mean of parameters
  theta <- apply(posterior[-seq_len(burnin),], 2, mean)
  p <- matrix(theta[1:length(p)],nrow(p),ncol(p))
  q <- theta[length(p)+1]
  r <- theta[length(p)+2]
  # mean of logliklihood
  ll <- theta[length(theta)]

  d_bar <- -2*ll
  d_hat <- -2*log_lik(p, q, r)
  pD = d_bar - d_hat
  d_hat + 2 * pD
}

dic_nb(posterior, burnin=50)

p_fit <- function(i, posterior) {
  # pull out the parameters from the posterior
  pp <- array(posterior[,1:length(p)], dim=c(nrow(posterior), nrow(p), ncol(p)))
  qp <- posterior[,length(p)+1];
  rp <- posterior[,length(p)+2];
  # OK, now for each one we have to do this vile computation...
  P = pp[,WhichP[i],]
  
  lambda <- (1-q) * (P%*%N[i,]) + q*(m[i] + P%*%M[i,])
  prob <- rp / (rp + lambda)

  # now E ~ NegBin(r, p)
  c(lambda=mean(lambda), r = mean(rp), P = mean(pnbinom(E[i]-1, size=rp, prob = prob, lower.tail = FALSE)))
}

# NOTE: This takes ages and chews up memory. It would be more efficient
#       to compute the fit per item instead.
aa <- lapply(1:nrow(dat), p_fit, posterior=posterior[-c(1:50),])
aaa <- data.frame(dat, do.call(rbind, aa))
aaa %>% top_n(-100, P) %>% arrange(P)

p <- array(posterior[,1:length(p)], dim=c(nrow(posterior), nrow(p), ncol(p)))

a <- lapply(51:1000, function(x) { data.frame(gST=dat$gST, Control=dat$Control, p_fit(posterior[x,])) })
ad <- do.call(rbind, a) %>% group_by(gST, Control) %>%
  summarize(lambda=mean(lambda), r=mean(r), P=mean(p))
out <- dat %>% left_join(ad)
out %>% top_n(-100, P) %>% arrange(P)


lambda <- rowMeans(apply(posterior[20:1000,], 1, lambda_fit))
out <- data.frame(gST = dat$gST, E=E, N, m, M, WP=WhichP, lambda=lambda, P=dpois(E, lambda))
out %>% top_n(-100, P) %>% arrange(P)
write.csv(out, paste0("model_fit_", qa, ".csv"), row.names=FALSE)

# check the weirdos
out %>% filter(WP == 1) %>% top_n(-10, P) %>% arrange(P)




# There's a few here that are way away from anything with any count. e.g.
# id0010100 is 4 away from 4 other IDs with low abundance, and 5 away from O169, abundance 2104.
# id0023628 is 4 away from O145B with abundance 4727... (This is in Ctrl2, so p1 error rate...)
# id0060422 is 3 awat from O111 with abundance 796... (This is in Ctrl2, so p1 error rate...)
d10100 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib1, by=c('gST2'='gST')) %>% head(n=10)
d10100 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib2, by=c('gST2'='gST')) %>% head(n=10)

d23628 <- calc_single_dist("id0023628")
d23628 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib1, by=c('gST2'='gST')) %>% head(n=10)
d60422<- calc_single_dist("id0060422")
d60422 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib1, by=c('gST2'='gST')) %>% head(n=10)

d10100 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib2, by=c('gST2'='gST')) %>% head(n=10)

d11966 <- calc_single_dist("id0011966")

d11966 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib2, by=c('gST2'='gST')) %>% head(n=10)
d11966 %>% mutate(Distance = as.numeric(as.character(Distance))) %>%
  arrange(Distance) %>% left_join(altlib1, by=c('gST2'='gST')) %>% head(n=10)



# OK, now cluster the data based on kmedoids or whatever
load("temp/sero_dist3.Rda")

# do a clustering
library(cluster)

fa.dist <- as.matrix(fa.dist)
# cluster the O-serogroups into 2...
wch <- substring(rownames(fa.dist),1,1)=="O"
clust2 <- pam(as.dist(fa.dist), 2)

clustO <- pam(as.dist(fa.dist[wch,wch]),2)
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

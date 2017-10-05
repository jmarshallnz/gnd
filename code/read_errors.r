#' Some stuff for figuring out read errors.
#' 
#' To do this, we assume that read errors arise (independently) with probability e per
#' base, which will be small. We assume e is the same for each of the 284 bases and that
#' errors are independent, and that at most one error occurs at each base (and all possible
#' state transistions at the base have the same error rate).
#' 
#' In addition, we assume that all 'true' sequences in the sample have been observed.
#' 
#' Then we get that the reads we get arise through a multinomial process based on the true
#' proportion of serogroups in the sample (and total number of reads).
#' 
#' Each of the reads then has a chance of being a bit wrong, so what ends up in the sample
#' is the sum of the correct number of reads plus the error ones which are supposed to be
#' something else.
#' 
#' The idea is that hopefully some of the serogroups that are really rare are more likely to
#' have arisen from the small error rate than be actually present in really small abundance in
#' the sample.
#' 
#' To do this, we put a quite restrictive prior on the prevalences, so that most of the 'unknown'
#' serogroups have a really low chance of cropping up apriori. This will then be dominated by the
#' data with any luck, to give us information that we can use.
#' 
#' Note that this system is under-determined. We have e plus the n (number of serogroups) proportions
#' to estimate, from just the n data points. Strong-ish priors on e and p will help!
#' 
#' The maths turns out to be quite nice if you condition the likelihood on the true number of reads
#' in serogroup i that are actually supposed to be serogroup j with error (or j = i with no error).
#' It turns out that we can Gibbs sample then - basically it's like the latent class problem for
#' testing in the absence of a gold standard. but using multinomial and dirichlet rather than
#' binomial and beta. We'll need to watch the mixing, and without strong priors we may not explore
#' the posterior space very well if the mixing is poor.

#' NOTE: Run create_distance_matrix.R first

source("code/read_abundance.R")

#' The sample. Can optionally include controls here (e.g. for MDS plot in the paper)
appendix <- ""
#y = rowSums(read_abundance())
#d = as.matrix(read.csv("temp/sero_dist15.csv", row.names=1))
#d = d[names(y), names(y)]

#' A dummy sample. Consists of a control column and a second column
#' with 1 known truth, and 4 others
y = matrix(0, 4, 2)
y[,1] = c(100, 10, 5, 5)
y[,2] = c(5, 5, 100, 100)
d <- 2 - diag(2, 4, 4)
d[1,2] <- d[2,1] <- 1
d[1,3] <- d[3,1] <- 1
d[3,4] <- d[4,3] <- 1

#' setup the structure of stuff based on genetic distance
dn <- apply(d, 2, function(x) { which(x == 1) })

#' Number of observed serogroups
n = nrow(y)

#' Dirichlet prior on the unknown prevalences p
prior_p = rep(0.00001, n)

#' Beta prior on the error rate e
prior_e = c(40000,1000000)

#' MCMC control
iters  = 5
burnin = 1
thin   = 1

#' Storage for the posteriors
post_e = rep(NA, iters/thin)
post_f = rep(NA, iters/thin)
post_p = matrix(NA, iters/thin, n)
post_x = list() # control latent vars
post_z = list()  # 'other' latent vars (which in future will be some over all other cols)

#' Latent variables for the truth (assume our sample is truth to start)
nx <- lengths(dn)
x <- lapply(seq_len(nrow(y)), function(n) { c(y[n,1], rep(0,nx[n])) })
wx <- lapply(seq_len(nrow(y)), function(n) { c(n, dn[[n]])})
sx <- y[,1] #' sum of x

z <- lapply(seq_len(nrow(y)), function(n) { c(y[n,2], rep(0,nx[n])) })
sz <- y[,2] #' sum of z

#' Current values
e = 0.01 # error rate within library
f = 0.01 # error rate between libraries (TODO: Does it need to be multiplicative with e?)
p = sweep(y, 2, FUN='/', colSums(y))

#' dirichlet function
rdirichlet<-function(n,a)
{
  l  <- length(a);
  x  <- matrix(rgamma(l*n,a), ncol=l, byrow=TRUE);
  sm <- x %*% rep(1,l);
  x / as.vector(sm);
}

#' MCMC loop
j = 1
for (l in 1:(iters + burnin)) {
  
  #' 1. Gibbs update x for the first column
  sx = rep(0, length(x))
  for (i in 1:length(x)) {
    #' x_{i.} = Multinomial(y_i, ...)
    xp = as.numeric(rmultinom(1, y[i,1], c(1-e, rep(e, nx[i]))*p[wx[[i]],1]))
    x[[i]] = xp
    # also update sx[i] (sum of x i)
    sx[wx[[i]]] = sx[wx[[i]]] + xp
  }
  #' 1b Gibbs update z for other columns (assumed to be truth ATM)
  sz = rep(0, length(z))
  for (i in 1:length(z)) {
    #' z_{i.} = Multinomial(y_i, ...)
    zp = as.numeric(rmultinom(1, y[i,1], c(1-f, rep(f, nx[i]))*p[wx[[i]],2]))
    z[[i]] = zp
    # also update sz[i] (sum of z i)
    sz[wx[[i]]] = sz[wx[[i]]] + zp
  }


  #' 2. Gibbs update p
  p[,1] = rdirichlet(1, sx + prior_p)
  p[,2] = rdirichlet(1, sz + prior_p)

  #' 2b TODO: Gibbs update p for other columns (assumed to be truth ATM)

  #' 3. Gibbs update e
  s0x = sum(unlist(lapply(x, function(y) { y[1] })))
  s1x = sx - s0x

  e = rbeta(1, s0x + prior_e[1], s1x + prior_e[2])
  f = rbeta(1, s0x + prior_e[1], s1x + prior_e[2])

  #' sample...
  if (l > burnin && (l-burnin) %% thin == 0) {
    cat("Done iteration", l, "of", iters+burnin, "\n")
    post_p[j,] = p
    post_e[j] = e
    post_x[j,,] = x
    j = j + 1;
  }
}
apply(post_x,2:3,sum)
apply(post_p,2,mean)*sum(y)
mean(post_e)
hist(post_e)

#' Compute posterior mean of latent variables, scaled to be a proportion from each type
post_x_means = apply(post_x,2:3,mean)
rownames(post_x_means) = colnames(post_x_means) = rownames(d)
post_x_means = post_x_means / rowSums(post_x_means)

#' Find those types that are less than 50% chance of being real (diagonal will
#' be those that arise correctly as this type..)
rows = diag(post_x_means) < 0.5

#' Figure out how many types we'd have if we adjust the 50% cut-off above.
conf = seq(0,1,length.out=100)
count_errors = function(px, conf) {
  cat("doing conf=", conf, "\n")
  sum(diag(px/rowSums(px)) <= conf)
}

num_conf = lapply(conf, function(y) { apply(post_x, 1, count_errors, y) })
num_conf = as.data.frame(num_conf)
names(num_conf) <- 1:100
pdf(paste0("error_rates_by_cutoff", appendix, ".pdf"), width=8, height=6)
plot_conf = apply(num_conf, 2, quantile, c(0.025, 0.5, 0.975))
plot(100*(1-conf), plot_conf[2,], ylim=c(0,403), type="l", ylab="Number of serogroups", xlab="Cut-off for arising from some other serogroup in error (%)")
lines(100*(1-conf), plot_conf[1,], lty="dashed")
lines(100*(1-conf), plot_conf[3,], lty="dashed")
dev.off()

heat_ma_map = function(mat) {
  library(RColorBrewer)
  pal = brewer.pal(9, "YlGnBu")
  par(mar = c(6,6,2,2))
  image(1:nrow(mat), 1:ncol(mat), mat, col=pal, xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, 1:nrow(mat), rownames(mat), las=2, cex.axis=0.8)
  axis(2, 1:ncol(mat), colnames(mat), las=2, cex.axis=0.8)
  mtext("Probable source", side=2, line=5)
  mtext("In sample", side=1, line=5)
}
mat_row = post_x_means[rows,]
mat_red = mat_row[,colSums(mat_row) > 0.5]

pdf(paste0("read_errors_maybe", appendix, ".pdf"), width=24, height=5)
heat_ma_map(mat_red)
dev.off()

# most likely source
most_likely_source <- colnames(mat_row)[apply(mat_row, 1, which.max)]

mapping = data.frame(sample = rownames(mat_row), probable_source = most_likely_source)
write.csv(mapping, paste0("error_mapping", appendix, ".csv"), row.names=FALSE)


mapping = read.csv(paste0("error_mapping", appendix, ".csv"))

fa15 <- read_fasta()

#' compute difference maps and higlight them
map_source <- mapping %>% left_join(fa15 %>% select(-md5), by=c('sample' = 'serogroup'))
map_dest   <- mapping %>% left_join(fa15 %>% select(-md5), by=c('probable_source' = 'serogroup'))

# now check the difference between them and highlight it...
diff = map_source[,-(1:2)] != map_dest[,-(1:2)]
map_diff = cbind(map_source[,1:2], diff)

map_diff = map_diff %>% arrange(probable_source)
map_matrix = as.matrix(map_diff[,-(1:2)])
# add in an alternating thingee
map_matrix = map_matrix + 2*(as.numeric(map_diff[,2]) %% 2)
colnames(map_matrix) = 1:284
rownames(map_matrix) = map_diff[,1]
pdf(paste0("error_maps_by_base_pairs", appendix, ".pdf"), width=10, height=20)
par(mar=c(4,6,2,6))
image(1:ncol(map_matrix), 1:nrow(map_matrix), t(map_matrix), col=c("white", "black", "grey80", "black"), xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, 1:nrow(map_matrix), rownames(map_matrix), las=2, cex.axis=0.6)
axis(4, 1:nrow(map_matrix), map_diff[,2], las=2, cex.axis=0.6)
dev.off()

# purity computation
impurity <- function(x) {
  p = prop.table(table(x))
  1 - sum(p^2)
}
# need to compute purity based on those that aren't errors
no_error_serogroups <- fa15 %>% filter(!(serogroup %in% map_diff$sample))
error_rates_by_base <- data.frame(errors = colSums(map_diff[,-c(1:2)]), impurity = apply(no_error_serogroups[,1:284], 2, impurity), triple=c(rep(letters[1:3],284/3),letters[1:2]))
plot(jitter(errors) ~ jitter(impurity, factor=20), col=triple, data=error_rates_by_base)

write.csv(error_rates_by_base, paste0("error_rates_by_base", appendix, ".csv"), row.names=FALSE)

# now compute the true abundances. This is standalone

abundance = read_abundance(removed=c(97,98,120))
abundance$serogroup = rownames(abundance)

mapping = read.csv(paste0("error_mapping", appendix, ".csv"), stringsAsFactors = FALSE)

new_abund = abundance %>%
  left_join(mapping, by = c("serogroup" = "sample")) %>%
  mutate(probable_source = factor(ifelse(is.na(probable_source), serogroup, probable_source)))

new_abund2 = apply(new_abund %>% select(-serogroup, -probable_source), 2, function(x) { tapply(x, new_abund$probable_source, sum) })
new_abund2 = data.frame(serogroup = levels(new_abund$probable_source), new_abund2)
write.csv(new_abund2, paste0("no_error_abundance", appendix, ".csv"), row.names=FALSE)

# do some clustering for fun
abund_per_sample = t(new_abund2[,-1]) / colSums(new_abund2[,-1])
abund_dist = dist(abund_per_sample)

abund.hc1 = hclust(abund_dist, method='complete')
plot(abund.hc1)

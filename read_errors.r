#' Some stuff for figuring out read errors.
#' 
#' To do this, we assume that read errors arise (independently) with probability e,
#' which will be small. We assume e is the same for each of the 284 bases, which is
#' possibly complete rubbish, and may be able to be improved upon later.
#' 
#' Then we assume that the reads we get arise through a multinomial process based on the true
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
#' to estimate, from just the n data points. Hopefully strong-ish priors on e and p will help!
#' 
#' The maths turns out to be quite nice if you condition the likelihood on the true number of reads
#' in serogroup i that are actually supposed to be serogroup j with error (or j = i with no error).
#' It turns out that we can Gibbs sample then - basically it's like the LCA problem but using multinomial
#' and dirichlet rather than binomial and beta. The mixing might be a bit sucky though (and without
#' strong priors we may not explore the posterior space very well perhaps?)

#' the sample
y = c(1000, 200, 10, 1)

#' The distances
d = 10*(1 - diag(1, 4))
d[1,3] = d[3,1] = 3
d[2,4] = d[4,2] = 1

#' The real sample
d = as.matrix(round(read.csv("sero_dist15.csv", row.names=1) * 284))
y = as.matrix(read.csv("sero_abundance.csv", row.names=1))

#' TODO: y needs to be not the abundance across everything, rather the abundance per sample

n_serogroups = nrow(y)
n_samples    = ncol(y)
K = max(d)

#' Dirichlet prior on the unknown prevalences p
prior_p = matrix(0.00001, n_serogroups, n_samples)

#' Beta prior on the error rate e, constant for all libraries/lanes
prior_e = c(1,1)

#' MCMC control
iters  = 10000
burnin = 1000
thin   = 50

#' Storage for the posteriors
post_e = rep(NA, iters/thin)
post_p = array(NA, c(iters/thin, n_samples, n_serogroups))
post_x = array(NA, dim=c(iters/thin, n_samples, n_serogroups, n_serogroups))

#' Latent variables for the truth (assume our sample is truth to start)
x = array(NA, dim=c(n_samples, n_serogroups, n_serogroups))
for (i in 1:n_samples) {
  x[i,,] = diag(y[,i]) # TODO: hmm, the way around we have it might need some work...
}

#' Current values
e = 0.01
p = y / colSums(y) # TODO hmm, I think we need Y around the other way for this?
ed = (1-e)^(K-d) * e^d

#' dirichlet function
rdirichlet<-function(n,a)
{
  l  <- length(a);
  x  <- matrix(rgamma(l*n,a), ncol=l, byrow=TRUE);
  sm <- x %*% rep(1,l);
  x / as.vector(sm);
}

#' TODO: we're going to extend x to be per-library (i.e. per 96-cell sample)
#'       and work out prob of error per genotype per cell (they may be correlated across samples, but let's ignore for now!)
#'       
#'       So, p needs to extend to be a matrix, and x to be an array
#'       
#' MCMC loop
j = 1
for (l in 1:(iters + burnin)) {

  #' 1. Gibbs update x
  for (n in 1:n_samples) {
    for (i in 1:n_serogroups) {
      #' x_{i.} = Multinomial(y_i, ...)
      x[n,i,] = rmultinom(1, y[i,n], ed[i,]*p[,n])
    }
  }

  #' 2. Gibbs update p
  for (n in 1:n_samples) {
    p[,n] = rdirichlet(1, colSums(x[n,,]) + prior_p[,n])
  }

  #' 3. Gibbs update e
  sdx = 0;
  skdx = 0;
  for (n in 1:n_samples) {
    sdx = sdx + sum(d*x[n,,])
    skdx = skdx + sum((K-d)*x[n,,])
  }
  e = rbeta(1, sdx + prior_e[1], skdx + prior_e[2])
  ed = (1-e)^(K-d) * e^d

  #' sample...
  cat("Done iteration", l, "of", iters+burnin, "\n")
  if (l > burnin && (l-burnin) %% thin == 0) {
    post_p[j,,] = p
    post_e[j] = e
    post_x[j,,] = x
    j = j + 1;
  }
}
#apply(post_x,2:3,sum)
#apply(post_p,2,mean)*sum(y)
mean(post_e)

# work out which ones are likely errors (i.e. which ones have x switching all the time?)
post_x_mean <- apply(post_x,2:4,mean)[1,,]
rownames(post_x_mean) = colnames(post_x_mean) = rownames(d)
post_x_mean = post_x_mean / rowSums(post_x_mean)

# filter out those that are different
rows = diag(post_x_mean) < 0.5

conf = seq(0,1,length.out=100)
num_conf = unlist(lapply(conf, function(x) { sum(diag(post_x_mean) <= x) }))
plot(100*(1-conf), num_conf, type="l", ylab="Number of serogroups", xlab="Cut-off for arising from some other serogroup in error (%)")
hist(100*diag(post_x_mean), xlab="Probability of error", ylab="Number of serogroups", main="Likelihood of errors in serogroups", col="gray70")

# repeat this for each posterior iteration instead
px = post_x[,1,,]

conf = seq(0,1,length.out=100)
count_errors = function(px, conf) {
#  cat(dim(px))
#  print(rowSums(px))
  sum(diag(px/rowSums(px)) <= conf)
}

num_conf = lapply(conf, function(y) { apply(px, 1, count_errors, y) })
num_conf = as.data.frame(num_conf)
names(num_conf) <- 1:100
plot_conf = apply(num_conf, 2, quantile, c(0.025, 0.5, 0.975))
plot(100*(1-conf), plot_conf[2,], ylim=c(0,403), type="l", ylab="Number of serogroups", xlab="Cut-off for arising from some other serogroup in error (%)")
lines(100*(1-conf), plot_conf[1,], lty="dashed")
lines(100*(1-conf), plot_conf[3,], lty="dashed")

hist(100*diag(post_x_mean), xlab="Probability of error", ylab="Number of serogroups", main="Likelihood of errors in serogroups", col="gray70")

heat_ma_map = function(mat) {
  library(RColorBrewer)
  pal = brewer.pal(9, "YlGnBu")
  par(mar = c(6,6,2,2))
  image(1:ncol(mat), 1:nrow(mat), t(mat), col=pal, xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, 1:ncol(mat), colnames(mat), las=2, cex.axis=0.8)
  axis(2, 1:nrow(mat), rownames(mat), las=2, cex.axis=0.8)
  mtext("Probable source", side=1, line=5)
  mtext("In sample", side=2, line=5)
}
mat_row = post_x_mean[rows,]
mat_red = mat_row[,colSums(mat_row) > 0]

pdf("read_errors_maybe.pdf", width=12, height=10)
heat_ma_map(mat_red)
dev.off()

# most likely source
most_likely_source <- colnames(mat_row)[apply(mat_row, 1, which.max)]

mapping = data.frame(sample = rownames(mat_row), probable_source = most_likely_source)

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
par(mar=c(4,6,2,6))
image(1:ncol(map_matrix), 1:nrow(map_matrix), t(map_matrix), col=c("white", "black", "grey80", "black"), xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, 1:nrow(map_matrix), rownames(map_matrix), las=2, cex.axis=0.8)
axis(4, 1:nrow(map_matrix), map_diff[,2], las=2, cex.axis=0.8)


# Plot posteriors
par(mfrow=c(2,2))
plot(post_p[,1], type="l", ylim=c(0,1), main="Abundance traces")
for (i in 2:n)
  lines(post_p[,i], col=i)
d <- list(); xlim = NULL; ylim = NULL
for (i in 1:n) {
  d[[i]] <- density(post_p[,i])
  if (d[[i]]$bw < 0.001)
    d[[i]] <- density(post_p[,i], bw=0.001)
  xlim = range(xlim, d[[i]]$x)
  ylim = range(ylim, d[[i]]$y)
}
plot(d[[1]], ylim=ylim, xlim=xlim, main="Abundance density")
for (i in 2:n) {
  lines(d[[i]], col=i)
  # and the marginal prior
  x = seq(xlim[1], xlim[2], length.out=100)
  y = dbeta(x, prior_p[i], sum(prior_p) - prior_p[i]) * nrow(post_p)
  lines(x, y, col=i, lty="dotted")
}
plot(post_e, type="l", main="Error rate traces")
plot(density(post_e), main="Error rate density")

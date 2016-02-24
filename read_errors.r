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
y = read.csv("sero_abundance.csv", row.names=1)$abundance

n = length(y)
K = max(d)

#' Dirichlet prior on the unknown prevalences p
prior_p = rep(0.00001, n)

#' Beta prior on the error rate e
prior_e = c(1,1)

#' MCMC control
iters  = 10000
burnin = 1000
thin   = 50

#' Storage for the posteriors
post_e = rep(NA, iters/thin)
post_p = matrix(NA, iters/thin, n)
post_x = array(NA, dim=c(iters/thin, n, n))

#' Latent variables for the truth (assume our sample is truth to start)
x = diag(y)

#' Current values
e = 0.01
p = y / sum(y)
ed = (1-e)^(K-d) * e^d

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
  
  #' 1. Gibbs update x
  for (i in 1:nrow(x)) {
    #' x_{i.} = Multinomial(y_i, ...)
    x[i,] = rmultinom(1, y[i], ed[i,]*p)
  }

  #' 2. Gibbs update p
  p = rdirichlet(1, colSums(x) + prior_p)

  #' 3. Gibbs update e
  e = rbeta(1, sum(d*x) + prior_e[1], sum((K-d)*x) + prior_e[2])
  ed = (1-e)^(K-d) * e^d

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

mat = apply(post_x,2:3,mean)
rownames(mat) = colnames(mat) = rownames(d)
mat = mat / rowSums(mat)

conf = seq(0,1,length.out=100)
num_conf = unlist(lapply(conf, function(x) { sum(diag(mat) <= x) }))
plot(conf, num_conf, type="l")

# filter out those that are different
rows = diag(mat) < 0.5

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
mat_row = mat[rows,]
mat_red = mat_row[,colSums(mat_row) > 0]
pdf("read_errors_maybe.pdf", width=12, height=10)
heat_ma_map(mat_red)
dev.off()


mapping = cbind(sample = rownames(mat)[rows], probable_source = colnames(mat)[mat.source[rows]])
write.csv(mapping, "read_errors.csv", row.names=FALSE)

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

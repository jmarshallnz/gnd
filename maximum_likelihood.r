#' exploration of maximum likelihood estimate for p and e
#' 
#' the basic idea is that y ~ Multinom(n, q)
#' 
#' where Ep = q, E = (1-e)^(max(d) - d)*e^d
#'
#' Thus, we can just use \hat{q} = y/n and then \hat{p} = E^{-1} \hat{q} for any given e

d = as.matrix(round(read.csv("sero_dist15.csv", row.names=1) * 284))
y = rowSums(read.csv("sero_abundance.csv", row.names=1))

e = seq(0, 0.01, by = 0.0001)

max_lik_est <- function(e, d, y) {
  E = (1-e)^(max(d) - d)*e^d
  p = solve(E) %*% y/sum(y)
  p
}

non_pos_prev <- function(e, d, y) {
  p = max_lik_est(e, d, y)
  wch = p <= 0
  wch
}

neg_prev = simplify2array(lapply(e, non_pos_prev, d, y))[,1,]
colnames(neg_prev) <- e
plot(e, nrow(neg_prev) - colSums(neg_prev), type="l", xlab="Error rate", ylab="Number of serogroups with p>0")

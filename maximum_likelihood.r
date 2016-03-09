#' exploration of maximum likelihood estimate for p and e
#' 
#' the basic idea is that y ~ Multinom(n, q)
#' 
#' where Ep = q, E = (1-e)^(max(d) - d)*e^d
#'
#' Thus, we can just use \hat{q} = y/n and then \hat{p} = E^{-1} \hat{q} for any given e

#' Hmm, maybe we can think of it as:
#' 
#' p(sequence X) = p(sequence X | no read errors)p(no read errors) +
#'                 p(sequence X | 1 read error)p(1 read error) + 
#'                 ..
#'                 p(sequence X | 62 read error)p(62 read error) + 
#'                 p(sequence X | >=63 read error)p(>=63 read error)
#'                 
#'               = p(X) (1-e)^284 +
#'                 p(1 diff from X) choose(284,1) (1-e)^283 e^1
#'                 ..
#'                 p(62 diff from X) (1-e)^284-62 e^62
#'                 p(>=63 diff from X) \sum (1-e)
#'
#' NOTE: Neither of these are right, as we've observed no other types
#'       which is impossible under this model, as p(other type exists) < 0
#'       for any e > 0. This is because solve(E) is negative for the other
#'       column(s) (TODO: all the time? Can we prove this? Do we need to extend
#'       to 4 bases or does 2 suffice?)
#'       
#'       I suspect we can show it just by having a single, unobserved SLV
#'       Then we show that solve(E) has it's largest +ve entry for the 0 and
#'       a negative entry for the observed SLV.
#'       
#'       Thus, we have to normalise in some way, but that will also
#'       mean that eventually those with low abundance will go negative
#'       as well. (TODO: all the time? Can we prove this? Might it be
#'       able to be shown that as e->0.5, p(max) -> inf?)
#'
#'       So, our probability model of independent read errors is basically
#'       wrong in at least _some_ way. We just don't get enough other stuff.
#'       Maybe as n -> inf though we would? i.e. is it a small-sample problem?
#'       No! That means e is really small, and we know that e isn't that
#'       small, as we observe at least _some_ errors. Maybe this is an
#'       issue with the '10' minimum thingee? If we included _everything_
#'       then perhaps it wouldn't be quite so bad? We would expect a shitload
#'       of extras though given the abundance and number of supposed
#'       errors we're making at the moment...
#'
#'       Some interesting things... If p_i = 1/M for all i then hat(p_i) = 1/M
#'                                  also happens if a suitably chosen subset is the same...
#'      
#'       I guess we could consider trying to maximise the original log-likelihood
#'       given e for p (on the appropriate scale).
#'       
#'       This can be done by maximising log-likelihood which will only
#'       require sums of log(q_i) which isn't really nice for optimisation
#'       as q_i is a matrix product (i.e. a sum)
#'       
#'       
#' p(type 1) = p(type 1 | really type 1)p(really type 1) + 
#'             p(type 1 | really type 2)p(really type 2) +
#'             ..
#'             p(type 1 | really type 403)p(really type 403) + 
#'             p(type 1 | some other type)p(some other type)
#'           = (1-e)^284 (??) * p_1 + 
#'             (1-e)^283*e (??) * p_2 +
#'             ..
#'             (1-e)^222*e^62 * p_403 +
#'             (1 - sum(above_errors)) * p_other
#'             
#'             cbind(E, 1-rowSums(E)) * c(p, p_other)
#'             
#'             (sum(above_errors) - 1)*p_1 +
#'             (sum(above_errors) - 1)*p_2 +
#'             ..
#'             (sum(above_errors) - 1)*p_403
#'             
#'           = (E + rowSums(E) - 1) * p = q - (1-rowSums(E))

#' p(type 1) = p(type 1 | really type 1)p(really type 1)/(sum(p(type_k | really type 1))) + 
#'             p(type 1 | really type 2)p(really type 2)/(sum(p(type_k | really type 2))) +
#'             ..
#'             p(type 1 | really type M)p(really type 2)/(sum(p(type_k | really type M)))

#' We want to scale each row and column by a_i
#' So each entry i,j gets scaled by a_i a_j
#' such that the row sum and column sum are 1
#' how do we do that?
#'  sum_i(a_i a_j x_ij) = 1
#'  sum_j(a_i a_j x_ij) = 1
#'  
#' An 'easy' alternate is just doing 1-rowSum on diagonal...
#' but then I'm not sure the probabilities work out properly.
#' e.g. in the case we observe lots of SLVs we'd have lots of errors
#'      so then rowsum would be greater than 1?
#'      I think it would be biasing to error?
#'
max_lik_est_improved <- function(e, d, y) {
  E = (1-e)^(max(d)-d)*e^d
  Ehat = t(t(E) / colSums(E)) #' This is not really right though, right??
  #' equiv to Ehat = E / colSums(E)
#  Ehat = (E + rowSums(E) - 1)
#  q = y/sum(y) - (1 - rowSums(E))
  q = y/sum(y)
  p = solve(Ehat) %*% q
  p
}

d = as.matrix(round(read.csv("sero_dist15.csv", row.names=1) * 284))
y = rowSums(read.csv("sero_abundance.csv", row.names=1))

e = seq(0, 0.0001, by = 0.000001)*100

max_lik_est <- function(e, d, y) {
  E = (1-e)^(284 - d)*e^d
#  E = E / rowSums(E)
#  E_other = 1 - rowSums(E)
#  E = cbind(E, E_other)
#  E = rbind(E, c(E_other, (2^max(d)-nrow(E))/2^max(d)))
  
  p = solve(E) %*% y/sum(y)
  p
}

non_pos_prev <- function(e, d, y) {
  p = max_lik_est_improved(e, d, y)
  wch = p <= 0
  wch
}

neg_prev = simplify2array(lapply(e, non_pos_prev, d, y))[,1,]
p_hat = simplify2array(lapply(e, max_lik_est_improved, d, y))[,1,]
colnames(neg_prev) <- e
plot(e, nrow(neg_prev) - colSums(neg_prev), type="l", xlab="Error rate", ylab="Number of serogroups with p>0")

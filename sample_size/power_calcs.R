library(ggplot2)

powa <- function(n1=300, n2=300, p1=0.45, p2=0.55, deff=NA, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  if (is.na(deff))
    deff <- 1 + ((cv^2 + 1) * cs - 1) * icc
  z <- qnorm(1 - alpha/2)
  q1 <- 1 - p1
  q2 <- 1 - p2
  pm <- (n1 * p1 + n2 * p2)/(n1 + n2)
  ds <- z * sqrt(deff*(1/n1 + 1/n2) * pm * (1 - pm))
  ex <- abs(p1 - p2)
  sd <- sqrt(deff*(p1 * q1/n1 + p2 * q2/n2))
  1 - pnorm((ds - ex)/sd) + pnorm((-ds - ex)/sd)
}

ss1d <- function(p0=0.45, p_diff=0.1, deff=NA, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  if (is.na(deff))
    deff <- 1 + ((cv^2 + 1) * cs - 1) * icc
  z <- qnorm(1 - alpha/2)

  n <- p0 * (1 - p0) * deff * z^2 / p_diff^2
  n
}

effsize1d <- function(n=300, p0=0.45, deff=NA, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  if (is.na(deff))
    deff <- 1 + ((cv^2 + 1) * cs - 1) * icc
  z <- qnorm(1 - alpha/2)
  
  p_diff = z*sqrt(p0*(1-p0)*deff/n)
  p_diff
}

# now a computation to figure out minimal delta_p detected given power=0.8
powa_dp <- function(delta_p, n1=300, n2=300, p1=0.45, power=0.8, deff=1, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  powa(n1, n2, p1, p1-delta_p, deff, cs, icc, cv, alpha) - power
}

delta_pw <- function(n1=300, n2=300, p1=0.45, power=0.8, deff=NA, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  tryCatch(
    uniroot(powa_dp, interval=c(1e-3, p1), n1=n1, n2=n2, p1=p1, power=power, deff=deff, cs=cs, icc=icc, cv=cv, alpha=alpha)$root,
    error=function(e) { NA })
}

pow_range <- expand.grid(n=seq(200,1000,by=200), p1=seq(0.10,0.40,by=0.01), power=c(0.6,0.8),
                         Deff=seq(2,5,by=1), alpha=0.05)
delta_p <- numeric(nrow(pow_range))
for (i in 1:nrow(pow_range)) {
  pr <- pow_range[i,]
  delta_p[i] <- delta_pw(n1=pr$n, n2=pr$n, p1=pr$p1, power=pr$power, deff=pr$Deff, alpha=pr$alpha)
}
pow_range$delta_p = delta_p
#write.csv(pow_range, "data/sample_size_test.csv", row.names=FALSE)
pow_range$Deff = paste0("Deff=", pow_range$Deff)
pow_range$num = factor(paste0("n=", pow_range$n), levels=paste0("n=", seq(200,1000, by=200)))

# add on prevalences from Patricia's work.
# These are very much best-guess. They're actually probably lower than this...
pre_val <- 13/60
post_val <- 6/60
pre_est <- broom::tidy(binom.test(13, 60))
pre_est2 <- broom::tidy(binom.test(13, 60, conf.level=0.5))

pre_est$conf.low
pre_est$conf.high
post_est <- broom::tidy(binom.test(6, 60))
broom::tidy(prop.test(x=c(13, 6), n=c(60,60)))

band <- data.frame(x=c(pre_est$conf.low, pre_est$conf.high), ymin=0, ymax=max(pow_range$delta_p, na.rm=TRUE)+0.02)
band2 <- data.frame(x=c(pre_est2$conf.low, pre_est2$conf.high), ymin=0, ymax=max(pow_range$delta_p, na.rm=TRUE)+0.02)

ggplot(pow_range) + geom_line(aes(x=p1,y=delta_p,col=factor(power))) +
  facet_grid(Deff~num) +
  theme_bw() +
  scale_color_manual(name='Power', values=c('black', 'red')) +
  ylab("Detectable reduction in prevalence") +
  xlab("Pre-intervention prevalence") + 
  geom_vline(xintercept=13/60, linetype='dotted') + 
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), data=band, alpha=0.2) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), data=band2, alpha=0.2) +
  geom_hline(yintercept=7/60, linetype='dotted') +
  scale_y_continuous(expand=c(0,0))


# now 1-d stuff. We want to work out minimal difference to detect given sample sizes etc.
pow_range <- expand.grid(n=seq(200,400,by=100), p1=seq(0.01,0.4,by=0.01),
                         Deff=seq(1.2,2.4,by=0.4), alpha=c(0.01, 0.05))
delta_p <- numeric(nrow(pow_range))
for (i in 1:nrow(pow_range)) {
  pr <- pow_range[i,]
  delta_p[i] <- effsize1d(n=pr$n, p0=pr$p1, deff=pr$Deff, alpha=pr$alpha)
}
pow_range$delta_p = delta_p
pow_range$alpha = paste0(100 - 100*pow_range$alpha, '%')
write.csv(pow_range, "data/sample_size_prop.csv", row.names=FALSE)

ggplot(pow_range) + geom_line(aes(x=p1,y=delta_p, lty=alpha)) +
  facet_grid(Deff~n) +
  theme_bw() +
  guides(lty=guide_legend('Confidence level')) +
  ylab("Uncertainty in prevalence") +
  xlab("Baseline prevalence")

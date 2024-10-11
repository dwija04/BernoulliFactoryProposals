set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")

ns <- xn[[1]]
x <- xn[[2]]
c <- xn[[3]]
t <- xn[[4]]
cov <- xn[[5]]
m <- xn[[6]]
delta_m <- xn[[7]]
N0 <- xn[[8]]
mu <- xn[[9]]


#The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

#Running MCMC using exact proposal
N <- 1e6

eta_bf <- 0.004 #step size
bf_chain <- cox_bf(N, init = rep(1,  m), ns = ns, x, c, t, cov, eta_bf)
bf_chain[[3]]

eta_mh <- 0.005
mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh)
mh_chain[[2]]

save(bf_chain, mh_chain, file = "output_cox_single_run.RData")



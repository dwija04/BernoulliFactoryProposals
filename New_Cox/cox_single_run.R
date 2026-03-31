set.seed(1)
library(mcmcse)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")

load("proposal_covariance.RData")
ns <- xn[[1]]
x <- xn[[2]]
c <- xn[[3]]
t <- xn[[4]]
cov <- xn[[5]]
m <- xn[[6]]
delta_m <- xn[[7]]
N0 <- xn[[8]]
mu <- xn[[9]]
init <- xn[[10]]
sqrt.cov <- xn[[11]]
inv.cov <- xn[[12]]
cov.svd <- xn[[13]]

N <- 1e6

#Running MCMC using exact proposal
eta_bf <- 0.06
bf_time <- system.time(bf_chain <- cox_bf(N, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_bf))
bf_chain[[3]]

# Running MCMC using RWMH
eta_rwmh <- 0.01
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = rep(1, m), ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_rwmh))
rwmh_chain[[2]]

# Running MCMC using inexact proposal
eta_mh <- 0.08
mh_time <- system.time(mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_mh))
mh_chain[[2]]
plot(density(mh_chain[[1]][, 100]))


save(bf_chain, rwmh_chain, mh_chain, file = "output_cox_single_run.RData")
save(bf_time, rwmh_time, mh_time, file = "output_cox_times.RData")



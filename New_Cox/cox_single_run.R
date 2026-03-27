set.seed(1)
library(mcmcse)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")
load("cox-proposal-matrix.RData")

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



N <- 1e6

#The posterior

cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

# Proposal matrix
prop.cov.svd <- svd(prop.cov)
prop.sqrt.cov <- prop.cov.svd$u %*% diag(prop.cov.svd$d^(1/2), m) %*% t(prop.cov.svd$v)
max(prop.cov - prop.sqrt.cov %*% prop.sqrt.cov)
prop.inv.cov <- qr.solve(prop.cov)


#Running MCMC using exact proposal
eta_bf <- 0.05
bf_time <- system.time(bf_chain <- cox_bf(N, init = init, ns = ns, x, c, t, cov = cov, prop.cov = prop.cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_bf))
bf_chain[[3]]

eta_rwmh <- 0.05
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = rep(1, m), ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_rwmh))
rwmh_chain[[2]]
plot(density(rwmh_chain[[1]][, 100]))


save(bf_chain, rwmh_chain,  file = "output_cox_single_run.RData")
save(bf_time, rwmh_time, file = "output_cox_times.RData")


set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")
# Rcpp::sourceCpp("cox_functions.cpp")

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
N <- 1e4

eta_bf <- 0.004 #step size
bf_time <- system.time(bf_chain <- cox_bf(N, init = rep(1,  m), ns = ns, x, c, t, cov, eta_bf))
bf_chain[[3]]

eta_bf <- 0.004 #step size
bf_time_new <- system.time(bf_chain_new <- cox_bf(N, init = rep(1,  m), ns = ns, x, c, t, cov, eta_bf))
bf_chain_new[[3]]

# eta_mh <- 0.005
# mh_time <- system.time(mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh))
# mh_chain[[2]]

eta_rwmh <- 0.0007
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = rep(1.5, m), ns = ns, x, c, t, cov, eta_rwmh))
rwmh_chain[[2]]

save(bf_chain, rwmh_chain, bf_chain_new, file = "output_cox_single_run.RData")
save(bf_time, rwmh_time, bf_time_new, file = "output_cox_times.RData")

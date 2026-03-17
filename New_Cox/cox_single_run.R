set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")
# load("cox-data.RData")
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

N <- 1e5
#The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

#Running MCMC using exact proposal

library(mcmcse)
eta_bf <- 0.004 #step size 

bf_time <- system.time(bf_chain <- cox_bf(N, init = rep(1,  m), ns = ns, x, c, t, cov, eta_bf))
save(bf_chain, file = "bf_chain_single.RData")
bf_chain[[3]]
mean(ess(bf_chain[[1]]))
min(ess(bf_chain[[1]]))/bf_time[3]
mean(bf_chain[[2]])
multiESS(bf_chain[[1]], r = 1)

# eta_mh <- 0.005
# mh_time <- system.time(mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh))
# mh_chain[[2]]

eta_rwmh <- 0.001
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = rep(1, m), ns = ns, x, c, t, cov, sqrt.cov, eta_rwmh))
save(rwmh_chain, file = "rwmh_chain_single.RData")

rwmh_chain[[2]]
mean(ess(rwmh_chain[[1]]))
min(ess(rwmh_chain[[1]]))/rwmh_time[3]
mean(rwmh_chain[[2]])
multiESS(rwmh_chain[[1]], r = 1)

plot(density(bf_chain[[1]][, 100]), col = "blue", main = "Posterior distribution of beta_1", xlab = "beta_1")
lines(density(rwmh_chain[[1]][, 100]), col = "red")
legend("topright", legend = c("RWMH", "Exact proposal"), col = c("blue", "red"), lty = 1)



save(bf_chain, rwmh_chain,  file = "output_cox_single_run-high-dim.RData")
save(bf_time, rwmh_time, file = "output_cox_times-high-dim.RData")

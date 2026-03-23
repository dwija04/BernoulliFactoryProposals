set.seed(1)
library(mcmcse)
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

N <- 1e4
#The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

# Initial MCMC

eta_rwmh_pilot <- 0.01
foo <- cox_rwmh(1e4, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop = prop.sqrt.cov, eta = eta_rwmh_pilot)
foo[[2]]
prop.cov <- cov(foo[[1]])

plot(density(foo[[1]][, 100]), main = "Posterior distribution of beta_100", xlab = "beta_1")
lines(density(rw_dummy[[1]][, 100]), col = "red")


prop.cov.svd <- svd(prop.cov)
prop.sqrt.cov <- prop.cov.svd$u %*% diag(prop.cov.svd$d^(1/2), m) %*% t(prop.cov.svd$v)
max(prop.cov - prop.sqrt.cov %*% prop.sqrt.cov)
prop.inv.cov <- qr.solve(prop.cov)




#Running MCMC using exact proposal
eta_bf <- 0.005
bf_time <- system.time(bf_chain <- cox_bf(1e4, init = init, ns = ns, x, c, t, cov = cov, prop.cov = prop.cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_bf))
save(bf_chain, file = "bf_chain_single.RData")
bf_chain[[3]]

plot(density(foo[[1]][, 23]), main = "Posterior distribution of beta_100", xlab = "beta_1")
lines(density(rw_dummy[[1]][, 23]), col = "red")
lines(density(bf_chain[[1]][, 23]), col = "blue")

mean(ess(bf_chain[[1]]))
min(ess(bf_chain[[1]]))/bf_time[3]
mean(bf_chain[[2]])
multiESS(bf_chain[[1]], r = 1)

# eta_mh <- 0.005
# mh_time <- system.time(mh_chain <- cox_mh(N, init = rep(1, m), ns = ns, x, c, t, cov, eta_mh))
# mh_chain[[2]]

eta_rwmh <- 0.01
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = rep(1, m), ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_rwmh))
save(rwmh_chain, file = "rwmh_chain_single.RData")

rwmh_chain[[2]]
mean(ess(rwmh_chain[[1]]))
min(ess(rwmh_chain[[1]]))/rwmh_time[3]
mean(rwmh_chain[[2]])
multiESS(rwmh_chain[[1]], r = 1)

plot(density(rwmh_chain[[1]][, 100]), col = "red", main = "Posterior distribution of beta_100", xlab = "beta_1")
lines(density(bf_chain[[1]][, 100]), col = "blue")
# plot(density(bf_chain[[1]][, 99]), col = "blue", main = "Posterior distribution of beta_100", xlab = "beta_1")
# lines(density(rwmh_chain[[1]][, 100]), col = "red")
# legend("topright", legend = c("RWMH", "Exact proposal"), col = c("blue", "red"), lty = 1)



save(bf_chain, rwmh_chain,  file = "output_cox_single_run-high-dim.RData")
save(bf_time, rwmh_time, file = "output_cox_times-high-dim.RData")

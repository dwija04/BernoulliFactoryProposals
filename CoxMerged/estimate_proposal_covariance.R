set.seed(10)
library(mcmcse)
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

# initial values
init <- lam1(t)

N <- 1e6

# The posterior covariance
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)

# Initial MCMC

eta_rwmh_pilot <- 0.001
foo <- cox_rwmh(N, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop = diag(1,m), eta = eta_rwmh_pilot)
foo[[2]]
prop.cov <- cov(foo[[1]])

prop.cov.svd <- svd(prop.cov)
prop.sqrt.cov <- prop.cov.svd$u %*% diag(prop.cov.svd$d^(1/2), m) %*% t(prop.cov.svd$v)
max(prop.cov - prop.sqrt.cov %*% prop.sqrt.cov)
prop.inv.cov <- qr.solve(prop.cov)

save(prop.cov.svd, prop.sqrt.cov, prop.inv.cov, file = "proposal_covariance.Rdata")

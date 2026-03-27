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


N <- 1e4

# Iniital MCMC
eta_rwmh_pilot <- 0.001
init <- lam1(t)

foo <- cox_rwmh(N, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop = diag(1, m), eta = eta_rwmh_pilot)
foo[[2]]
prop.cov <- cov(foo[[1]])
plot(density(foo[[1]][, 100]))

save(prop.cov, file = "cox-proposal-matrix.RData")

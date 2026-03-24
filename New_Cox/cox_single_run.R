set.seed(1)
library(mcmcse)
source("cox_functions.R")
load("estimated-cov.RData")
load("proposal_covariance.Rdata")
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

init <- lam1(t)
N <- 1e6

# The posterior
cov.svd <- svd(cov)
sqrt.cov <- cov.svd$u %*% diag(cov.svd$d^(1/2), m) %*% t(cov.svd$v)
max(cov - sqrt.cov %*% sqrt.cov)
inv.cov <- qr.solve(cov)
 

eta_rwmh_pilot <- 0.01
rwmh_time <- system.time(rwmh_chain <- cox_rwmh(N, init = init, ns = ns, x, c, t, 
                cov = cov, sqrt.prop = prop.sqrt.cov, 
                eta = eta_rwmh_pilot)
)                
rwmh_chain[[2]]


#Running MCMC using exact proposal
eta_bf <- 0.06
bf_time <- system.time(bf_chain <- cox_bf(N, init = init, ns = ns, x, c, t, cov = cov, prop.cov = prop.cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_bf))
bf_chain[[3]]

save(rwmh_chain, bf_chain, bf_time, rwmh_time, file = "cox_single.RData")

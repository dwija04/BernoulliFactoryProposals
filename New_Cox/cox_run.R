set.seed(1)
source("cox_functions.R")
load("estimated-cov.RData")
load("cox-data.RData")
load("proposal_covariance.RData")

library(foreach) 
library(doParallel)
library(mcmcse)

N <- 2e5
eta_bf <- 0.06
eta_rwmh <- 0.01

#extracting data
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


output_ram <- list()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

#Number of repetitions
reps <- 100

output_cox <- foreach(b = 1:reps) %dopar% {
  bf_time <- system.time(bf <- cox_bf(N, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_bf))
  rwmh_time <- system.time(rwmh <- cox_rwmh(N, init = init, ns = ns, x, c, t, cov = cov, sqrt.prop.cov = prop.sqrt.cov, eta = eta_rwmh))
  bf_chain <- bf[[1]]
  rwmh_chain <- rwmh[[1]]
  
  bf_loops_avg <- mean(bf[[2]], round = 2)
  bf_loops_max <- max(bf[[2]], round = 2)
  
  bf_multi_ess <- multiESS(bf_chain)
  rwmh_multi_ess <- multiESS(rwmh_chain)
  
  bf_time <- bf_time[3]
  rwmh_time <- rwmh_time[3]
  
  bf_ess <- ess(bf_chain)
  rwmh_ess <- ess(rwmh_chain)
  
  print(paste('Replication:', b))
  
  list(bf_time, rwmh_time, bf_loops_avg, bf_loops_max, bf_multi_ess, rwmh_multi_ess, bf_ess, rwmh_ess)
}

save(N, output_cox, file = "output_Cox.RData")


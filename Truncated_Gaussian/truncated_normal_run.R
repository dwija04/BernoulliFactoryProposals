source("truncated_normal_functions.R")
library(doParallel)
library(foreach)
library(mcmcse)
set.seed(1)

m <- 1e6 #Length of chain

output_trunc_gaussian <- list()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)
reps <- 100


output_trunc_gaussian <- foreach(b = 1:reps) %dopar% {
  print(b)
  bf_time <- system.time(bf <- BF_trunc_norm(n = m, var = 30, init = 1))
  mh_time <- system.time(mh <- MH_truncnorm(n = m, var = 17, init = 1))
  bf_chain <- bf[[1]]
  mh_chain <- mh[[1]]
  
  bf_loops_avg <- mean(bf[[2]])
  bf_loops_max <- max(bf[[2]])
  bf_ess <- ess(bf_chain)
  
  mh_ess <- ess(mh_chain)
  
  list(bf_time, mh_time, bf_loops_avg, bf_loops_max, bf_ess, mh_ess)
}

save(output_trunc_gaussian, file = "output_trunc_gaussian.RData")




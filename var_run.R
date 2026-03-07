set.seed(1)
library(doParallel)
source("var_functions.R")

nu <- 5
m <- 1e3
n <- 1e6
y <- 2
X <- generate_samples(nu, m*n)
X <- matrix(X, nrow = m, ncol = n)
# # foo <- var_kde(y, nu, X, m, n)
# var_IS <- foo[[1]]
# var_sn_IS <- foo[[2]]


output_var_kde <- list()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)
reps <- 1000

est_IS_n <- numeric(reps)
est_sn_IS_n <- numeric(reps)
est_IS_n_10 <- numeric(reps)
est_sn_IS_n_10 <- numeric(reps)
est_IS_n_100 <- numeric(reps)
est_sn_IS_n_100 <- numeric(reps)

output_var_kde <- foreach(b = 1:reps, .combine = 'rbind') %dopar% {
  print(b)
  x <- X[b, ]
  
  time <- system.time(estimates <- var_kde(y, nu, x, n))
  
  list(
    est_IS_n = estimates$est_IS_n,
    est_sn_IS_n = estimates$est_sn_IS_n,
    est_IS_n_10 = estimates$est_IS_n_10,
    est_sn_IS_n_10 = estimates$est_sn_IS_n_10,
    est_IS_n_100 = estimates$est_IS_n_100,
    est_sn_IS_n_100 = estimates$est_sn_IS_n_100
  )
}
save(output_var_kde, file = "output_var_kde.RData")



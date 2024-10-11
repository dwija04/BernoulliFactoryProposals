source("truncated_normal_functions.R")
set.seed(1)
m <- 1e6

bf_chain <- BF_trunc_norm(n = m, var = 30, init = 1)
mh_chain <- MH_truncnorm(n = m, var = 17, init = 1)

save(bf_chain, mh_chain, file = "output_trunc_gaussian_single_run.RData")


set.seed(1)

# Load data
load("estimated-cov.RData")
load("cox-data.RData")

# Source Rcpp MCMC functions
Rcpp::sourceCpp("cox_functions.cpp")

# Extract objects from xn
ns <- as.integer(xn[[1]])   # ensure integer vector
x <- xn[[2]]
c <- xn[[3]]
t <- xn[[4]]
cov <- xn[[5]]
m <- xn[[6]]
delta_m <- xn[[7]]
N0 <- xn[[8]]  # used in target if needed
mu <- xn[[9]]

# Number of MCMC iterations
N <- 1e6

# Step sizes
eta_bf <- 0.004
eta_mh <- 0.005
eta_rwmh <- 0.0007


bf_time <- system.time (bf_chain <- cox_bf(
  N = N,
  init = rep(1, m),
  n_r = ns,
  x = x,
  c_r = c,
  t_knots = t,
  cov_r = cov,
  eta = eta_bf,
  delta_m = delta_m
) )
bf_chain$accept_rate
# bf_chain$bernoulli_loops


# mh_time <- system.time(mh_chain <- cox_mh(
#   N = N,
#   init = rep(1, m),
#   n_r = ns,
#   x = x,
#   c_r = c,
#   t_knots = t,
#   cov_r = cov,
#   eta = eta_mh,
#   delta_m = delta_m,
#   B = 200
# ))
# mh_chain$accept_rate


rwmh_time <- system.time(rwmh_chain <- cox_rwmh(
  N = N,
  init = rep(1.5, m),
  n_r = ns,
  x = x,
  c_r = c,
  t_knots = t,
  cov_r = cov,
  eta = eta_rwmh,
  delta_m = delta_m
))
rwmh_chain$accept_rate

# -------------------------
# Save MCMC outputs
# -------------------------
save(bf_chain, rwmh_chain, bf_time, rwmh_time, file = "output_cox_single_run_temp.RData")


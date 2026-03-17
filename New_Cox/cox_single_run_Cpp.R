set.seed(1)

library(mcmcse)
# Load data
load("estimated-cov.RData")
load("cox-data.RData")

# Source Rcpp MCMC functions
Rcpp::sourceCpp("cox_functions.cpp")
source("cox_functions.R")

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
N <- 1e4

# Step sizes
eta_bf <- 0.003
# eta_mh <- 0.005
eta_rwmh <- 0.001

cov_raw <- cov(rwmh_chain$chi)
bf_time <- system.time (bf_chain <- cox_bf(
  N = N,
  init = rep(1, m),
  n_r = ns,
  x = x,
  c_r = c,
  t_knots = t,
  cov_r = cov_raw,
  eta = eta_bf,
  delta_m = delta_m
) )
bf_chain$accept_rate
mean(ess(bf_chain$chi))
mean(ess(bf_chain$chi))/bf_time[3]
mean(bf_chain$bernoulli_loops)
multiESS(bf_chain$chi, r = 1)

plot.ts(bf_chain$chi[, 91:100])
dim(bf_chain$chi)

plot(density(bf_chain$chi[, 99]), col = "blue", ylab = "Estimated Density", xlab = "x", main = "", lwd = 2)
lines(density(rwmh_chain$chi[, 99]), col = "green")
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


# rwmh_time <- system.time(rwmh_chain <- cox_rwmh( 
#   N = 1e5,
#   init = rep(1, m),
#   n_r = ns,
#   x = x,
#   c_r = c,
#   t_knots = t,
#   cov_r = cov,
#   eta = eta_rwmh,
#   delta_m = delta_m
# ))

rwmh_time <- system.time(rwmh_chain <- cox_rwmh (N, init = rep(1, m), ns, x, c, t, cov, eta))
rwmh_chain$accept_rate
mean(ess(rwmh_chain$chi))
mean(ess(rwmh_chain$chi))/rwmh_time[3]

# -------------------------
# Save MCMC outputs
# -------------------------
save(bf_chain, rwmh_chain, bf_time, rwmh_time, file = "output_cox_single_run_temp.RData")


load("output_cox_single_run_temp.RData")
# load("output_cox_times-high-dim.RData")


delta_m <- 50/(m-1)
grid <- seq(0, 50, length = 500)


est_fun1 <- numeric(length = length(grid))
est_fun2 <- numeric(length = length(grid))

bf_samps <- bf_chain[[1]]
rwmh_samps <- rwmh_chain[[1]]


ess_bf <- min(ess(bf_samps))
ess_rwmh <- min(ess(rwmh_samps))


#mESS_bf <- multiESS(bf_samps)
#mESS_rwmh <- multiESS(rwmh_samps)

time_bf <- bf_time[3]
time_rwmh <- rwmh_time[3]

ess_per_time_bf <- ess_bf/time_bf
ess_per_time_rwmh <- ess_rwmh/time_rwmh


#log posterior
log_post_bf <- bf_chain[[4]]
log_post_rwmh <- rwmh_chain[[3]]

print(paste("Min ESS BF: ", ess_bf))
print(paste("Min ESS RWMH: ", ess_rwmh))

ess_per_time_bf <- ess_bf/time_bf
ess_per_time_rwmh <- ess_rwmh/time_rwmh 

print(paste("ESS per unit time BF : ", round(ess_per_time_bf, 4)))
print(paste("ESS per unit time RWMH: ", round(ess_per_time_rwmh, 4)))

bern_loops_avg <- mean(bf_chain[[2]])
# summary(bf_chain[[2]])

print(paste("Average number of mean loops BF old bounds: ", round(bern_loops_avg, 4)))


#True density
temp <- seq(0, 50, length = 1e4)
# y_temp <- (lam1(temp))


# pdf("plots/cox-component-density.pdf")
j <- 100
plot(density(bf_samps[, j]), col = "blue", ylab = "Estimated Density", xlab = "x", main = "", lwd = 2)
lines(density(rwmh_samps[, j]), col = "green")
legend("topright", legend = c("Bernoulli factory MCMC", "RWMH"), col = c("blue", "green"), cex = 1.2, lty = 1, lwd = 2, bty = "n")

# dev.off()




library(mcmcse)
load("output_cox_single_run_temp.RData")


#PLOTS
delta_m <- 50/(m-1)
grid <- seq(0, 50, length = 100)
est_fun1 <- numeric(length = length(grid))
est_fun2 <- numeric(length = length(grid))

bf_samps <- bf_chain[[1]]
rwmh_samps <- rwmh_chain[[1]]

ESS_bf <- min(ess(bf_samps))
ESS_rwmh <- min(ess(rwmh_samps))

mESS_bf <- multiESS(bf_samps)
mESS_rwmh <- multiESS(rwmh_samps)

time_bf <- bf_time[3]
time_rwmh <- rwmh_time[3]

ess_per_time_bf <- ESS_bf/time_bf
ess_per_time_rwmh <- ESS_rwmh/time_rwmh

mESS_per_time_bf <- mESS_bf/time_bf
mESS_per_time_rwmh <- mESS_rwmh/time_rwmh

print(paste("Minimum ESS BF: ", round(ESS_bf, 1)))
print(paste("Minimum ESS RWMH: ", round(ESS_rwmh, 1)))

print(paste("Multivariate ESS BF: ", round(mESS_bf, 1)))
print(paste("Multivariate ESS RWMH: ", round(mESS_rwmh, 1)))

print(paste("ESS per second BF: ", round(ess_per_time_bf, 3)))
print(paste("ESS per second RWMH: ", round(ess_per_time_rwmh, 3)))

print(paste("Multivariate ESS per second BF: ", round(mESS_per_time_bf, 3)))
print(paste("Multivariate ESS per second RWMH: ", round(mESS_per_time_rwmh, 3)))


#log posterior
log_post_bf <- bf_chain[[4]]
log_post_rwmh <- rwmh_chain[[3]]


#True density
temp <- seq(0, 50, length = 1e4)
y_temp <- (lam1(temp))



pdf("plots/cox-component-density.pdf")
j <- 100
plot(density(bf_samps[-c(1:1000), j]), col = "blue", ylab = "Estimated Density", xlab = "x", main = "")
lines(density(rwmh_samps[-c(1:1000), j]), col = "green")
legend("topright", legend = c("Bernoulli factory MCMC", "Inexact Metropolis-Hastings", "RWMH"), col = c("blue", "red", "green"), cex = 1.2, lty = 1, lwd = 2, bty = "n")
dev.off()

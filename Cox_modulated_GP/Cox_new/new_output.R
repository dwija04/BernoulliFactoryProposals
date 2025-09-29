##########################################
# Single chain output (plots)
##########################################
library(mcmcse)
source("cox_functions.R")
load("output_cox_single_run.RData")
load("output_cox_times.RData")

#PLOTS
delta_m <- 50/(m-1)
grid <- seq(0, 50, length = 100)
est_fun1 <- numeric(length = length(grid))
est_fun2 <- numeric(length = length(grid))

bf_samps <- bf_chain[[1]]
rwmh_samps <- rwmh_chain[[1]]


ess_bf <- min(ess(bf_samps))
ess_rwmh <- min(ess(rwmh_samps))


# mESS_bf <- multiESS(bf_samps)
# mESS_rwmh <- multiESS(rwmh_samps)
 
time_bf <- bf_time[3]
time_rwmh <- rwmh_time[3]

ess_per_time_bf <- ESS_bf/time_bf
ess_per_time_rwmh <- ESS_rwmh/time_rwmh


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
y_temp <- (lam1(temp))


pdf("plots/cox-component-density.pdf")
j <- 100
plot(density(bf_samps[, j]), col = "blue", ylab = "Estimated Density", xlab = "x", main = "", lwd = 2)
lines(density(rwmh_samps[, j]), col = "green")
legend("topright", legend = c("Bernoulli factory MCMC", "RWMH"), col = c("blue", "green"), cex = 1.2, lty = 1, lwd = 2, bty = "n")
dev.off()

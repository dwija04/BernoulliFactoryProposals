load("output_trunc_gaussian.RData")

##########################################
# Replications output
##########################################
#Extracting the output
reps <- length(output_trunc_gaussian)

bf_time <- numeric(reps)
mh_time <- numeric(reps)
bf_loops_mean <- numeric(reps)
bf_loops_max <- numeric(reps)
bf_ess <- numeric(reps)
mh_ess <- numeric(reps)

for(i in 1:reps)
{
  bf_time[i] <- output_trunc_gaussian[[i]][[1]][[3]]
  mh_time[i] <- output_trunc_gaussian[[i]][[2]][[3]]
  bf_loops_mean[i] <- output_trunc_gaussian[[i]][[3]]
  bf_loops_max[i] <- output_trunc_gaussian[[i]][[4]]
  bf_ess[i] <- output_trunc_gaussian[[i]][[5]]
  mh_ess[i] <- output_trunc_gaussian[[i]][[6]]
}

#Loops
avg_bf_loops_mean <- round(mean(bf_loops_mean), 2)
avg_bf_loops_max <- round(mean(bf_loops_max), 2)

bernoulli_df <- data.frame(
  Metric = c("Average Loops", "Maximum Loops"),
  Bernoulli_factory = c(avg_bf_loops_mean, avg_bf_loops_max)
)
print(bernoulli_df)

#Effective Sample Size

avg_bf_ess <- round(mean(bf_ess))
avg_mh_ess <- round(mean(mh_ess))

avg_bf_ess_time <- round(mean(bf_ess/bf_time))
avg_mh_ess_time <- round(mean(mh_ess/mh_time))

ess_df <- data.frame(
  Metric = c("Effective Sample Size (ESS)", "ESS per unit time", "Average computing time"),
  Exact_proposal = c(avg_bf_ess, avg_bf_ess_time, round(mean(bf_time), 2)),
  Approximate_proposal = c(avg_mh_ess, avg_mh_ess_time, round(mean(mh_time), 2))
)
print(ess_df)

#Average computing times
print(mean(round(bf_time, 2)))
print(mean(round(mh_time, 2)))
##########################################
## Single chain plots
##########################################
load("output_trunc_gaussian_single_run.RData")
#Density Plots


pdf("plots/trunc_gauss_density.pdf")
x <- seq(0, 15, 0.01)
plot(density(bf_chain[[1]]), ylim = c(0, .4), col = "blue", main = " ", xlab = "x", ylab = "Estimated Density")
lines(density(mh_chain[[1]]), col = "red")
lines(x, dgamma(x, shape = 2, rate = 1), col = "black")
legend("topright", legend = c("Target density", "Bernoulli factory MCMC", "Approximate Metropolis-Hastings") ,
  col = c("black","blue","red"), lty = c(1, 1, 1), cex = 1.2, bty = "n")
dev.off()

pdf("plots/trunc_gauss_acf.pdf")
lag.max <- 40
mh_acf <- acf(mh_chain[[1]], lag.max = lag.max , plot = FALSE)$acf
bf_acf <- acf(bf_chain[[1]], lag.max = lag.max , plot = FALSE)$acf
plot(0:lag.max , mh_acf, type = 'l', col = "red",
  ylab = "Estimated Autocorrelation Function", xlab = "Lags")
lines(0:lag.max , bf_acf, col = "blue")
legend("topright", legend = c("Bernoulli factory MCMC", "Approximate Metropolis-Hastings") ,
  col = c("blue","red"), lty = c(1, 1, 1),
  cex = 1.2, bty = "n")
dev.off()



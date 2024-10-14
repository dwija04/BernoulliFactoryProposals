##########################################
## Code for reproducing the results from
## the RAM example
##########################################

load("output_RAM.RData")

reps <- length(output_ram_bern)
p <- 8 #dimension

bf_time <- numeric(reps)
aux_time <- numeric(reps)

bf_down_loops_avg <- matrix(nrow = reps, ncol = p/2)
bf_down_loops_max <- matrix(nrow = reps, ncol = p/2)
bf_up_loops_avg <- matrix(nrow = reps, ncol = p/2)
bf_up_loops_max <- matrix(nrow = reps, ncol = p/2)
bf_bern_loops_avg <- matrix(nrow = reps, ncol = p/2)
bf_bern_loops_max <- matrix(nrow = reps, ncol = p/2)

aux_down_loops_avg <- matrix(nrow = reps, ncol = p/2)
aux_down_loops_max <- matrix(nrow = reps, ncol = p/2)
aux_up_loops_avg <- matrix(nrow = reps, ncol = p/2)
aux_up_loops_max <- matrix(nrow = reps, ncol = p/2)
aux_z_loops_avg <- matrix(nrow = reps, ncol = p/2)
aux_z_loops_max <- matrix(nrow = reps, ncol = p/2)

bf_multi_ess <- numeric(reps)
aux_multi_ess <- numeric(reps)

bf_ess <- matrix(nrow = reps, ncol = p)
aux_ess <- matrix(nrow = reps, ncol = p)

for(i in 1:reps)
{
  foo1 <- output_ram_bern[[i]]
  foo2 <- output_ram_aux[[i]]
  
  bf_time[i] <- foo1[[1]]
  aux_time[i] <- foo2[[1]]
  
  bf_down_loops_avg[i, ] <- foo1[[2]]
  bf_down_loops_max[i, ] <- foo1[[3]]
  bf_up_loops_avg[i, ] <- foo1[[4]]
  bf_up_loops_max[i, ] <- foo1[[5]]
  bf_bern_loops_avg[i, ] <- foo1[[6]]
  bf_bern_loops_max[i, ] <- foo1[[7]]
  
  aux_down_loops_avg[i, ] <- foo2[[2]]
  aux_down_loops_max[i, ] <- foo2[[3]]
  aux_up_loops_avg[i, ] <- foo2[[4]]
  aux_up_loops_max[i, ] <- foo2[[5]]
  aux_z_loops_avg[i, ] <- foo2[[6]]
  aux_z_loops_max[i, ] <- foo2[[7]]
  
  bf_multi_ess[i] <- foo1[[8]]
  aux_multi_ess[i] <- foo2[[8]]
  
  bf_ess[i, ] <- foo1[[9]] 
  aux_ess[i, ] <- foo2[[9]]
  
}

#Tables for the loops

bernoulli_loops_df <- data.frame(
  Metric = c("Average Downward Loops", "Average Upward Loops", "Average Bernoulli Loops", 
             "Max Downward Loops", "Max Upward Loops", "Max Bernoulli Loops"),
  Sensor_1 = round(c(colMeans(bf_down_loops_avg)[1], colMeans(bf_up_loops_avg)[1], colMeans(bf_bern_loops_avg)[1], 
               colMeans(bf_down_loops_max)[1], colMeans(bf_up_loops_max)[1], colMeans(bf_bern_loops_max)[1]), 2),
  Sensor_2 = round(c(colMeans(bf_down_loops_avg)[2], colMeans(bf_up_loops_avg)[2], colMeans(bf_bern_loops_avg)[2],
               colMeans(bf_down_loops_max)[2], colMeans(bf_up_loops_max)[2], colMeans(bf_bern_loops_max)[2]), 2),
  Sensor_3 = round(c(colMeans(bf_down_loops_avg)[3], colMeans(bf_up_loops_avg)[3], colMeans(bf_bern_loops_avg)[3],
               colMeans(bf_down_loops_max)[3], colMeans(bf_up_loops_max)[3], colMeans(bf_bern_loops_max)[3]), 2),
  Sensor_4 = round(c(colMeans(bf_down_loops_avg)[4], colMeans(bf_up_loops_avg)[4], colMeans(bf_bern_loops_avg)[4],
               colMeans(bf_down_loops_max)[4], colMeans(bf_up_loops_max)[4], colMeans(bf_bern_loops_max)[4]), 2)
)

print(bernoulli_loops_df)


auxiliary_loops_df <- data.frame(
  Metric = c("Average Downward Loops", "Average Upward Loops", "Average Auxiliary Loops", 
             "Max Downward Loops", "Max Upward Loops", "Max Auxiliary Loops"),
  Sensor_1 = round(c(colMeans(aux_down_loops_avg)[1], colMeans(aux_up_loops_avg)[1], colMeans(aux_z_loops_avg)[1], 
               colMeans(aux_down_loops_max)[1], colMeans(aux_up_loops_max)[1], colMeans(aux_z_loops_max)[1]), 2),
  Sensor_2 = round(c(colMeans(aux_down_loops_avg)[2], colMeans(aux_up_loops_avg)[2], colMeans(aux_z_loops_avg)[2],
               colMeans(aux_down_loops_max)[2], colMeans(aux_up_loops_max)[2], colMeans(aux_z_loops_max)[2]), 2),
  Sensor_3 = round(c(colMeans(aux_down_loops_avg)[3], colMeans(aux_up_loops_avg)[3], colMeans(aux_z_loops_avg)[3],
               colMeans(aux_down_loops_max)[3], colMeans(aux_up_loops_max)[3], colMeans(aux_z_loops_max)[3]), 2),
  Sensor_4 = round(c(colMeans(aux_down_loops_avg)[4], colMeans(aux_up_loops_avg)[4], colMeans(aux_z_loops_avg)[4],
               colMeans(aux_down_loops_max)[4], colMeans(aux_up_loops_max)[4], colMeans(aux_z_loops_max)[4]), 2)
)

print(auxiliary_loops_df)


#Effective Sample Size

bf_mESS_per_unit_time <- bf_multi_ess / bf_time
aux_mESS_per_unit_time <- aux_multi_ess / aux_time

Multi_ESS_df <- data.frame(
  Method = c("Bernoulli Factory", "Auxiliary Variable"),
  MultiESS = c(round(mean(bf_multi_ess), 0), round(mean(aux_multi_ess), 0)),
  MultiESS_by_time = c(round(mean(bf_mESS_per_unit_time), 4), round(mean(aux_mESS_per_unit_time), 4)),
  Average_compute_time = c(mean(round(bf_time, 2)), mean(round(aux_time, 2)))
)

print(Multi_ESS_df)

# Compute ESS for Bernoulli factory method
ESS <- data.frame(
  Sensor = c("Sensor 1 (x)", "Sensor 1 (y)", "Sensor 2 (x)", "Sensor 2 (y)", "Sensor 3 (x)", "Sensor 3 (y)", 
             "Sensor 4 (x)", "Sensor 4 (y)"),
  Bernoulli_ESS = round(c(colMeans(bf_ess)[1], colMeans(bf_ess)[2], colMeans(bf_ess)[3], colMeans(bf_ess)[4], 
                    colMeans(bf_ess)[5], colMeans(bf_ess)[6], colMeans(bf_ess)[7], colMeans(bf_ess)[8]), 0),
  Auxiliary_ESS = round(c(colMeans(aux_ess)[1], colMeans(aux_ess)[2], colMeans(aux_ess)[3], colMeans(aux_ess)[4],
                    colMeans(aux_ess)[5], colMeans(aux_ess)[6], colMeans(aux_ess)[7], colMeans(aux_ess)[8]), 0)
)

print(ESS)


##########################################
# Plots from single run
##########################################
load("RAM-samples-single-run.RData")
set.seed(11)


#Taking a random subset of the samples for density plots
sample_indices <- sample(1:nrow(bf_chain_single[[1]]), size = 1e4, replace = FALSE)
bf_chain_single[[1]] <- bf_chain_single[[1]][sample_indices, ]
aux_chain_single[[1]] <- aux_chain_single[[1]][sample_indices, ]

sensor_1.bern <- bf_chain_single[[1]][, 1:2]
sensor_2.bern <- bf_chain_single[[1]][, 3:4]
sensor_3.bern <- bf_chain_single[[1]][, 5:6]
sensor_4.bern <- bf_chain_single[[1]][, 7:8]

pdf("plots/sensor_bernoulli.pdf")
par(mfrow = c(2, 2))
plot(sensor_1.bern[, 1], sensor_1.bern[, 2], xlab = expression(x[11]), ylab = expression(x[12]), main = "Sensor 1")
plot(sensor_2.bern[, 1], sensor_2.bern[, 2], xlab = expression(x[21]), ylab = expression(x[22]), main = "Sensor 2")
plot(sensor_3.bern[, 1], sensor_3.bern[, 2], xlab = expression(x[31]), ylab = expression(x[32]), main = "Sensor 3")
plot(sensor_4.bern[, 1], sensor_4.bern[, 2], xlab = expression(x[41]), ylab = expression(x[42]), main = "Sensor 4")
dev.off()

sensor_1.aux <- aux_chain_single[[1]][, 1:2]
sensor_2.aux <- aux_chain_single[[1]][, 3:4]
sensor_3.aux <- aux_chain_single[[1]][, 5:6]
sensor_4.aux <- aux_chain_single[[1]][, 7:8]

pdf("plots/sensor_aux.pdf")
par(mfrow = c(2, 2))
plot(sensor_1.aux[, 1], sensor_1.aux[, 2], xlab = expression(x[11]), ylab = expression(x[12]), main = "Sensor 1")
plot(sensor_2.aux[, 1], sensor_2.aux[, 2], xlab = expression(x[21]), ylab = expression(x[22]), main = "Sensor 2")
plot(sensor_3.aux[, 1], sensor_3.aux[, 2], xlab = expression(x[31]), ylab = expression(x[32]), main = "Sensor 3")
plot(sensor_4.aux[, 1], sensor_4.aux[, 2], xlab = expression(x[41]), ylab = expression(x[42]), main = "Sensor 4")
dev.off()








source("sensor_network_functions.R")
set.seed(1)
library(foreach)
library(doParallel)
library(mcmcse)

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)


# Observation indicators among the first four sensors. 
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)


# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)



# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0, 
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)

loc <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8) #initial location

#Length of chain
m <- 2e5
#variance of proposal
var = 1.08 

output_ram_bern <- list()
output_ram_aux <- list()
num_cores <- 50
doParallel::registerDoParallel(cores = num_cores)

#Number of repetitions
reps <- 100


output_ram_bern <- foreach(b = 1:reps) %dopar% {
  bf_time <- system.time( bf <- MHwG.RAM.bernoulli(initial.loc = loc, jump.scale = rep(var, 4), Ob, Os, Xb, Xs, Yb, Ys, n.sample = m, n.burn = 0))

  bf_chain <- bf[[1]]
  
  bf_down_loops_avg <- colMeans(bf$N.d)
  bf_down_loops_max <- apply(bf$N.d, 2, max)
  bf_up_loops_avg <- colMeans(bf$N.u)
  bf_up_loops_max <- apply(bf$N.u, 2, max)
  bf_bern_loops_avg <- colMeans(bf$N.bern)
  bf_bern_loops_max <- apply(bf$N.bern, 2, max)
  
  
  bf_multi_ess <- multiESS(bf_chain)
  
  bf_time <- bf_time[3]
  
  bf_ess <- ess(bf_chain)
  

  list(bf_time, bf_down_loops_avg, bf_down_loops_max, bf_up_loops_avg, bf_up_loops_max, 
       bf_bern_loops_avg, bf_bern_loops_max, bf_multi_ess, bf_ess)
}


output_ram_aux <- foreach(b = 1:reps) %dopar% {
  aux_time <- system.time(aux <- MHwG.RAM.auxiliary(initial.loc = loc, initial.aux = loc, jump.scale = rep(var, 4), Ob, Os, Xb, Xs, Yb, Ys, n.sample = m, n.burn = 0))
  
  aux_chain <- aux[[1]]
  
  
  aux_down_loops_avg <- colMeans(aux$N.d)
  aux_down_loops_max <- apply(aux$N.d, 2, max)
  aux_up_loops_avg <- colMeans(aux$N.u)
  aux_up_loops_max <- apply(aux$N.u, 2, max)
  aux_z_loops_avg <- colMeans(aux$N.z)
  aux_z_loops_max <- apply(aux$N.z, 2, max)
  
  aux_multi_ess <- multiESS(aux_chain)
  aux_time <- aux_time[3]
  aux_ess <- ess(aux_chain)
  

  list(aux_time, aux_down_loops_avg, aux_down_loops_max, aux_up_loops_avg, 
       aux_up_loops_max, aux_z_loops_avg, aux_z_loops_max, aux_multi_ess, aux_ess)
}

save(output_ram_bern, output_ram_aux, file = "output_RAM.RData")


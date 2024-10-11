source("sensor_network_functions.R")
set.seed(1)
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
m <- 1e6
scale <- 1.08
bf_chain_single <- MHwG.RAM.bernoulli(initial.loc = loc, jump.scale = rep(scale, 4), Ob, Os, Xb, Xs, Yb, Ys, n.sample = m, n.burn = 0)
aux_chain_single <- MHwG.RAM.auxiliary(initial.loc = loc, initial.aux = loc, jump.scale = rep(scale, 4), Ob, Os, Xb, Xs, Yb, Ys, n.sample = m, n.burn = 0)

save(bf_chain_single, aux_chain_single, file = "RAM-samples-single-run.RData")


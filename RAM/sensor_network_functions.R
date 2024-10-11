library(MASS)
library(mcmcse)

#To find Euclidean distance between 2 locations
norm2 <- function(loca, locb) 
{
  sqrt(sum((loca - locb)^2))
}

#evaluating the density of the target distribution
pi <- function(loc, R = 0.3, sigma = 0.02, Ob, Os, Xb, Xs, Yb, Ys) 
{
  
  First.term <- NULL
  for (i in 1 : 2) 
  {
    TEMP <- sapply(1 : 4, function(j) {
      exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Ob[j, i]) *
        (1 - exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Ob[j, i]) 
    })
    First.term <- c(First.term, TEMP)
  }
  
  
  Second.term <- NULL
  for (i in 1 : 3) 
  {
    TEMP <- sapply((i + 1) : 4, function(j) {
      exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                 loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Os[i, j]) *
        (1 - exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                        loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Os[i, j]) 
    })
    Second.term <- c(Second.term, TEMP)
  }
  
  
  First.obs.term <- NULL
  for (i in 1 : 2) 
  {
    TEMP <- sapply(1 : 4, function(j) {
      dnorm(Yb[j, i], mean = norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Ob[j, i]
    })
    First.obs.term <- c(First.obs.term, TEMP)
  }
  
  
  Second.obs.term <- NULL
  for (i in 1 : 3) 
  {
    TEMP <- sapply((i + 1) : 4, function(j) 
    {
      dnorm(Ys[i, j], mean = norm2(loc[(2 * i -1) : (2 * i)], 
                                   loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Os[i, j]
    })
    Second.obs.term <- c(Second.obs.term, TEMP)
  }
  
  
  log.lik <- sum(log(c(First.term, Second.term, First.obs.term, Second.obs.term)))
  #log likelihood + log prior
  post <- log.lik + sum(dnorm(loc, mean = rep(0, 8), sd = rep(10, 8), log = TRUE)) 
  return(post)
  
}



#RAM transition kernel using the auxiliary variable method
ram.kernel.auxiliary <- function(current.location, current.aux, loc.number, scale) 
{
  
  eps <- 10^(-308) # for stability of the acceptance ratio
  accept <- 0
  x.c <- current.location #tracks the current state
  log.x.c.den <- pi(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  z.c <- current.aux #tracks the auxiliary variable
  log.z.c.den <- pi(z.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.c.den <- exp(log.z.c.den)
  index <- (2 * loc.number - 1) : (2 * loc.number)
  
  # downhill transition
  x.p1 <- x.c
  x.p1[index] <- x.p1[index] + 
    rnorm(2, 0, scale)
  log.x.p1.den <- pi(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1  # number of downhill loops
  #Running Accept-Reject 
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) 
  {
    x.p1 <- x.c
    x.p1[index] <- x.p1[index] + 
      rnorm(2, 0, scale) #proposal
    log.x.p1.den <- pi(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }
  
  # uphill transition
  x.p2 <- x.p1 
  x.p2[index] <- x.p2[index] + 
    rnorm(2, 0, scale)
  log.x.p2.den <- pi(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  N.u <- 1 #number of uphill loops
  #Running Accept Reject 
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) 
  {
    x.p2 <- x.p1
    x.p2[index] <- x.p2[index] + 
      rnorm(2, 0, scale) #proposal
    log.x.p2.den <- pi(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }
  
  # downhill transition for the auxiliary variable z
  N.dz <- 1   #downhill loops for z
  z <- x.p2 
  z[index] <- z[index] + 
    rnorm(2, 0, scale) #proposal
  log.z.den <- pi(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.den <- exp(log.z.den)
  #Running Accept Reject for z
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) 
  {
    z <- x.p2
    z[index] <- z[index] + 
      rnorm(2, 0, scale) #proposal
    log.z.den <- pi(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    z.den <- exp(log.z.den)
    N.dz <- N.dz + 1
  }
  
  # accept or reject the proposal using MH ratio
  min.nu <- min(1, (x.c.den + eps) / (z.c.den + eps))
  min.de <- min(1, (x.p2.den + eps) / (z.den + eps))
  
  l.mh <- log.x.p2.den - log.x.c.den + log(min.nu) - log(min.de) #calculating the ratio
  
  if (l.mh > -rexp(1)) 
  {
    x.c <- x.p2 
    z.c <- z
    accept <- 1
  }
  
  c(x.c, z.c, N.d, N.u, N.dz, accept)
}



#RAM transition kernel using the Bernoulli factory method
ram.kernel.bernoulli <- function(current.location, loc.number, scale, beta = 1) 
{
  
  eps <- 10^(-308) # for stability of the acceptance ratio
  accept <- 0
  x.c <- current.location #tracks the current state
  log.x.c.den <- pi(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  
  index <- (2 * loc.number - 1) : (2 * loc.number)
  # downhill transition
  x.p1 <- x.c
  x.p1[index] <- x.p1[index] + 
    rnorm(2, 0, scale) 
  log.x.p1.den <- pi(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1 #number of downhill loops
  
  #Running Accept-Reject
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) 
  {
    x.p1 <- x.c
    x.p1[index] <- x.p1[index] + 
      rnorm(2, 0, scale) #proposal
    log.x.p1.den <- pi(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }
  
  # uphill transition
  x.p2 <- x.p1
  x.p2[index] <- x.p2[index] + 
    rnorm(2, 0, scale)
  log.x.p2.den <- pi(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  N.u <- 1 #number of uphill loops  
  
  #Running Accept-Reject
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) 
  {
    x.p2 <- x.p1
    x.p2[index] <- x.p2[index] + rnorm(2, 0, scale) #proposal
    log.x.p2.den <- pi(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }
  y <- x.p2 #The proposal obtained using RAM 
  
  #The bounds for using Bernoulli factories
  cy <-  x.p2.den #exp(pi(y, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys))
  cx <- x.c.den #exp(pi(current.location, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys))
  C <- cy / (cy + cx)
  
  accept <- 0
  N.bern <- 0 #number of Bernoulli loops
  
  current.x <- current.location
  
  #Running Bernoulli factory
  while(!accept)
  {
    N.bern <- N.bern + 1
    C1 <- rbinom(1, 1,  C )
    if(C1 == 1)
    { 
      m1 <- current.x
      m1[index] <- m1[index] + rnorm(2, 0, scale) 
      p_y <- min(1, (cx + eps) /  ( exp(pi(m1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps) )
      C2 <- rbinom(1, 1, p_y)
      if(C2 == 1)
      {
        ret <- 1
        accept <- 1
      }
    }
    else
    {
      m2 <- y
      m2[index] <- m2[index] +  rnorm(2, 0, scale) 
      p_x <- min(1, ( cy + eps)/ ( exp(pi(m2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps) )
      C2 <- rbinom(1, 1, p_x)
      
      if(C2 == 1)
      {
        ret <- 0
        accept <- 1
      }
    }
  }
  
  if(ret == 1)
  {
    x.c <- y 
    accept <- 1
  }
  else
  {
    x.c <- current.x
    accept <- 0
  }
  
  c(x.c, N.d, N.u, N.bern, accept)
}


#Running the MCMC chain using RAM kernel with auxiliary variables
MHwG.RAM.auxiliary <- function(initial.loc, initial.aux, jump.scale, 
                     Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {
  
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 4) 
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  Nz <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) 
  {
    for (j in 1 : 4) 
    {
      TEMP <- ram.kernel.auxiliary(loc.t, aux.t, j, jump.scale[j])
      #Storing the results
      loc.t <- TEMP[1 : 8] #storing the values of Xs
      aux.t <- TEMP[9 : 16] #storing the auxiliary variables
      Nd[i, j] <- TEMP[17] #downhill loops
      Nu[i, j] <- TEMP[18] #uphill loops
      Nz[i, j] <- TEMP[19] #downhill loops for z
      accept[i, j] <- TEMP[20] #1 if proposal is accepted, and 0 otherwise
    }
    out[i, ] <- loc.t
    if(i%%1e4 == 0) print(i)
  }
  list(x = out[-c(1 : n.burn), ], 
       accept = accept[-c(1 : n.burn), ],
       N.d = Nd[-c(1 : n.burn), ],
       N.u = Nu[-c(1 : n.burn), ],
       N.z = Nz[-c(1 : n.burn), ])
  
}


#Running the MCMC chain using RAM kernel with Bernoulli factories1`  q`
MHwG.RAM.bernoulli <- function(initial.loc, jump.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample, n.burn = 0) {
  
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  Nd <- matrix(NA, nrow = n.total, ncol = 4)
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  N_bern <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) 
  {
    for (j in 1 : 4) 
    {
      TEMP <- ram.kernel.bernoulli(loc.t, j, jump.scale[j])
      #Storing the results
      loc.t <- TEMP[1 : 8] #storing the values of Xs
      Nd[i, j] <- TEMP[9] #downhill loops
      Nu[i, j] <- TEMP[10] #uphill loops
      N_bern[i, j] <- TEMP[11] #Bernoulli loops
      accept[i, j] <- TEMP[12] #1 if proposal is accepted, and 0 otherwise
    }
    out[i, ] <- loc.t
    if(i%%1e4 == 0) print(i)
  }
  list(x = out[-c(1 : n.burn), ], 
       accept = accept[-c(1 : n.burn), ],
       N.d = Nd[-c(1 : n.burn), ],
       N.u = Nu[-c(1 : n.burn), ],
       N.bern = N_bern[-c(1 : n.burn), ])
  
}

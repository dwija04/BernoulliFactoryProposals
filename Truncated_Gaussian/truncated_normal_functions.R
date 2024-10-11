library(truncnorm)
library(mcmcse)


#Gamma (2, 1) target distribution
pi <- function(x, alpha, beta) 
{
  ret <- x^(alpha-1) * exp(-beta*x)
  return(ret)
}

#function to calculate the bounds
C_f <- function(x, alpha, beta)
{
  ret <- x^(alpha - 1) * exp((-1)*beta*x)
  return(ret)
}

#Running the MCMC chain with exact proposals using Bernoulli factories
BF_trunc_norm <- function(n, alpha = 2, beta = 1, var, init)
{
  x <- numeric(n)
  x[1] <- init 
  loops_total <- numeric(n)
  accept_rate <- 0 
  
  for(i in 2:n)
  {
    y <- rtruncnorm(1, a = 0, b = Inf, mean = x[i-1], sd = sqrt(var))
    
    #calculating the bounds
    c_num <- C_f(y, alpha, beta)
    c_den <- C_f(x[i-1], alpha, beta)
    
    C <- c_num / (c_den + c_num)
    
    accept <-0 
    loops <- 0
    #running the Bernoulli factory
    while(!accept)
    {
      loops <- loops + 1
      C1 <- rbinom(1, 1, C)
      if(C1 == 1)
      {
        M <- rnorm(1, mean = x[i-1], sd = sqrt(var))

        p_num <- (M > 0) 
        if(p_num == 1)
        {
          x[i] <- y
          accept <- 1
          accept_rate <- accept_rate + 1
        }
      }
      else
      {
        M <- rnorm(1, mean = y, sd = sqrt(var))
        p_den <- (M > 0) 
        if(p_den == 1)
        {
          x[i] <- x[i-1]
          accept <- 1
        }
      }
    }
    loops_total[i] <- loops
  }
  accept_rate <- accept_rate/n
  return(list(x, loops_total, accept_rate))
}

#Running MCMC chain with approximate proposals using Metropolis Hastings
MH_truncnorm <- function(n, alpha = 2, beta = 1,var, init) 
{
  x <- numeric(n)
  x[1] <- init
  accept <- 0
  for(i in 2:n) 
  {
    y <- rtruncnorm(1, a = 0, b = Inf, mean = x[i-1], sd = sqrt(var))
    #Evaluating MH ratio
    ratio <- (dgamma(y, shape = alpha, rate = beta) * dtruncnorm(x[i-1], a = 0, b = Inf, mean = y, sd = sqrt(var))) /
      (dgamma(x[i-1], shape = alpha, rate = beta) * dtruncnorm(y, a = 0, b = Inf, mean = x[i-1], sd = sqrt(var)))
    acceptance_prob <- min(1, ratio)
    
    if(runif(1) < acceptance_prob) {
      x[i] <- y
      accept <- accept + 1
    } else {
      x[i] <- x[i-1]
    }
  }
  accept_rate <- accept/n
  return(list(x, accept_rate))
}

#install.packages(c("nloptr", "quadprog", "mvtnorm", "TruncatedNormal", "plot3D"))
# library("lineqGPR")
library(MASS)
library(TruncatedNormal)

lam1 <- function(x)
{
  return( 2*exp(-x/15) + exp(-((x-25)/10)^2) )
}

phi <- function(del, x, t)
{
  val <- abs((x-t)/del)
  val <- (1-val) * (val <= 1) + 0 *(val > 1)
  return(val)
}

target <- function(chi, x, c, t, cov, n)
{
  n0 <- length(n)
  mat <- -0.5 %*% t(chi)%*%inv.cov%*%chi - n0*t(c)%*%chi
  ret <- mat
  
  for(j in 1:N0)
  {
    track <- 0
    temp <- x[[j]]
    for(i in 1:n[j])
    {
      track <- track + log(t(phi(delta_m, temp[i], t)) %*% chi)
    }
    ret <- ret + track
  }
  return(ret)
}

is_positive <- function(x)
{
  if(!prod(x > 0)) return(0)
  return(1)
}

#Running MCMC using exact proposal
cox_bf <- function(N, init, ns, x, c, t, cov, eta)
{
  p <- length(init)
  chi <- matrix(0, N, p)
  log.post <- numeric(length = N)
  chi[1,] <- init
  accept_rate <- 0
  log.post[1] <- target(chi[1, ], x, c, t, cov, ns)
  bernoulli_loops <- numeric(N)
  for(i in 2:N)
  {
    accept <- 0
    if(i%%(N/10) == 0) print(i)
    
    #Accept-Reject step
    while(!accept)
    {
      z <- chi[i-1, ] + as.numeric(sqrt.cov %*%rnorm(p, sd = sqrt(eta)))
      if(is_positive(z)) 
      {
        y <- z
        accept <- 1
      }
    }
    #Bernoulli factory
    c1 <- target(y, x, c, t, cov, ns)
    c2 <- target(chi[i-1, ], x, c, t, cov, ns)
    C <- exp(c1 - c2)/(1 + exp(c1 - c2))
    bern_loops <- 0
    accept <- FALSE
    while(!accept) 
    {
      bern_loops <- bern_loops + 1
      C1 <- rbinom(1, 1, C)
      if(C1 == 1) 
      {
        m1 <- chi[i-1, ] + as.numeric(sqrt.cov %*%rnorm(m, sd = sqrt(eta)))
        if(is_positive(m1)) p_x <- 1
        else p_x <- 0
        C2 <- rbinom(1, 1, p_x)
        if(C2 == 1) 
        {
          chi[i, ] <- y
          log.post[i] <- c1
          accept <- TRUE
        }
      } 
      else 
      {
        m2 <- y + as.numeric(sqrt.cov %*%rnorm(m, sd = sqrt(eta)))
        if(is_positive(m2)) p_y <- 1
        else p_y <- 0
        C2 <- rbinom(1, 1, p_y)
        if(C2 == 1) 
        {
          chi[i, ] <- chi[i-1, ]
          log.post[i] <- c2
          accept <- TRUE
        }
      }
    }
    bernoulli_loops[i] <- bern_loops
    if( chi[i, 1] == y[1] ) accept_rate <- accept_rate + 1
  }
  bern_loops <- bern_loops/N
  accept_rate <- accept_rate/N
  return(list(chi, bernoulli_loops, accept_rate, log.post))
}

#Running MCMC using approximate proposal
cox_mh <- function(N, init, ns, x, c, t, cov, eta)
{
  p <- length(init)
  chi <- matrix(0, N, p)
  log.post <- numeric(length = N)
  chi[1,] <- init
  accept_rate <- 0
  log.post[1] <- target(chi[1, ], x, c, t, cov, ns)
  
  for(i in 2:N)
  {
    accept <- 0
    if(i%%(N/10) == 0) print(i)
    
    #Accept-Reject step
    while(!accept)
    {
      z <- chi[i-1, ] + as.numeric(sqrt.cov %*%rnorm(p, sd = sqrt(eta)))
      if(is_positive(z)) 
      {
        y <- z
        accept <- 1
      }
    }
    #Acceptance ratio
    num <- target(y, x, c, t, cov, ns) 
    denom <- target(chi[i-1, ], x, c, t, cov, ns)
    log.alpha.k <- num - denom
    log.beta.k <- pmvnorm(mu = chi[i-1, ], sigma = cov, lb = rep(0, p), B = 200, type = "mc") 
    log.beta.k <- log.beta.k - pmvnorm(mu = y, sigma = cov, lb = rep(0, p), B = 200, type = "mc")      
    log.ratio <- log.alpha.k + log.beta.k
    if(runif(1) < min(1, exp(log.ratio)))
    {
      chi[i, ] <- y
      log.post[i] <- num
    }
    else
    {
      chi[i, ] <- chi[i-1, ]
      log.post[i] <- denom
    }
    if( chi[i, 1] == y[1] ) accept_rate <- accept_rate + 1
  }
  accept_rate <- accept_rate/N
  return(list(chi, accept_rate, log.post))
}

# To smooth the output
smooth <- function(delta_m, t, samp, val)
{
  t(phi(delta_m, val, t)) %*% samp
}



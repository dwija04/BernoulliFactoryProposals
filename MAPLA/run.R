### Generating data
set.seed(123)

n <- 1000
d <- 500
p <- 0.3

theta <- numeric(d)
for(j in 1:d)
{
  if(runif(1) < p)
    theta[j] <- rnorm(1, mean = 0, sd = 1)
  else
    theta[j] <- 1-p
}

X <- matrix(rnorm(n*d), nrow = n, ncol = d)
Y <- X %*% theta + rnorm(n)

### Functions

f.grad <- function(theta, X, Y)
{
  grad <- t(X) %*% (X %*% theta - Y) 
  return(as.vector(grad))
}

g.func <- function(lam, theta) return(lam * sum(abs(theta)))
f.func <- function(theta, X, Y) return(sum((X %*% theta - Y)^2))

log.pi.func <- function(theta, X, Y, lam)
{
  return(-f.func(theta, X, Y) - g.func(lam, theta))
}



oracle <- function(eta, g, X, Y, theta, lam)
{
  grad <- f.grad(theta, X, Y)
  theta_new <- theta - eta * grad
  theta_new <- sign(theta_new) * pmax(0, abs(theta_new) - eta * lam)
  return(theta_new)
}


MAPLA <- function(n.samps, theta.init, eta)
{
  samps <- matrix(0, nrow = n.samps, ncol = d)
  samps[1, ] <- theta.init
  for(i in 2:n)
  {
    x <- samps[i-1, ]
    z <- oracle(eta, g.func, X, Y, samps[i-1, ], lam)
    log.temp <- f.func(x, X, Y) - f.func(z, X, Y) 
                - 0.25*(1/eta) * sum((z - x + eta * f.grad(z, X, Y))^2)
                + 0.25*(1/eta) * sum((z - x + eta * f.grad(x, X, Y))^2)
    
  }
}
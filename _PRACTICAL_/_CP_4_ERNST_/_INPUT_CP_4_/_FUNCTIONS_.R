#-------------------------------------------------------------------------------

if (!require("sde", quietly = TRUE)) {
  install.packages("sde")
  library("sde")
} else {
  library("sde")
}

#-------------------------------------------------------------------------------
# 1. Function aimed at reproducing the asymptotic behaviour of MR process ####
BB.single <- function(dim.k, n.sim=5000, n=2000){
  Bridge.data <- replicate(n.sim, 
                           {
                             BB.data <- replicate(dim.k, BBridge(0,0,N=n-1))
                             process <- apply(BB.data, 1, function(x) crossprod(x,x))
                             return(process)
                           })
  BB.stat <- apply(Bridge.data, 2, function(x) max(abs(x)))
  return(list(Bridge.data, BB.stat))
}

#-------------------------------------------------------------------------------
# 2. Approximating Kolmogorov Distribution ####
f <- function(x,i) {exp(-(2*i-1)^2*pi^2/(8*x^2)) }
kolm <- function(x) {sqrt(2*pi)/x*(f(x,1)+f(x,2)+f(x,3)+f(x,4)+
                                     f(x,5)+f(x,6)+f(x,7)+f(x,8)+
                                     f(x,9)+f(x,10))}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. GOF ####

## 3.1. GOF - Univariate ####
GOF_univariate <- function(gam.fit, index){
  
  # Extract the design matrix from the GAM fit
  X <- model.matrix(gam.fit)
  n.e <- nrow(X)
  
  # Calculate the mean response values
  mu <- gam.fit$fitted.values
  
  # Compute the working residuals
  w <- gam.fit$family$mu.eta(gam.fit$linear.predictors)*
    (gam.fit$y - mu)/(gam.fit$sig2*gam.fit$family$variance(mu))
  
  # Psi process 
  Psi <- X[,index] * w
  W <- cumsum(Psi)/sqrt(n.e)
  # Inverse Variance of Psi
  Vs <- ginv(crossprod(Psi)/(n.e))
  # Squared Inverse Variance of Psi
  E <- eigen(Vs)
  P <- E$vectors
  D.val <- pmax(E$value, 0)
  D.sq <- (sqrt(D.val))
  B.inv.sq <- (P%*%D.sq%*%t(P))
  
  # Standardized Martingale-Residual Process
  efp <- t(B.inv.sq %*% t(W))
  
  # Compute the p-value
  pvalue = 1 - kolm(max(abs(efp)))
  
  return(list(pvalue, efp))
}

## 3.2. GOF - Multivariate ####
GOF_multivariate <- function(gam.fit, index, BB.stat){
  
  # Extract the design matrix from the GAM fit
  X <- model.matrix(gam.fit)
  n.e <- nrow(X)
  
  # Calculate the mean response values
  mu <- gam.fit$fitted.values
  
  # Compute the working residuals
  w <- gam.fit$family$mu.eta(gam.fit$linear.predictors)*
    (gam.fit$y - mu)/(gam.fit$sig2*gam.fit$family$variance(mu))
  
  # Compute the weights for the process
  wt <- sqrt(mu*(1-mu))
  
  # Determine the dimension of the covariate set
  dim.k = length(index)
  
  # Extract the columns corresponding to the covariate set from the design
  dim.spl = index
  
  # Psi for multivariate covariates is a matrix
  Psi <- X[,dim.spl] * w * wt
  Psi.last <- matrix(apply(Psi, 2, sum), ncol=dim.k)
  Psi.tilde <- Psi - apply(Psi.last, 2, 
                           function(x) x * rep(1/n.e, n.e))
  W <- apply(Psi.tilde, 2, cumsum)/sqrt(n.e)
  # Inverse Variance of Psi
  Vs <- solve(crossprod(Psi.tilde)/(n.e))
  # Squared Inverse Variance of Psi
  E <- eigen(Vs)
  P <- E$vectors
  D.val <- pmax(E$value, 0)
  D.sq <- diag(sqrt(D.val))
  B.inv.sq <- (P%*%D.sq%*%t(P))
  
  # Standardized Martingale-Residual Process
  efp <- t(B.inv.sq %*% t(W))
  
  # Compute the maximum statistics
  U.n <- apply(efp, 1, function(x) crossprod(x,x))
  stat <- max(abs(U.n))
  
  # Compute the p-value
  pvalue = mean(BB.stat>=stat)
  
  return(list(pvalue, efp))
}
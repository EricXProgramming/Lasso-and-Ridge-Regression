### Load some required packages and extra functions. 
library('MASS') # Needed for the mvrnorm function. 

### Function to implement Algorithm 1: the Zhang and Oles Cyclic Coordinate ### Descent (CLG) optimization algorithm using
### a particular choice of trust region update. 

## Function Inputs: 
# X is a n x d design matrix for our model.
# y is a n x 1 vector of our observed response variables. 
# tau is a d x 1 vector consisting of the variances from the 
#     univariate Gaussian priors with mean 0 on each of the d parameters. 
# e is tolerance condition to stop the algorithm once the desired 
#   amount of error is achieved, defaulted to 0.01.  

## Function Outputs:
# A d x 1 vector of the estimated MLE parameters for our model.

CLG <- function(tau, X, y, e = 0.01){
  n <- dim(X)[1] # The number of observations for the data.
  d <- dim(X)[2] # The number of parameters for your model. 
  
  # Initialize some values for our algorithm. 
  Beta <- rep(0, times = d)
  Delta <- rep(1, times = d)
  r <- rep(0, times = n)
  Convergence <- FALSE
  
  while(Convergence == FALSE){
    Beta.old <- Beta # Save the previous values of your estimated parameters. 
    for(j in 1:d){
      # Compute the tentative step. 
      Delta.v.j <- TentativeStep(j, Beta[j], tau[j], Delta[j], r, X, y)
      
      # Update r as well as Beta.
      Delta.Beta.j <- min(max(Delta.v.j, -Delta[j]), Delta[j]) # Limit step to trust region. 
      Delta.r <- Delta.Beta.j * X[, j] * y
      r <- r + Delta.r
      Beta[j] <- Beta[j] + Delta.Beta.j
      
      # Update the size of the trust region. 
      Delta[j] <- max(2 * abs(Delta.Beta.j), Delta[j] / 2)
    }
    
    # Checks to see if the stopping condition is achieved. 
    if(norm(Beta - Beta.old, type = '2') < e){
      Convergence = TRUE
    }
    
  }
  return(Beta)
}

### Function to compute the jth (corresponding to the jth parameter) 
### tentative step [equation 18] for a Gaussian prior. 

## Function Inputs:
# j represents the jth tentative step for the jth parameter. 
# Beta.j is the jth parameter. 
# tau.j is the variance from a univariate Gaussian prior with mean 0 on the jth parameter. 
# Delta.j is the jth step size. 
# X is a n x d Design matrix for our model.
# y is a n x 1 Response variable. 

## Function Outputs:
# The jth tentative step. 

TentativeStep <- function(j, Beta.j, tau.j, Delta.j, r, X, y){
  # Compute F for each observation which is "essentially" 
  # an upper bound on the second derivative of f.
  delta <- Delta.j * abs(X[, j])
  F <- 1 / (2 + exp(abs(r) - delta) + exp(delta - abs(r)))
  F[abs(r) <= delta] <- .25
  
  # Compute the jth tentative step. 
  top <- sum((X[, j] * y) / (1 + exp(r))) - (Beta.j / tau.j)
  bottom <- (t(X[, j] ^ 2) %*% as.matrix(F)) + (1 / tau.j)
  Delta.v.j <- as.numeric(top / bottom) # Convert data type to a number. 
  return(Delta.v.j)
}

### Function to implement the Zhang and Oles Cyclic Coordinate 
### Descent-lasso (CLG-lasso) optimization algorithm using
### a particular choice of trust region update. 

## Function Inputs: 
# X is a n x d design matrix for our model.
# y is a n x 1 vector of our observed response variables. 
# lambda is a d x 1 is the scale parameter from a Laplace prior with mean 0 
#        (center at 0) on the jth parameter.
# e is tolerance condition to stop the algorithm once the desired 
#   amount of error is achieved, defaulted to 0.01.  

## Function Outputs:
# A d x 1 vector of the estimated MLE parameters for our model.

CLG.lasso <- function(lambda, X, y, e = 0.01){
  n <- dim(X)[1] # The number of observations for the data.
  d <- dim(X)[2] # The number of parameters for your model. 
  
  # Initialize some values for our algorithm. 
  Beta <- rep(0, times = d)
  Delta <- rep(1, times = d)
  r <- rep(0, times = n)
  Convergence <- FALSE
  
  while(Convergence == FALSE){
    Beta.old <- Beta # Save the previous values of your estimated parameters. 
    for(j in 1:d){
      # Compute the tentative step. 
      Delta.v.j <- TentativeStep.Lasso(j, Beta[j], lambda[j],
                                       Delta[j], r, X, y)
      
      # Update r as well as Beta.
      Delta.Beta.j <- min(max(Delta.v.j, -Delta[j]), Delta[j]) # Limit step to trust region. 
      Delta.r <- Delta.Beta.j * X[, j] * y
      r <- r + Delta.r
      Beta[j] <- Beta[j] + Delta.Beta.j
      
      # Update the size of the trust region. 
      Delta[j] <- max(2 * abs(Delta.Beta.j), Delta[j] / 2)
    }
    
    # Checks to see if the stopping condition is achieved. 
    if(norm(Beta - Beta.old, type = '2') < e){
      Convergence = TRUE
    }
    
  }
  return(Beta)
}

### Function to implement Algorithm 2: Computation of the tentative step
### in Cyclic Coordinate Descent-lasso (CLG-lasso). 

## Function Inputs:
# j represents the jth tentative step for the jth parameter. 
# Beta.j is the jth parameter. 
# lambda.j is the scale parameter from a Laplace prior with mean 0 
#          (center at 0) on the jth parameter. 
# Delta.j is the jth step size. 
# X is a n x d Design matrix for our model.
# y is a n x 1 Response variable. 

## Function Outputs:
# The jth tentative step used for CLG-lasso. 

TentativeStep.Lasso <- function(j, Beta.j, lambda.j, Delta.j, r, X, y){
  if(Beta.j == 0){
    s <- 1 # Try positive direction. 
    
    # Compute the tentative step (Delta.v.j.) using formula 20.
    delta <- Delta.j * abs(X[, j])
    F <- 1 / (2 + exp(abs(r) - delta) + exp(delta - abs(r)))
    F[abs(r) <= delta] <- .25
    top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda.j * sign(Beta.j))
    bottom <- (t(X[, j] ^ 2) %*% as.matrix(F))
    Delta.v.j <- as.numeric(top / bottom) 
    
    if(Delta.v.j <= 0){ # Positive direction failed.
      s <- -1 # Try negative direction.
      
      # Compute the tentative step (Delta.v.j.) using formula 20.
      delta <- Delta.j * abs(X[, j])
      F <- 1 / (2 + exp(abs(r) - delta) + exp(delta - abs(r)))
      F[abs(r) <= delta] <- .25
      top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda.j * sign(Beta.j))
      bottom <- (t(X[, j] ^ 2) %*% as.matrix(F))
      Delta.v.j <- as.numeric(top / bottom)
      
      if(Delta.v.j >= 0){ # Negative direction failed.
        Delta.v.j <- 0
      }
    }
    
  } else {
    
    s <- Beta.j / abs(Beta.j)
    
    # Compute the tentative step (Delta.v.j.) using formula 20.
    delta <- Delta.j * abs(X[, j])
    F <- 1 / (2 + exp(abs(r) - delta) + exp(delta - abs(r)))
    F[abs(r) <= delta] <- .25
    top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda.j * sign(Beta.j))
    bottom <- (t(X[, j] ^ 2) %*% as.matrix(F))
    Delta.v.j <- as.numeric(top / bottom)
    
    if ((s * (Beta.j + Delta.v.j)) < 0){ # Cross over 0.
      Delta.v.j <- -Beta.j
    }
  } 
  return(Delta.v.j)
} 

### Function to implement the Hosmer-Lemeshow Goodness-of-Fit 
### (GOF) Test for binary response data. 

## Function Inputs: 
# y is the true observed response values.
# phat is the estimated response values from your logistic regression
#      model. 
# ngrp is the number of groups to consider, defaulted to 10. 
# print.table is an indicator for whether or not you want
#             the fitted probability table printed out,
#             defaulted to TRUE. 

## Function Outputs:
# Fitted probability table. 
# Hosmer-Lemeshow GOF Test Chi-squared test statistic, 
# degrees of freedom, and p-value. 

## Notes: 
# A higher chi-squared statistic gives you a lower p-value. 

binary.gof <- function(y, phat, ngrp = 10, print.table = FALSE){
  # y <- fit$y
  # phat <- fitted(fit)
  fittedgrps <- cut( phat, quantile( phat, seq(0,1,by=1/ngrp) ), include.lowest=TRUE )
  n <- aggregate( y, list( fittedgrps ), FUN=length )[,2]
  Obs <- aggregate( y, list( fittedgrps ), FUN=sum )[,2]
  Exp <- aggregate( phat, list( fittedgrps ), FUN=sum )[,2]
  if( print.table==TRUE ){
    cat( "\nFitted Probability Table:\n\n" )
    rslt <- as.data.frame( cbind( 1:ngrp, n, Obs, Exp ) )
    names( rslt )[1] <- "group"
    print( rslt )
  }
  chisqstat <- sum( (Obs - Exp)^2 / ( Exp*(1-Exp/n) ) )
  df <- ngrp-2
  pVal <- pchisq( chisqstat, df, lower.tail=FALSE )
  # cat( "\n Hosmer-Lemeshow GOF Test:\n\n" )
  # cbind(chisqstat, df, pVal)
  results <- cbind(chisqstat, df, pVal)
  return(results)
}

### Function to generate logistic regression simulated data where 
### explanatory variables are sampled from a multivariate
### normal distribution. 

## Function Inputs: 
# n.obs is the number of observations you would like your
#       dataset to have (i.e. the number of rows of your
#       design matrix X). 
# betas is a vector of the true logistic regression parameter values #       for your simulated data set. 

## Function Outputs:
# Simulated data set given as a data frame object. 

## Notes: 
# If your explanatory variables have a particular distribution, then
#    you want to consider those situations. 
# We set the covariance structure for the multivariate random
#    normal variable draws to be the identity matrix in order to
#    obtain independence between each of the covariates.
# Set some of your true parameter values to be equal to 0 in
#     order to guarantee that not all of your responses will be 1.
# Try to keep your true logistic regression model parameter values 
#     to be somewhat small. 
simulateBinaryData <- function(n.obs, betas){
  # Initialize some values. 
  n.params <- length(betas)
  mu <- rep(0, times = n.params)
  Sigma <- diag(x = 1, nrow = n.params, ncol = n.params)
  
  # Generate your explanatory variables using a 
  # multivariate random normal distribution assuming that
  # covariates are independent from one another. 
  X <- mvrnorm(n = n.obs, mu, Sigma)
  
  # Generate your binary response variable.
  z <- X %*% betas
  p <- 1 / (1 + exp(-z))
  Y  <- rbinom(n = n.obs, size = 1, prob = p)
  
  # Return your simulated data set.
  data <- as.data.frame(cbind(X, Y))
  return(data)
}


set.seed(123) # Guarantee consistent simulation results.
data1 <- simulateBinaryData(n.obs = 100000,
                            betas = c(1, 0.3, 0.4, 0, 0))
mod1 <- glm(Y ~ . - 1, family = "binomial", data = data1)
summary(mod1)
binary.gof(mod1$y, fitted(mod1))
CLG(tau = c(.01, .01, .01, .01, .01),
    X = data1[, -dim(data1)[2]],
    y = data1[, dim(data1)[2]]  , e = 0.01)
CLG.lasso(lambda = c(.01, .01, .01, .01, .01),
          X = data1[, -dim(data1)[2]],
          y = data1[, dim(data1)[2]]  , e = 0.01)


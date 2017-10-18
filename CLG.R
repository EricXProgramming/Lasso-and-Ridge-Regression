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
  Delta <- rep(0, times = d)
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
  top <- sum((X[, j] * y) / (1 + exp(r))) - (Beta[j] / tau[j])
  bottom <- (t(X[, j] ^ 2) %*% as.matrix(F)) + (1 / tau[j])
  Delta.v.j <- as.numeric(top / bottom) # Convert data type to a number. 
  return(Delta.v.j)
  }

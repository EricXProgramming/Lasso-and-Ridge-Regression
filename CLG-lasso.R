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
  Delta <- rep(0, times = d)
  r <- rep(0, times = n)
  Convergence <- FALSE
  
  while(Convergence == FALSE){
    Beta.old <- Beta # Save the previous values of your estimated parameters. 
    for(j in 1:d){
      # Compute the tentative step. 
      Delta.v.j <- TentativeStep.lasso(j, Beta[j], lambda[j], Delta[j], r, X, y)
      
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
    top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda[j] * sign(Beta[j]))
    bottom <- (t(X[, j] ^ 2) %*% as.matrix(F))
    Delta.v.j <- as.numeric(top / bottom) 
    
    if(Delta.v.j <= 0){ # Positive direction failed.
      s <- -1 # Try negative direction.
      
      # Compute the tentative step (Delta.v.j.) using formula 20.
      delta <- Delta.j * abs(X[, j])
      F <- 1 / (2 + exp(abs(r) - delta) + exp(delta - abs(r)))
      F[abs(r) <= delta] <- .25
      top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda[j] * sign(Beta[j]))
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
      top <- sum((X[, j] * y) / (1 + exp(r))) - (lambda[j] * sign(Beta[j]))
      bottom <- (t(X[, j] ^ 2) %*% as.matrix(F))
      Delta.v.j <- as.numeric(top / bottom)
      
      if ((s * (Beta.j + Delta.v.j)) < 0){ # Cross over 0.
        Delta.v.j <- -Beta.j
      }
    } 
  return(Delta.v.j)
} 

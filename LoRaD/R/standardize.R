#' This function takes a data frame y, that contains the log transform parameter values
#' The log kernel is always the last column 
#'
#' @y data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @return A new data frame consisting of standardize parameter values
#' @export 
#'
standardize <- function(y) {
  cat("\nStandardizing Params:\n")
  
  #Getting rid of log kernel because it's not a parameter
  last_col_num <- ncol(y)
  x <- as.matrix(y[,-last_col_num])
  
  #Getting dimensions of matrix
  n <- nrow(x)
  p <- ncol(x)
  
  rmax <- 0.0
  #Getting rmax
  for (i in 1:n) {
    #Making Vector of parameters
    v <- x[i,]
    #Squaring elements
    r <- sqrt(sum(v^2))
    if (r > rmax){
      rmax <- r
    }
  }
  
  #Calculating the mean vector
  m <- matrix(rep(colMeans(x), n), byrow = TRUE, ncol = p)
  
  #Calculating the Variance-Covariance Matrix
  s <-cov(x)
  
  #Getting Eigen vectors and values
  lamb <- eigen(s)$value
  vec <- eigen(s)$vector
  
  #Calculating sqrt of s
  sqrts <- vec%*%diag(sqrt(lamb))%*%t(vec)
  invsqrts <- vec%*%diag(1/sqrt(lamb))%*%t(vec)
  
  #Performing the standardization
  z <- (x-m)%*%invsqrts
  
  #Computing log Jacobian
  logJ <- as.numeric(determinant(invsqrts,logarithm = TRUE)$modulus)
  
  #Converting z to dataframe
  zdf <- as.data.frame(z)
  
  #Adding log Jacobian to log kernel to complete transformation
  logkern <- names(y)[last_col_num]
  zdf[logkern] <- y[,last_col_num]+logJ
  
  #Returning important info
  list(logJ, invsqrts, colMeans(x), rmax)
}





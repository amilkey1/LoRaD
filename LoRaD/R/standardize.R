#' This function takes a data frame y that contains the log transformed parameter values
#' The log kernel is always the last column 
#' 
#' @param df data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @param coverage Fraction of training sample used to compute working parameter space
#' @return A list of standardize info of logJ, invsqrts, colMeans(x), rmax
#' @export 
#'
standardize <- function(df, coverage) {
  #cat("\nStandardizing Params:\n")
  
  #Getting dimensions of matrix
  n <- nrow(df)
  p <- ncol(df) - 1  # last column is not a parameter

  # Create matrix x from df by dropping log_kernel column
  last_col_num <- ncol(df)
  x <- as.matrix(df[,-last_col_num])
  
  # Calculating the mean vector
  meanvect <- colMeans(x)
  m <- matrix(rep(colMeans(x), n), byrow = TRUE, ncol = p)

  # Calculating the Variance-Covariance matrix
  s <-stats::cov(x)
  
  # Getting Eigen vectors and values
  lamb <- eigen(s)$value
  vec <- eigen(s)$vector
  
  # Calculating sqrt of s
  if (p == 1) {
    sqrts <- matrix(sqrt(lamb), ncol=1)
    invsqrts <- matrix(1/sqrt(lamb), ncol=1)
    logJ <- 0.5*log(lamb)
  }
  else {
    sqrts <- vec%*%diag(sqrt(lamb))%*%t(vec)
    invsqrts <- vec%*%diag(1/sqrt(lamb))%*%t(vec)
    logJ <- as.numeric(determinant(invsqrts,logarithm = TRUE)$modulus)
  }
    
  # Performing the standardization
  z <- (x-m)%*%invsqrts
  
  # make a vector of norms
  norms <- numeric()
  for (i in 1:n) {
    # vector of params
    v <- z[i,]
    
    # compute the radius (norm)
    r <- sqrt(sum(v^2))
    
    # append r to norms vector
    norms[i] <- r
  }
  
  # Converting z to dataframe
  zdf <- as.data.frame(z)
  
  # Append norms vector to z
  zdf["norm"] <- norms
  
  # Append log kernel vector to z
  zdf["logkernel"] <- df[,last_col_num] + logJ
  
  # sort z by norm
  zdfsorted <- zdf[order(zdf$norm),]

  # Retain only first n*coverage elements of sorted data frame
  n <- ceiling(n*coverage)
  zdf_cropped <- zdfsorted[1:n,]
  
  # rmax is the largest radius (norm) found in that portion of the
  # training sample retained
  rmax <- zdf_cropped$norm[[n]]
  
  #Returning important info
  list(logJ, invsqrts, meanvect, rmax)
}





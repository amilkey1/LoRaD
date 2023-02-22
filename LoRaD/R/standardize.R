#' This function takes a data frame y, that contains the log transform parameter values
#' The log kernal is always the last column 
#'
#' @y data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @return A new data frame consisting of standardize parameter values
#' @export 
#'
standardize <- function(y) {
  cat("\nStandardizing Params:\n")
  
  #Getting rid of log kernal bc not a param
  last_col_num <- ncol(y)
  x <- as.matrix[y[,-last_col_num]]
  #Getting dimensions of matrix
  n <- nrow(x)
  p <- ncol(x)
  
  #Mean time
  m <- matrix(rep(colMeans(x), n), byrow = TRUE, ncol = p)
  
  #Var-Cov Matrix
  s <-cov(x)
  
  #Getting Eigen vectors and values
  Lamb <- eigen(s)$value
  Vec <- eigen(s)$vector
  
  #Calculating sqrt of s
  sqrts <- Vec%*%diag(sqrt(Lamb))%*%t(Vec)
  invsqrts <- Vec%*%diag(1/sqrt(Lamb))%*%t(Vec)
  
  #Standardizing Stuff
  z <- (x-m)%*%invsqrts
  
  #Computing Log Jacobian
  logJ <- as.numeric(determinant(invsqrts,logarithm = TRUE)$modulus)
  
  #Converting z to dataframe
  zdf <- as.data.frame(z)
  
  #Adding Jacobian to log kernal
  logkern <- names(y)[last_col_num]
  zdf[logkern] <- y[,last_col_num]+logJ
  
  #Returning dataframe
  zdf
}





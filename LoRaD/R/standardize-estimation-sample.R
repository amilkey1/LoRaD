#' This function standardized the estimation sample
#'
#' \code {standardize_estimation_sample} A function that tajes the parameters standardinfo and y
#' @param standardinfo Contains list of logJ, invsqrts, colMeans(x), rmax
#' @param y  data frame containing a column for each transformed model parameter sampled 
#' @return A new data frame consisting of standardized estimation sample with log kernal in last column
#' @export 
#'
standardize_estimation_sample <- function(standardinfo, y) {
  #cat("\nStandardizing estimation sample:\n")
  
  #Getting rid of log kernel because it's not a parameter
  last_col_num <- ncol(y)
  x <- as.matrix(y[,-last_col_num])
  
  #Getting dimensions of matrix
  n <- nrow(x)
  p <- ncol(x)
  
  logJ <- standardinfo[[1]]
  invsqrts <- standardinfo[2]
  inversesqrts_matrix <- matrix(unlist(invsqrts), nrow=p, ncol=p)
  #Getting Mean Vec
  mean_vec <- as.matrix(unlist(standardinfo[3]), nrow=p, ncol=1)
  
  m <- matrix(rep(mean_vec, n), byrow = TRUE, ncol = p)
  
  
  #Performing the standardization
  z <- (x-m)%*%inversesqrts_matrix
  
  
  #Converting z to dataframe
  zdf <- as.data.frame(z)
  
  #Adding log Jacobian to log kernel to complete transformation
  logkern <- names(y)[last_col_num]
  zdf[logkern] <- y[,last_col_num]+logJ
  
  #Returning important info
  zdf 
  
}





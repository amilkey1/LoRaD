#' This function calculates the Lorad estimate of the marginal likelihood
#'
#' @params Parameters, log kernel is always the last column
#' @colspec Identification of columns
#' @trainingfrac Fraction of Samples Used in Training
#' @trainingmode Used Random, Left or Right
#' @return Lorad estimate of marginal likelihood
#' @export 
#'
lorad <- function(params, colspec, trainingfrac, trainingmode) {
  transformdf <- transform(params, colspec)
  nsamples <- nrow(transformdf)
  tmode <- tolower(trainingmode)
  
  # Specified training fraction must lie between 0 and 1
  if (training_frac <= 0 || training_frac >= 1.0) {
  	warning(sprintf("training fraction must be between 0 and 1 (%g)",training_frac))
  	stop()
  }
	
  y <- floor(trainingfrac*nsamples)
  if (tmode == "random") {
    x <- 1:nsamples
    z <- sample(x, y)
  }
   else if (tmode == "left") {
    z <- 1:y
   }
  else if (tmode == "right") {
    z <- y:nsamples
  }
  else {
    warning(sprintf("Unknown training mode (%s)", trainingmode))
    stop()
  }
  
  # Partition transformed samples into training and estimation samples
  trainingdf <- transformdf[z,]
  estimationdf <- transformdf[-z,] 
  standardinfo <- standardize(trainingdf)
  
  # standardinfo contains list(logJ, invsqrts, colmeans(x), rmax)
  
  # Printing out important info
  cat("\nPartitioning Samples Into Training and Estimation:\n")
  cat(sprintf("   Sample Size Is %d\n",nrow(transformdf)))
  # Printing out Training Info
  cat(sprintf("   Training Sample Size Is %d\n",nrow(trainingdf)))
  # Printing out Estimation Sample Size
  cat(sprintf("   Estimation Sample Size %d\n",nrow(estimationdf)))
  
  df <- standardize_estimation_sample(standardinfo, estimationdf)
  df
  
  # Store all but last column from df
  # Store last column as vector log posterior kernel
  
  last_col_num <- ncol(df)
  x <- as.matrix(df[,-last_col_num])
  logpostkern <- as.matrix(df[,last_col_num])
  
  # p is the number of parameters
  p <- ncol(x)
  
  # rmax is the maximum radius of any point in the training sample
  rmax <- standardinfo[[4]]
  
  cat("\nProcessing Training Sample:\n")
  cat(sprintf("   Lowest Radial Distance is %.5f\n",rmax))
  
  # sigma_squared = 1 for standard normal
  sigmasqr <- 1
  
  # Compute Delta, which is the integral from 0 to rmax of the marginal
  # distribution of radial vector lengths from a p-dimensional multivariate
  # standard normal distribution
  s <- p/2.0
  t <- rmax^2/(2.0*sigmasqr)
  logdelta <- pgamma(s,t,sigmasqr)
  cat(sprintf("   Log Delta %.5f\n",logdelta))
  
  #Calculating normalizing constant for reference function (multivariate std normal)
  sigma <- sqrt(sigmasqr)
  logmvnormconstant <- .5*p*log(2.0*pi)+1.0*p*log(sigma)
  # cat(sprintf("   Normalizing Constant for Reference Function %.5f\n",logmvnormconstant))
  
  
  #Initialize Variables
  logratios <- numeric()
  nestimation <- nrow(estimationdf)
  
  
  #Calculate sum of ratios, Using multivariate standard normal reference function 
  j <- 0.0
  for (i in 1:nestimation) {
    #Calculating norm for the ith sample
    v <- x[i,]
    r <- sqrt(sum(v^2))
    if (r<=rmax) {
      j <- j+1 
      logkernel <- logpostkern[i,]
      logreference <- .5*sigmasqr*r^2-logmvnormconstant
      logratio <- logreference-logkernel
      logratios[j] <- logratio
    }
  }
  logratios
}

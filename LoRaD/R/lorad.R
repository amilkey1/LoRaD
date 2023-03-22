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
  trainingdf <- transformdf[z,]
  estimationdf <- transformdf[-z,] 
  standardinfo <- standardize(trainingdf)
  #Printing out important info
  cat("\nPartitioning Samples Into Training and Estimation:\n")
  cat(sprintf("   Sample Size Is %d\n",nrow(transformdf)))
  #Printing out Training Info
  cat(sprintf("   Training Sample Size Is %d\n",nrow(trainingdf)))
  #Printing out Estimation Sample Size
  cat(sprintf("   Estimation Sample Size %d\n",nrow(estimationdf)))
  #list(logJ, invsqrts, colMeans(x), rmax)
  #standardinfo
  df <- standardize_estimation_sample(standardinfo, estimationdf)
  df
  
  #Store all but last column from Dataframe
  #Store last column as vector log posterior kernel
  last_col_num <- ncol(df)
  x <- as.matrix(df[,-last_col_num])
  logpostkern <- as.matrix(df[,last_col_num])
  #Getting number of params
  p <- ncol(x)
  #Getting max radius of any point in training sample
  rmax <- standardinfo[[4]]
  
  
  cat("\nProcessing Training Sample:\n")
  cat(sprintf("   Lowest Radial Distance is %.5f\n",rmax))
  
  #St Dev For St normal is one
  sigmasqr <- 1
  
  #Calculating Delta
  S <- p/2.0
  T <- rmax^2/(2.0*sigmasqr)
  logdelta <- pgamma(S,T,sigmasqr)
  cat(sprintf("   Log Delta %.5f\n",logdelta))
  #Calculating normalizing constant for reference function (St Normal)
  sigma <- sqrt(sigmasqr)
  logmvnormconstant <- .5*p*log(2.0*pi)+1.0*p*log(sigma)
  cat(sprintf("   Normalizing Constant for Reference Function %.5f\n",logmvnormconstant))
  
  
  #Initialize Variables
  logratios <- numeric()
  nestimation <- nrow(estimationdf)
  
  
  #Calculate Sum of Ratios, Using Multinormal reference function 
  J <- 0.0
  for (i in 1:nestimation) {
    #Calculating Norm for the Sample
    v <- x[i,]
    r <- sqrt(sum(v^2))
    if (r<=rmax) {
      J <- J+1 
      logkernel <- logpostkern[i,]
      
      logreference <- .5*sigmasqr*r^2-logmvnormconstant
      logratio <- logreference-logkernel
      logratios[J] <- logratio
    }
  }
  logratios
}

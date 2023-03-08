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
}

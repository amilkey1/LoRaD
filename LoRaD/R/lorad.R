#' This function takes a data frame y, that contains the log transform parameter values
#' The log kernel is always the last column 
#'
#' @y data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @return A new data frame consisting of standardize parameter values
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
  z
}
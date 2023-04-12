#' This function calculates the Lorad estimate of the marginal likelihood
#'
#' @params Parameters, log kernel is always the last column
#' @colspec Identification of columns
#' @training_frac Fraction of Samples Used in Training
#' @trainingmode Used Random, Left or Right
#' @coverage Fraction of training sample used to compute working parameter space
#' @file_name Name of parameter file
#' @return Lorad estimate of marginal likelihood
#' @export 
#'
lorad <- function(params, colspec, training_frac, training_mode, coverage, file_name) {
  cat("This is loRad (ver. 1.0):\n")
  cat(sprintf("   Parameter Sample File is : %s\n", file_name))
  transform_df <- transform(params, colspec)
  
  cat(sprintf("   Traning fraction is %g\n", training_frac))
  cat(sprintf("   Coverage specified is %g\n", coverage))
  
  nsamples <- nrow(transform_df)
  tmode <- tolower(training_mode)
    
  
  # Specified training fraction must lie between 0 and 1
  if (training_frac <= 0 || training_frac >= 1.0) {
  	warning(sprintf("training fraction must be between 0 and 1 (%g)",training_frac))
  	stop()
  }
	
  y <- floor(training_frac*nsamples)
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
    warning(sprintf("Unknown training mode (%s)", training_mode))
    stop()
  }
  
  # Partition transformed samples into training and estimation samples
  training_df <- transform_df[z,]
  estimation_df <- transform_df[-z,] 
  
   cat("Reading parameter sample file ...\n")
   cat(sprintf("   Processed %d column specifications \n", length(colspec)))
   cat(sprintf("   Found %d parameters \n", ncol(estimation_df)-1))
   cat(sprintf("   Found %d columns \n", length(params)))
   cat(sprintf("   File has %d lines \n", nrow(params)+1))
   cat(sprintf("   Found %d values for each column \n", nrow(params)))
  
  # Printing out important info
  cat("\nPartitioning samples into training and estimation:\n")
  cat(sprintf("   Sample size is %d\n",nrow(transform_df)))
  # Printing out Training Info
  cat(sprintf("   Training sample Size is %d\n",nrow(training_df)))
  # Printing out Estimation Sample Size
  cat(sprintf("   Estimation sample size %d\n",nrow(estimation_df)))
  
  # Extract just the parameters from estimation_df
  # Leave out last column (log posterior kernel values)
  
  last_col_num <- ncol(estimation_df)
  log_post_kern <- as.matrix(estimation_df[,last_col_num])
  
  # compute mean vector and inverse square root of covariance matrix
  # needed for transforming the estimation sample
  
  standard_info <- standardize(training_df, coverage)
  # standardinfo contains list(logJ, invsqrts, colmeans(x), rmax)
  logj <- standard_info[[1]]
  # cat(sprintf("   logj = %.5f\n",logj))
  df <- standardize_estimation_sample(standard_info, estimation_df)
 
  last_col_num <- ncol(estimation_df)
  x <- as.matrix(df[,-last_col_num])
  log_post_kern <- as.matrix(df[,last_col_num])
      
  # p is the number of parameters
  p <- ncol(x)
  
  # rmax is the maximum radius of any point in the training sample
  rmax <- standard_info[[4]]
  
  cat("\nProcessing training sample...\n")
  cat(sprintf("   Lowest radial distance is %.5f\n",rmax))
  
  # sigma_squared = 1 for standard normal
  sigma_sqr <- 1
  
  # Compute Delta, which is the integral from 0 to rmax of the marginal
  # distribution of radial vector lengths from a p-dimensional multivariate
  # standard normal distribution
  s <- p/2.0
  t <- rmax^2/(2.0*sigma_sqr)
  log_delta <- log(stats::pgamma(t, shape=s, scale=1))
  # cat(sprintf("   Log Delta %.5f\n",log_delta))
  
  #Calculating normalizing constant for reference function (multivariate std normal)
  sigma <- sqrt(sigma_sqr)
  log_mvnorm_constant <- .5*p*log(2.0*pi)+1.0*p*log(sigma)
  # cat(sprintf("   Normalizing Constant for Reference Function %.5f\n",log_mvnorm_constant))
  
  
  #Initialize Variables
  log_ratios <- numeric()
  nestimation <- nrow(estimation_df)
  
  #Calculate sum of ratios, Using multivariate standard normal reference function 
  j <- 0.0
  for (i in 1:nestimation) {
    #Calculating norm for the ith sample
    v <- x[i,]
    r <- sqrt(sum(v^2))
    if (r<=rmax) {
      j <- j+1 
      log_kernel <- log_post_kern[i]
      log_reference <- -.5*sigma_sqr*r^2 - log_mvnorm_constant
      log_ratio <- log_reference - log_kernel
      log_ratios[j] <- log_ratio
    }
  }
 
  # log_ratios
  cat("\nProcessing estimation sample...\n")
  cat(sprintf("   Number of samples used is %d\n",j))
  cat(sprintf("   Nominal coverage specified is %f\n",coverage))
  cat(sprintf("   Actual coverage is %f\n",j/nestimation))
  
  if (length(log_ratios)==0){
    warning(sprintf("No estimation samples were within the working parameter space (rmax=%g)",rmax))
    stop()
  }
  
  log_sum_ratios <- calcLogSum(log_ratios)
  #Calculate LoRaD estimate of maximum likelihood
  log_marginal_likelihood <- log_delta - (log_sum_ratios - log(nestimation))
  cat(sprintf("   Log marginal likelihood is %f\n",log_marginal_likelihood))
  
}

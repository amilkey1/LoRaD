#' Calculates the LoRaD estimate of the marginal likelihood
#'
#' @param params Data frame in which rows are sample points and columns are parameters, except that last column holds the log posterior kernel
#' @param colspec Named character vector associating column names in params with column specifications
#' @param training_frac Number between 0 and 1 specifying the training fraction
#' @param training_mode One of random, left, or right, specifying how training fraction is chosen
#' @param coverage Number between 0 and 1 specifying fraction of training sample used to compute working parameter space
#' @return LoRaD estimate of log marginal likelihood
#' @export 
#'
lorad_estimate <- function(params, colspec, training_frac, training_mode, coverage) {
    cat("This is loRad (ver. 1.0):\n")

    # Transform any parameters that are constrained and consolidate log kernel components
    # into a single column named log_kernel that includes Jacobian terms for transformations
    transform_df <- lorad_transform(params, colspec)
    
    nsamples <- nrow(transform_df)
    nparams <- ncol(transform_df) - 1
    tmode <- tolower(training_mode)

    # Specified training fraction must lie between 0 and 1
    if (training_frac <= 0 || training_frac >= 1.0) {
        warning(sprintf("training fraction must be between 0 and 1 (%g)",training_frac))
        stop()
    }
	
	# Determine sites included in training sample and place remainder in estimation sample
    y <- floor(training_frac*nsamples)
    training_included_sites_str <- "randomly chosen"
    estim_included_sites_str <- "complement of training site set"
    if (tmode == "random") {
        x <- 1:nsamples
        z <- sample(x, y)
    }
    else if (tmode == "left") {
        z <- 1:y
        training_included_sites_str <- sprintf("sites 1 to %d", y)
        estim_included_sites_str <- sprintf("sites %d to %d", y+1, nsamples)
    }
    else if (tmode == "right") {
        z <- y:nsamples
        training_included_sites_str <- sprintf("sites %d to %d", y, nsamples)
        estim_included_sites_str <- sprintf("sites 1 to %d", y)
    }
    else {
        warning(sprintf("Unknown training mode (%s)", training_mode))
        stop()
    }

    # Provide feedback to user  
    cat(sprintf("   Parameter sample comprises %d sampled points\n", nsamples))
    cat(sprintf("   Each sample point is a vector of %d parameter values\n", nparams))
    cat(sprintf("   Training fraction is %g\n", training_frac))
    cat(sprintf("   Training fraction mode is %s\n", tmode))
    cat(sprintf("   Coverage specified is %g\n", coverage))

    # Partition transformed samples into training and estimation samples
    training_df <- transform_df[z,]
    estimation_df <- transform_df[-z,] 

    # Print out information about training and estimation samples
    cat("\nPartitioning samples into training and estimation:\n")
    cat(sprintf("   Sample size is %d\n",nrow(transform_df)))
    cat(sprintf("   Training sample size is %d (%s)\n", nrow(training_df), training_included_sites_str))
    cat(sprintf("   Estimation sample size %d (%s)\n", nrow(estimation_df), estim_included_sites_str))
    
    # Store the log kernel values in a vector
    last_col_num <- ncol(estimation_df)
    log_post_kern <- as.array(estimation_df[,last_col_num])

    # Compute mean vector and inverse square root of covariance matrix
    # needed for transforming the estimation sample. standard_info is 
    # a list with elements logJ, invsqrts, colmeans(x), and rmax
    standard_info <- lorad_standardize(training_df, coverage)
    logj <- standard_info[[1]]
    
    # Now standardize the estimation sample using mean vector and 
    # standard deviation matrix estimated from the training sample
    df <- lorad_standardize_estimation_sample(standard_info, estimation_df)
 
    # Extract posterior kernel
    last_col_num <- ncol(estimation_df)
    x <- as.matrix(df[,-last_col_num])
    log_post_kern <- df[,last_col_num]
  
    # The number of parameters
    p <- ncol(x)

    # The maximum radius of any point in the training sample
    rmax <- standard_info[[4]]

    cat("\nProcessing training sample...\n")
    cat(sprintf("   Lowest radial distance is %.9f\n",rmax))

    # The variance is 1 for standard normal
    sigma_sqr <- 1
  
    # Compute Delta, which is the integral from 0 to rmax of the marginal
    # distribution of radial vector lengths from a p-dimensional multivariate
    # standard normal distribution
    s <- p/2.0
    t <- rmax^2/(2.0*sigma_sqr)
    log_delta <- log(stats::pgamma(t, shape=s, scale=1))

    cat(sprintf("   Log Delta %.9f\n",log_delta))

    # Calculate the normalizing constant for reference function (multivariate standard normal)
    sigma <- sqrt(sigma_sqr)
    log_mvnorm_constant <- .5*p*log(2.0*pi)+1.0*p*log(sigma)
    
    # Initialize variables
    log_ratios <- numeric()
    nestimation <- nrow(estimation_df)
  
    # Calculate sum of ratios, Using multivariate standard normal reference function 
    j <- 0.0
    for (i in 1:nestimation) {
        # Calculate norm for the ith sample
        v <- x[i,]
        r <- sqrt(sum(v^2))
        
        # If norm is in working parameter space, use it
        if (r <= rmax) {
            j <- j+1 
            log_kernel <- log_post_kern[i]
            log_reference <- -.5*sigma_sqr*r^2 - log_mvnorm_constant
            log_ratio <- log_reference - log_kernel
            log_ratios[j] <- log_ratio
        }
    }
 
    cat("\nProcessing estimation sample...\n")
    cat(sprintf("   Number of samples used is %d\n",j))
    cat(sprintf("   Nominal coverage specified is %f\n",coverage))
    cat(sprintf("   Actual coverage is %f\n",j/nestimation))
  
    if (length(log_ratios)==0){
        warning(sprintf("No estimation samples were within the working parameter space (rmax=%g)",rmax))
        stop()
    }
  
    log_sum_ratios <- lorad_calc_log_sum(log_ratios)

    #Calculate LoRaD estimate of maximum likelihood
    log_marginal_likelihood <- log_delta - (log_sum_ratios - log(nestimation))
    cat(sprintf("   Log marginal likelihood is %f\n",log_marginal_likelihood))
}

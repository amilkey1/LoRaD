#' Log (or log-ratio) transform parameters having constrained support
#' 
#' Log-transforms parameters with support (0,infinity), log-ratio transforms 
#' K-dimensional parameters with support a (K-1)-simplex, logit transforms 
#' parameters with support [0,1], and leaves unchanged parameters with 
#' unconstrained support (-infinity, infinity).
#' 
#' @param params Data frame containing a column for each model parameter 
#'     sampled as well as one or more columns that, when summed, constitute 
#'     the log joint posterior kernel
#' @param colspec Named character vector matching each column name in 
#'     params with a column specification
#' @returns A new data frame comprising transformed parameter values
#'     with a final column holding the log joint posterior kernel
#' @examples
#' # A data frame with a column of normally-distributed parameter values,
#' # a column of beta-distributed parameter values, a column of gamma-distributed
#' # parameter values, and three columns representing Dirichlet-distributed
#' # parameter values, is log transformed. The columns containing the 
#' # log densities of each parameter are summed and the total becomes the last
#' # column in the returned data frame. If this were a real example, the columns
#' # specified as `posterior` would represent components of the log posterior  
#' # kernel (i.e. the unnormalized posterior, likelihood*prior).
#' normals <- rnorm(5,0,10)
#' prob_normals <- dnorm(normals,0,10,log=TRUE) 
#' proportions <- rbeta(5,1,2)
#' prob_proportions <- dbeta(proportions,1,2,log=TRUE)
#' lengths <- rgamma(5, 10, 1)
#' prob_lengths <- dgamma(lengths,10,1,log=TRUE)
#' freqA <- c(0.6, 0.2, 0.1, 0.5, 0.2)
#' freqB <- c(0.3, 0.1, 0.8, 0.1, 0.5)
#' freqC <- c(0.1, 0.7, 0.1, 0.4, 0.3)
#' prob_freqs <- 0.0  # assuming Dirichlet(1,1,1)
#' paramsdf <- data.frame(
#'     normals,prob_normals,
#'     proportions, prob_proportions,
#'     lengths, prob_lengths,
#'     freqA, freqB, freqC, prob_freqs)
#' columnkey <- c(
#'     "normals"="unconstrained", 
#'     "prob_normals"="posterior", 
#'     "proportions"="proportion", 
#'     "prob_proportions"="posterior", 
#'     "lengths"="positive", 
#'     "prob_lengths"="posterior", 
#'     "freqA"="simplex", 
#'     "freqB"="simplex", 
#'     "freqC"="simplexfinal",
#'     "prob_freqs"="posterior")
#' print(lorad_transform(paramsdf, columnkey))

lorad_transform <- function(params, colspec) {
	num_rows <- nrow(params)
	num_columns <- ncol(params)
    col_names <- colnames(params)
    if (length(colspec)>num_columns){
      stop("colspec does not match column names - colspec is too long")
    }
    else if (length(colspec)<num_columns){
      stop("colspec does not match column names - colspec is too short")
    }

    # Create an empty data frame with same number of rows as params 
    # in which to hold transformed parameter values
    df_col_names <- c()
	df <- as.data.frame(matrix(nrow=nrow(params),ncol=0))
	
    # Initialize log_kernel as a vector of zeros
	log_kernel <- vector(mode="numeric", length=num_rows) # params$log.kernel
	
	# OK will only be false if user has failed to specify a column type for some column
    OK <- TRUE
    UNKNOWN_COLSPEC <- ""
    
    #U1 <- 0.0
    #FIRSTSIMPLEX <- FALSE
    #sumsimplex <- 0.0 
    #FINALSIMPLEX <- FALSE
    
    U1 <- NULL
    sumsimplex <- NULL
    INSIMPLEX <- FALSE

    # Transform (if positive) or copy (if unconstrained) each column into df
    for (i in 1:num_columns) {
    	column_type <- colspec[col_names[i]]
    	if (is.na(column_type)) {
    	    # User failed to supply a column spec for this column
    		OK <- FALSE
    	}
    	else {
    	    if (column_type == "positive") {
    	        # Create vector of log-transformed parameter values
				transformed <- log(params[[i]])
				
				# Add the log-Jacobian for the transformation to log_kernel for each row
				# The log-Jacobian for a log transformation is just the transformed value
				log_kernel <- log_kernel + transformed
				
				# Append transformed to the growing data frame
				df <- cbind(df, transformed)
				
				# Append the column name to the stored column names
				df_col_names <- cbind(df_col_names, col_names[i])
            }
            else if (column_type == "proportion") {
                # Create vector of log-transformed parameter values
                transformed <- log(params[[i]])-log(1-params[[i]])
    	    
                log_J <- log(params[[i]])+log(1-params[[i]])
                # Add the log-Jacobian for the transformation to log_kernel for each row
                # The log-Jacobian for a log transformation is just the log(P)+log(1-P) value
                log_kernel <- log_kernel + log_J
    	    
                # Append transformed to the growing data frame
                df <- cbind(df, transformed)
    	    
                # Append the column name to the stored column names
                df_col_names <- cbind(df_col_names, col_names[i])
            }
            else if (column_type == "simplex") {
    	        logUi <- log(params[[i]])
    	        if (INSIMPLEX) {
                    # Create vector of log-transformed parameter values
                    transformed <- logUi - logU1

                    # Add the log-Jacobian for the transformation to log_kernel for each row
                    # The log-Jacobian for a log transformation is just the log(U) value
                    log_kernel <- log_kernel + logUi

                    # Append transformed to the growing data frame
                    df <- cbind(df, transformed)

                    # Append the column name to the stored column names
                    df_col_names <- cbind(df_col_names, col_names[i])
                    
                    # Add to sumsimplex
                    sumsimplex <- sumsimplex + params[[i]]
    	        }
    	        else {
    	            INSIMPLEX = TRUE
    	            logU1 <- logUi
    	            
                    # Add the log-Jacobian for the transformation to log_kernel for each row
                    # The log-Jacobian for a log transformation is just the log(U) value
                    log_kernel <- log_kernel + logUi

                    # Initialize sumsimplex
                    sumsimplex <- params[[i]]
    	        }
            }
            else if (column_type == "simplexfinal") {
    	        logUi <- log(params[[i]])
    	        if (INSIMPLEX) {
                    # Create vector of log-transformed parameter values
                    transformed <- logUi - logU1

                    # Add the log-Jacobian for the transformation to log_kernel for each row
                    # The log-Jacobian for a log transformation is just the log(U) value
                    log_kernel <- log_kernel + logUi

                    # Append transformed to the growing data frame
                    df <- cbind(df, transformed)

                    # Append the column name to the stored column names
                    df_col_names <- cbind(df_col_names, col_names[i])
                    
                    # Add to sumsimplex
                    sumsimplex <- sumsimplex + params[[i]]
                                        
                    # Check whether simplex columns add to 1 for each row
                    max_sumsimplex <- max(abs(1 - sumsimplex))
                    if (max_sumsimplex > 0.00001) {
                        stop("simplex columns do not add to 1 for all rows")
                    }
                    
                    INSIMPLEX <- FALSE
                    
                } else {
                    stop("simplexfinal column found with no preceding simplex columns")
                }    	        
            }
			else if (column_type == "unconstrained") {
				# Append the vector of original parameter values to the growing data frame
				# No transformation necessary since parameter is already unconstrained
				df <- cbind(df, params[[i]])
				
				# Append the column name to the stored column names
				df_col_names <- cbind(df_col_names, col_names[i])
			}
			else if (column_type == "posterior") {
			    # Add the values in this column to log_kernel since the column spec
			    # says that this column is a component of the log kernel
				log_kernel <- log_kernel + params[[i]]
			}
			else if (column_type == "iteration") {
			    # Ignore iteration columns
			}
			else if (column_type == "ignore") {
			    # Ignore ignore columns
			}
			else {
			  # User supplied an unknown column spec
			  stop(sprintf("unknown colspec found: %s", column_type))
			}
    	}
    }
    
    # Now that log_kernel takes account of the original log kernel as well as the
    # Jacobians from all transformations, add log_kernel as the final column in df
	df <- cbind(df, log_kernel)
	df_col_names <- cbind(df_col_names, "log_kernel")
	
	# Attach column names to df
	names(df) <- df_col_names
    
    if (!OK) {
        stop("colspec does not match column names")
    }
    
    # Return the transformed parameter data frame
    df
}

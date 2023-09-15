#' Takes a data frame consisting of sampled parameter values
#' and a named character vector associating a column name
#' in the data frame to a column specification and 
#' returns a new data frame containing parameters that 
#' have been log-transformed if originally constrained.
#' 
#' @param params Data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @param colspec Named character vector associating each column name in params to a column specification
#' @return A new data frame consisting of transformed parameter values
#' @export
#'
lorad_transform <- function(params, colspec) {
	num_rows <- nrow(params)
	num_columns <- ncol(params)
    col_names <- colnames(params)

    # Create an empty data frame with same number of rows as params 
    # in which to hold transformed parameter values
    df_col_names <- c()
	df <- as.data.frame(matrix(nrow=nrow(params),ncol=0))
	
    # Initialize log_kernel as a vector of zeros
	log_kernel <- vector(mode="numeric", length=num_rows) # params$log.kernel
	
	# OK will only be false if user has failed to specify a column type for some column
    OK <- TRUE

    # Transform (if nonnegative) or copy (if unconstrained) each column	into df
    for (i in 1:num_columns) {
    	column_type <- colspec[col_names[i]]
    	if (is.na(column_type)) {
    	    # User failed to supply a column spec for this column
    		OK <- FALSE
    	}
    	else {
    	    if (column_type == "nonnegative") {
    	        # Create vector of log-transformed parameter values
				transformed <- log(params[c(i)])
				
				# Add the log-Jacobian for the transformation to log_kernel for each row
				# The log-Jacobian for a log transformation is just the transformed value
				log_kernel <- log_kernel + transformed
				
				# Append transformed to the growing data frame
				df <- cbind(df, transformed)
				
				# Append the column name to the stored column names
				df_col_names <- cbind(df_col_names, col_names[i])
			}
			else if (column_type == "unconstrained") {
				# Append the vector of original parameter values to the growing data frame
				# No transformation necessary since parameter is already unconstrained
				df <- cbind(df, params[c(i)])
				
				# Append the column name to the stored column names
				df_col_names <- cbind(df_col_names, col_names[i])
			}
			else if (column_type == "posterior") {
			    # Add the values in this column to log_kernel since the column spec
			    # says that this column is a component of the log kernel
				log_kernel <- log_kernel + params[c(i)]
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
        warning("colspec does not match column names")
        stop()
        # print("error")
    }
    
    # Return the transformed parameter data frame
    df
}

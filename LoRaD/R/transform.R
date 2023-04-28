#' This function takes a params data frame consisting of sampled parameter values 
#' and a colspec data frame that specifies what each column in the params
#' data frame represents.
#'
#' \code {transform} Function that uses parameters params and colspec
#' 
#' @param params Data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @param colspec Dictionary classifying each column of the params data frame 
#' @return A new data frame consisting of transformed parameter values
#' @export
transform <- function(params, colspec) {
	num_columns <- ncol(params)
	#print(sprintf("%d",num_columns))
    col_names <- colnames(params)
    OK <- TRUE
	log_kernel <- params$log.kernel
	#print(length(colspec))
	df_col_names <- c()
	df <-as.data.frame(matrix(nrow=nrow(params),ncol=0))
    for (i in 1:num_columns) {
    	value <- colspec[col_names[i]]
    	if (is.na(value)) {
    		OK <- FALSE
    	}
    	else {
    	    if (value == "positive") {
				transformed <- log(params[c(i)])
				log_kernel = log_kernel+transformed
				df <- cbind(df, transformed)
				df_col_names <- cbind(df_col_names, col_names[i])
			}
			else if (value == "unconstrained") {
				df <- cbind(df, params[c(i)])
				df_col_names <- cbind(df_col_names, col_names[i])
			}
    	}
    	#print(sprintf("%20s %20s", col_names[i], value))
    }
    
	df <- cbind(df, log_kernel)
	df_col_names <- cbind(df_col_names, "log_kernel")
	names(df) <- df_col_names
    
    if (!OK) {
    	print("error")
    }
    
    df
}

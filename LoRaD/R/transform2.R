# this function will log-transform "positive" parameters

transform <- function(df){
	cat("hello")
	log_edge_len <- log(df$edgelen)
	log_kappa <- log(df$kappa)
	log_kernel <- df$log.kernel+log(df$edgelen)+log(df$kappa)
	transformed_df <- data.frame(log_kernel, log_edge_len, log_kappa)
	transformed_df
}
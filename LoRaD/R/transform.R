#' Calculate the log of the marginal likelihood
#'
#' This function takes a params data frame consisting of sampled parameter values 
#' and a colspec data frame that specifies what each column in the params
#' data frame represents.
#'
#' @params data frame containing a column for each model parameter sampled as well as columns that constitute the log posterior kernel
#' @colspec dictionary classifying each column of the params data frame
#' @return A new data frame consisting of transformed parameter values
#' @export
logmarglike <- function(params) {
    cat("Hello\n")
    log_kappa <- log(params$kappa)
    log_edgelen <- log(params$edgelen)
    log_kernel <- params$log.kernel + log(params$kappa) + log(params$edgelen)
    df <- data.frame(log_kernel, log_edgelen, log_kappa)
    test1 <- colspec[0]
    df
}

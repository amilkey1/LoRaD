#' This function calculates log of a sum without leaving the log scale
#'
#'
calcLogSum <- function(logx) {
  n <- length(logx)
  max_logx <- max(logx)
  sum_of_terms <- 0
  for (i in 1:n) {
    sum_of_terms <- sum_of_terms + exp(logx[[i]] - max_logx)
    
  }
  log_sum <- log(sum_of_terms) + max_logx
  log_sum
}





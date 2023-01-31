#' A Rabbit Function
#'
#' This function allows you to express your love of rabbits.
#' @param love Do you love rabbits? Defaults to TRUE.
#' @keywords rabbits
#' @export
#' @examples
#' rabbit_function()

rabbit_function <- function(love=TRUE){
    if(love==TRUE){
        print("I love rabbits!")
    }
    else {
        print("I am not a cool person.")
    }
}
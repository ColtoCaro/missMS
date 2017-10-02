#R functions for missMS

#' Generate a random number for Gaussian mean parameter
#'
#' @param n_ The number of random variates to create
#' @param xi esn parameter
#' @param omega esn parameter
#' @param alpha esn parameter
#' @param tau esn parameter
#'
#' @export

testfun <- function(n_, xi){
  n_ + xi
}

#' Calculation of the mode of a dataset
#'
#' @param x a vector
#'
#' @return mode
#' @export
#'
#' @examples
#' Mode(rnorm(80,20,3))
Mode <- function(x) {
  densite <- density(x)
  mode <- densite$x[which(densite$y==max(densite$y))]
  return(mode)}

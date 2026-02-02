#' Discrete Variable Classification Test
#'
#' This function tests whether a vector can be classified as discrete
#' using a criterion based on the number of unique values relative to
#' the total size of the vector. A variable is considered discrete if
#' the number of unique values is less than or equal to the square root
#' of the vector's length.
#'
#' @param vector A numeric, integer, factor, or character vector.
#'
#' @return A boolean:
#'   - `TRUE`: If the variable is classified as discrete.
#'   - `FALSE`: If the variable is classified as continuous.
#'
#' @details
#' The logic behind this function relies on the following criterion:
#' A variable is considered discrete if \code{length(unique(vector))} is less
#' than or equal to \code{floor(sqrt(length(vector)))}, where \code{length(vector)}
#' is the total size of the vector. This rule dynamically adapts to the vector's size.
#'
#' @examples
#' # Example with a discrete vector
#' vec1 <- c(1, 1, 2, 2, 3, 3, 3)
#' discret.test(vec1) # Returns TRUE
#'
#' # Example with a continuous vector
#' vec2 <- seq(0, 1, length.out = 100)
#' discret.test(vec2) # Returns FALSE
#'
#' # Example with a character vector
#' vec3 <- c("A", "B", "A", "C", "A", "B")
#' discret.test(vec3) # Returns TRUE
#'
#' @seealso \code{\link{unique}}, \code{\link{length}}, \code{\link{sqrt}}
#'
#' @export
discret.test <- function(vector) {
	if(length(unique(vector))<floor(sqrt(length(vector)))) {
		return(TRUE)
	}else {
		return(FALSE)
	}
}

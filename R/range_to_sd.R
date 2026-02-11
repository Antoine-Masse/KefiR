#' Estimate standard deviation from a confidence range
#'
#' @description
#' Computes an estimated standard deviation (`sd`) assuming that
#' the observed range `[min, max]` corresponds to a given confidence level
#' under a normal distribution.
#'
#' @param min Numeric. Lower bound of the observed range.
#' @param max Numeric. Upper bound of the observed range.
#' @param conf.level Numeric. Confidence level of the range (default = 0.99).
#'
#' @details
#' The function assumes the range represents approximately `±z·σ`,
#' where `z` is the quantile of the standard normal distribution
#' associated with the given confidence level.
#'
#' @return
#' A numeric value corresponding to the estimated standard deviation.
#'
#' @examples
#' range_to_sd(10, 20, conf.level = 0.95)
#'
#' @export
range_to_sd <- function(min, max, conf.level = 0.99) {
  z <- qnorm((1 + conf.level) / 2)
  (max - min) / (2 * z)
}

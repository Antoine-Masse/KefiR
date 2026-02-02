#' Generate a new variable correlated with the first column of a dataframe
#'
#' @description
#' Creates a new numeric variable correlated with the first column of a given dataframe.
#' The user can specify the desired correlation coefficient, mean, and standard deviation
#' of the new variable, as well as options to enforce exact or approximate correlation and statistics.
#'
#' @param dataframe A data frame containing at least one numeric column.
#' @param cor.f Numeric, between -1 and 1. Desired correlation coefficient with the first column.
#' @param meani Optional. Desired mean for the new variable. Default is 0.
#' @param sdi Optional. Desired standard deviation for the new variable. Default is 1.
#' @param strict Logical. If `TRUE`, enforces exactly the specified mean and sd;
#'               if `FALSE`, allows natural fluctuations. Default is `FALSE`.
#' @param strict_corr Logical. If `TRUE`, enforces exactly the desired correlation;
#'                    if `FALSE`, allows minor fluctuations. Default is `FALSE`.
#'
#' @return A numeric vector of the same length as the number of rows in the dataframe,
#' representing the generated correlated variable.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Converts the input to a matrix.
#'   \item Generates a normal random variable `z` with mean `meani` and sd `sdi`.
#'   \item Centers and standardizes the first column.
#'   \item If `strict_corr = TRUE`, orthogonalizes `z` relative to the first column for exact correlation.
#'   \item Combines both standardized vectors using `cor.f` and `sqrt(1 - cor.f^2)`.
#'   \item If `strict = TRUE`, rescales to exactly match `meani` and `sdi`.
#' }
#'
#' @examples
#' df <- data.frame(x = rnorm(100))
#' y <- cor.e(df, cor.f = 0.5, meani = 10, sdi = 2)
#' cor(df$x, y)
#' mean(y); sd(y)
#'
#' @export
cor.e <- function(dataframe, cor.f, meani = 0, sdi = 1, strict = FALSE, strict_corr = FALSE) {
  if (is.vector(dataframe)) dataframe <- data.frame(V1 = dataframe)
  if (!is.data.frame(dataframe)) stop("Input must be a dataframe or vector.")
  if (!is.numeric(dataframe[[1]])) stop("First column must be numeric.")
  if (cor.f < -1 || cor.f > 1) stop("Correlation must be between -1 and 1.")
  if (sdi <= 0) stop("Standard deviation must be positive.")

  mat.ini <- as.matrix(dataframe)
  z <- rnorm(nrow(mat.ini), mean = meani, sd = sdi)
  z_cent <- z - mean(z)
  x1 <- mat.ini[, 1]
  x1_cent <- x1 - mean(x1)
  x1_std <- x1_cent / sd(x1_cent)

  if (strict_corr) {
    proj_coef <- sum(z_cent * x1_cent) / sum(x1_cent^2)
    proj <- proj_coef * x1_cent
    z_orth <- z_cent - proj
    z_orth_std <- z_orth / sd(z_orth)
    new_var <- cor.f * x1_std + sqrt(1 - cor.f^2) * z_orth_std
  } else {
    z_std <- z_cent / sd(z_cent)
    new_var <- cor.f * x1_std + sqrt(1 - cor.f^2) * z_std
  }

  if (strict) {
    new_var <- (new_var - mean(new_var)) / sd(new_var)
    new_var <- new_var * sdi + meani
  } else {
    new_var <- new_var * sd(z) + mean(z)
  }

  new_var
}

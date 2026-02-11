#' Print method for MANOVA results
#'
#' Provides a concise summary of MANOVA results without displaying
#' raw data matrices which can be very verbose.
#'
#' @param x A MANOVA result object from m.test()
#' @param ... Additional arguments (not used)
#'
#' @details
#' This print method hides the raw data (x, g) and large objects
#' (manova_model, lda, predictions) to keep the output readable.
#'
#' To access full data: use str(result) or result$x, result$g directly.
#'
#' @export
print.manova_result <- function(x, ...) {

  cat("================================================================================\n")
  cat("MANOVA RESULT\n")
  cat("================================================================================\n\n")

  # Data summary (instead of full data)
  if (!is.null(x$x) && is.matrix(x$x)) {
    cat("Data:\n")
    cat("  Observations: ", nrow(x$x), "\n", sep = "")
    cat("  Dependent variables: ", ncol(x$x), " (",
        paste(colnames(x$x), collapse = ", "), ")\n", sep = "")
  }

  if (!is.null(x$g)) {
    cat("  Groups: ", nlevels(as.factor(x$g)), " levels\n", sep = "")
  }
  cat("\n")

  # Assumptions
  cat("Assumptions checked:\n")
  cat("  Multivariate normality: ",
      ifelse(isTRUE(x$check_multivariate_normality), "PASSED", "VIOLATED"), "\n", sep = "")
  cat("  Covariance homogeneity: ",
      ifelse(isTRUE(x$check_covariance_homogeneity), "PASSED", "VIOLATED"), "\n", sep = "")
  cat("  Multicollinearity: ",
      ifelse(isTRUE(x$check_multicollinearity), "PASSED", "VIOLATED"), "\n", sep = "")
  cat("  Multivariate outliers: ",
      ifelse(isTRUE(x$check_outliers_mv), "PASSED", "VIOLATED"), "\n", sep = "")
  cat("\n")

  # Test statistics (compact format)
  if (!is.null(x$test_statistics)) {
    cat("Test statistics:\n")

    if (!is.null(x$test_statistics$Wilks)) {
      cat("  Wilks' Lambda: ", round(x$test_statistics$Wilks$statistic, 4),
          " (p = ", .format_pval(x$test_statistics$Wilks$p_value), ")\n", sep = "")
    }

    if (!is.null(x$test_statistics$Pillai)) {
      cat("  Pillai's trace: ", round(x$test_statistics$Pillai$statistic, 4),
          " (p = ", .format_pval(x$test_statistics$Pillai$p_value), ")\n", sep = "")
    }

    cat("\n")
  }

  # Post-hoc summary
  if (!is.null(x$posthoc_manova)) {
    cat("Post-hoc analyses:\n")

    # Discriminant analysis
    if (!is.null(x$posthoc_manova$discriminant_analysis)) {
      da <- x$posthoc_manova$discriminant_analysis
      cat("  Discriminant analysis: ", da$n_functions, " function(s)\n", sep = "")
      if (!is.null(da$prop_variance) && length(da$prop_variance) > 0) {
        cat("    LD1 explains ", round(da$prop_variance[1] * 100, 1), "% variance\n", sep = "")
      }
    }

    # Protected ANOVAs
    if (!is.null(x$posthoc_manova$protected_anovas)) {
      n_sig <- sum(sapply(x$posthoc_manova$protected_anovas, function(a) a$significant))
      n_total <- length(x$posthoc_manova$protected_anovas)
      cat("  Protected ANOVAs: ", n_sig, "/", n_total, " DVs significant\n", sep = "")
    }

    cat("\n")
  }

  # Access instructions
  cat("Use str(result) to see full structure\n")
  cat("Access data: result$x, result$g\n")
  cat("Access post-hocs: result$posthoc_manova\n")
  cat("================================================================================\n")

  invisible(x)
}

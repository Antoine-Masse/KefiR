.format_pval <- function(p, digits = 2, threshold = 1e-6) {
  if (p < 0 || p > 1) {
    stop("Invalid p-value: must be between 0 and 1.")
  }
  if (abs(p) < threshold) {
    # Écriture scientifique pour les petites valeurs
    formatC(p, format = "e", digits = digits)
  } else {
    # Écriture classique avec arrondi
    signif(p, digits = digits)
  }
}
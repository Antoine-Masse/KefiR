# Définition d'une classe S3 avec un print() personnalisé
print.posthoc <- function(x, ...) {
  cat(.msg("\n=== Results ===\n","\n=== Résultats ===\n"))
  print(x$groups)
  #cat("\nP-value :", format(x$p.value, scientific = TRUE), "\n")

  if (!is.null(x$bootstrap)) {
    cat("\n=== Bootstrap ===\n")
    print(x$bootstrap$groups)
    #cat("\nP-value Bootstrap (95% CI) :", format(x$bootstrap$p.value["95%"], scientific = TRUE), "\n")
  }
}	
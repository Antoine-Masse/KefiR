# Définition d'une classe S3 avec un print() personnalisé
.print_table_with_title <- function(title, tbl) {
  cat(title, "\n")
  print(tbl)
}

#' @export
print.posthoc <- function(x, ...) {
  # Déterminer le nombre de groupes
  n_groups <- nrow(x$groups)

  # Pour 2 groupes : pas d'affichage du titre (résultats déjà présentés dans .posthoc())
  # Pour >2 groupes : afficher le titre
  if (n_groups > 2) {
    cat(.msg("\nPairwise post-hoc comparisons:\n","\nComparaisons post-hoc par paires :\n"))
  }

  print(x$groups)
  #cat("\nP-value :", format(x$p.value, scientific = TRUE), "\n")

  if (!is.null(x$bootstrap)) {
    .print_table_with_title("=== Bootstrap ===", x$bootstrap$groups)
    #cat("\nP-value Bootstrap (95% CI) :", format(x$bootstrap$p.value["95%"], scientific = TRUE), "\n")
  }
}	

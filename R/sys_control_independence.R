#' Contrôle de l'indépendance entre facteurs et ajout d'interactions
#'
#' Cette fonction analyse les paires de variables qualitatives (facteurs)
#' parmi les prédicteurs d'un modèle donné, afin de vérifier leur dépendance
#' statistique à l'aide du \code{\link[DescTools]{GTest}} (test du Khi-deux/G-test).
#' Pour chaque paire de facteurs dont la p-value (corrigée par la méthode de Holm)
#' est inférieure au seuil \code{alpha}, une interaction pure \code{Var_i:Var_j}
#' est ajoutée à la partie droite de la formule (si elle n'est pas déjà présente).
#'
#' @param formula Une formule R décrivant le modèle, par exemple \code{Note ~ Lycée*Sexe + Options}.
#' @param g Un \code{data.frame} ou \code{tibble} contenant les variables
#'          décrites dans \code{formula}.
#' @param alpha Seuil de signification pour les p-values corrigées. Défaut : \code{0.05}.
#' @param debug Affiche des informations de débogage lorsque \code{TRUE}. Défaut : \code{FALSE}.
#'
#' @return Une \code{\link[stats]{formula}} mise à jour, incluant les interactions pures
#'         détectées comme nécessaires.
#'
#' @details
#' \enumerate{
#'   \item Les prédicteurs (variables explicatives) sont extraits de \code{formula}
#'         et vérifiés dans \code{g}.
#'   \item Pour chaque paire de facteurs, un G-test est réalisé via \code{\link[DescTools]{GTest}}.
#'   \item Les p-values obtenues sont corrigées par la méthode de Holm.
#'   \item Si une p-value corrigée est inférieure à \code{alpha}, l'interaction pure
#'         entre ces deux facteurs est ajoutée à la formule.
#' }
#'
#' @importFrom DescTools GTest
#'
#' @examples
#' \dontrun{
#' # Exemple d'utilisation
#' data <- data.frame(
#'   Note = rnorm(100),
#'   Lycée = factor(sample(c("LycéeA", "LycéeB"), 100, replace = TRUE)),
#'   Sexe = factor(sample(c("M", "F"), 100, replace = TRUE)),
#'   Options = factor(sample(c("Musique", "Arts"), 100, replace = TRUE))
#' )
#'
#' # Formule de base
#' f <- Note ~ Lycée * Sexe + Options
#'
#' # Mise à jour de la formule avec les interactions nécessaires
#' new_f <- control_independence(f, data, alpha = 0.05, debug = TRUE)
#' new_f
#' }
#'
#' @export
control_independence <- function(formula, g, alpha = 0.05, debug = FALSE, verbose=FALSE, k=NULL) {
  # Vérifier la dépendance entre facteurs de 'g' présents dans 'formula'
  # et ajouter les interactions pures (X:Y) pour les paires détectées dépendantes,
  # sans retirer les interactions déjà présentes (X*Y).

  # 1) Récupération de la variable réponse et des prédicteurs (variables réelles)
  all_vars_in_formula <- all.vars(formula)
  response_var <- all_vars_in_formula[1]                   # e.g. "Note"
  explanatory_vars <- setdiff(all_vars_in_formula, response_var)

  # Vérifier que ces prédicteurs existent dans g
  if (!all(explanatory_vars %in% names(g))) {
    exit("Explanatory variables in formula not found in 'g'.\n","Les variables explicatives de la formule ne correspondent pas aux colonnes de g.\n")
  }

  # data.frame des prédicteurs
  g_subset <- g[, explanatory_vars, drop = FALSE]

  # Si verbose == TRUE, on annonce la méthode
  if (verbose) {
    cat(
      paste0(
        k,
        .msg(
          ") Method: Crossed GTests from {DescTools} to detect factor dependencies (Holm correction).\n",
          ") Méthode : GTests croisés du package {DescTools} pour détecter les dépendances entre facteurs (correction de Holm).\n"
        )
      )
    )
  }
  # 2) Calcul GTest paire par paire (uniquement si c'est factor vs factor)
  n <- ncol(g_subset)
  p_values <- matrix(NA, n, n, dimnames = list(explanatory_vars, explanatory_vars))

  p_list <- numeric(0)
  pair_positions <- list()

  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      if (is.factor(g_subset[[i]]) && is.factor(g_subset[[j]])) {
        # Effectuer le GTest
        test_result <- DescTools::GTest(g_subset[[i]], g_subset[[j]])
        p_list <- c(p_list, test_result$p.value)
        pair_positions[[length(pair_positions)+1]] <- c(i, j)
      }
    }
  }

  if (length(p_list) == 0) {
    if (debug) {cat(.msg("Debug: No pair of factors detected => formula unchanged.\n","Debug : Aucune paire de facteurs détectée => formule inchangée.\n"))}
    return(formula)
  }

  # 3) Correction de Holm sur les p-values brutes
  p_list_corrected <- p.adjust(p_list, method = "holm")

  # Remplir la matrice p_values corrigée
  idx <- 1
  for (pos in pair_positions) {
    i <- pos[1]
    j <- pos[2]
    p_values[i, j] <- p_list_corrected[idx]
    idx <- idx + 1
  }

  if (debug) {
    cat(.msg("Debug: Corrected p-values (Holm) matrix:\n","Debug : P-values corrigées (matrice) :\n"))
    print(p_values)
  }

  # 4) Déterminer les interactions pur `Var_i:Var_j` à ajouter
  new_inters <- character(0)
  for (i in seq_len(n - 1)) {
    for (j in seq(i+1, n)) {
      if (!is.na(p_values[i,j]) && p_values[i,j] < alpha) {
        # Interaction pure i:j
        inter_term <- paste0(explanatory_vars[i], ":", explanatory_vars[j])
        new_inters <- c(new_inters, inter_term)
      }
    }
  }

  if (length(new_inters) == 0) {
    if (debug) cat(.msg("Debug: No significant dependencies found => no new interactions.\n","Debug : Aucune dépendance détectée entre les facteurs.\n"))
	if (isTRUE(verbose)) {
	  cat(.msg("No significant factor dependencies were found via GTest() => no new interactions added.\n",
	           "Aucune dépendance significative entre facteurs n'a été détectée via GTest() => aucune interaction ajoutée.\n"))
	}
    return(formula)
  }

  # 5) Conserver la partie gauche (réponse) et la partie droite (déjà existante)
  lhs <- deparse(formula[[2]])                     # "Note"
  rhs <- deparse(formula[[3]])                     # e.g. "établissement + Lycée*Sexe + Options + Sport"

  if (debug) {
    cat(.msg("Debug: Original formula RHS:\n","Debug : Partie droite (formule initiale) :\n"), rhs, "\n")
  }

  # Pour chaque nouvelle interaction Var_i:Var_j, on l'ajoute si elle n'est pas déjà présente dans la RHS
  # (ainsi on ne perd pas Lycée*Sexe et on n'écrase rien)
  for (term in unique(new_inters)) {
    # Vérifier si term ("Sexe:Options") n'apparaît déjà pas dans la RHS
    if (!grepl(term, rhs, fixed = TRUE)) {
      # On ajoute
      rhs <- paste(rhs, "+", term)
    }
  }

  # 6) Construire la nouvelle formule
  formula_str <- paste0(lhs, " ~ ", rhs)
  updated_formula <- as.formula(formula_str)

  if (debug) {
    cat(.msg("Debug: Dependencies found => interactions added. Updated formula:\n","Debug : Dépendance détectée => ajout d'interactions. Nouvelle formule :\n"))
    print(updated_formula)
  }

  # Si de nouvelles interactions ont été détectées, on affiche un message (en une phrase).
	if (isTRUE(verbose)) {
	  cat(
		sprintf(
		  .msg("The following interactions were detected via GTest() and added to the model:\n","Les interactions suivantes ont été détectées via GTest() et ajoutées au modèle :\n"),
		  paste(new_inters, collapse = ", ")
		)
	  )
	}
  return(updated_formula)
}

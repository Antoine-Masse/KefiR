#' Filtrage des groupes avec effectifs insuffisants
#'
#' @description
#' Prétraitement qui détecte et filtre les groupes (ou combinaisons de facteurs)
#' ayant moins de `min_n` observations. Retourne les données filtrées et un message
#' informatif sur les groupes écartés.
#'
#' @param data Data frame contenant les données.
#' @param factor_vars Character vector. Noms des variables factorielles.
#' @param min_n Integer. Nombre minimum d'observations par groupe (défaut: 3).
#' @param verbose Logical. Afficher les messages (défaut: TRUE).
#' @param k Integer. Compteur de messages pour .vbse().
#' @param code Logical. Mode code.
#'
#' @return Liste contenant:
#' \itemize{
#'   \item \code{data}: Data frame filtré (sans les groupes insuffisants).
#'   \item \code{removed_groups}: Character vector des groupes écartés.
#'   \item \code{removed_n}: Integer. Nombre total d'observations supprimées.
#'   \item \code{original_n}: Integer. Nombre d'observations avant filtrage.
#'   \item \code{filtered}: Logical. TRUE si des données ont été filtrées.
#'   \item \code{k}: Compteur mis à jour.
#' }
#'
#' @details
#' Cette fonction est appelée en amont de l'analyse pour éviter les erreurs
#' ou avertissements liés aux groupes trop petits (bartlett.test, bootstrap, etc.).
#'
#' Pour les analyses avec plusieurs facteurs, c'est l'interaction croisée
#' de tous les facteurs qui est utilisée pour déterminer les effectifs.
#'
#' @examples
#' \dontrun{
#' # Données avec un groupe trop petit
#' dt <- data.frame(
#'   y = rnorm(20),
#'   g = factor(c(rep("A", 10), rep("B", 8), rep("C", 2)))
#' )
#'
#' # Filtrer groupes < 3 observations
#' result <- .filter_small_groups(dt, "g", min_n = 3)
#' # Le groupe C sera écarté
#' }
#'
#' @keywords internal
#' @export
.filter_small_groups <- function(data,
                                  factor_vars,
                                  min_n = 3,
                                  verbose = TRUE,
                                  k = 0,
                                  code = FALSE) {

  # Validation des entrées
  if (is.null(data) || nrow(data) == 0) {
    return(list(
      data = data,
      removed_groups = character(0),
      removed_n = 0,
      original_n = 0,
      filtered = FALSE,
      k = k
    ))
  }

  if (length(factor_vars) == 0) {
    return(list(
      data = data,
      removed_groups = character(0),
      removed_n = 0,
      original_n = nrow(data),
      filtered = FALSE,
      k = k
    ))
  }

  original_n <- nrow(data)


  # Créer le groupement (interaction si plusieurs facteurs)
  if (length(factor_vars) == 1) {
    group_var <- factor(data[[factor_vars[1]]])
  } else {
    # Interaction de tous les facteurs
    group_var <- interaction(data[, factor_vars, drop = FALSE], drop = TRUE, sep = ":")
  }

  # Calculer effectifs par groupe
  group_counts <- table(group_var)

  # Identifier groupes avec effectifs insuffisants
  small_groups <- names(group_counts)[group_counts < min_n]

  # Si aucun groupe insuffisant, retourner données inchangées

  if (length(small_groups) == 0) {
    return(list(
      data = data,
      removed_groups = character(0),
      removed_n = 0,
      original_n = original_n,
      filtered = FALSE,
      k = k
    ))
  }

  # Calculer nombre d'observations à supprimer
  removed_n <- sum(group_counts[small_groups])

  # Construire message informatif
  # Format: "groupe1 (n=X), groupe2 (n=Y)"
  groups_info <- sapply(small_groups, function(g) {
    paste0(g, " (n=", group_counts[g], ")")
  })

  # Message d'avertissement
  k <- .vbse(
    paste0("WARNING: Groups with insufficient observations detected (n < ", min_n, ").\n",
           "\tGroups removed from analysis: ", paste(groups_info, collapse = ", "), "\n",
           "\tObservations removed: ", removed_n, " / ", original_n, "\n",
           "\tAnalysis will continue with the remaining ", original_n - removed_n, " observations."),
    paste0("ATTENTION : Groupes avec effectifs insuffisants détectés (n < ", min_n, ").\n",
           "\tGroupes écartés de l'analyse : ", paste(groups_info, collapse = ", "), "\n",
           "\tObservations supprimées : ", removed_n, " / ", original_n, "\n",
           "\tL'analyse continuera avec les ", original_n - removed_n, " observations restantes."),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  # Filtrer les données
  keep_idx <- !(group_var %in% small_groups)
  data_filtered <- data[keep_idx, , drop = FALSE]

  # Recalculer les niveaux des facteurs pour supprimer les niveaux vides

  for (fvar in factor_vars) {
    if (is.factor(data_filtered[[fvar]])) {
      data_filtered[[fvar]] <- droplevels(data_filtered[[fvar]])
    }
  }

  return(list(
    data = data_filtered,
    removed_groups = small_groups,
    removed_n = as.integer(removed_n),
    original_n = original_n,
    filtered = TRUE,
    k = k
  ))
}

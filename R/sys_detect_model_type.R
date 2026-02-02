#' Détection du type de modèle (ANOVA vs. ANCOVA) - VERSION 3
#'
#' Cette fonction \code{.detect_model_type} détermine automatiquement, à partir d'une formule
#' R et d'un jeu de données, s'il s'agit d'un scénario d'ANOVA (uniquement des prédicteurs
#' catégoriques) ou d'ANCOVA (combinaison de prédicteurs catégoriques et numériques).
#'
#' Concrètement, la fonction détecte et convertit au besoin les variables binaires (0/1)
#' en facteurs, puis sépare les prédicteurs en deux catégories (\emph{factorielles} et
#' \emph{numériques}). Sur cette base, elle renvoie :
#' \itemize{
#'   \item \code{TRUE} si le modèle correspond à une ANCOVA (facteurs + covariables continues) ;
#'   \item \code{FALSE} si le modèle correspond à une ANOVA (exclusivement des facteurs) ;
#'   \item \code{exit(...)} (arrêt) si aucun facteur catégoriel ou aucune covariable
#'         continue n'est détectée.
#' }
#'
#' @param formula Une \code{\link[stats]{formula}}. Par exemple \code{y ~ x1 + x2 + ...}.
#' @param data Un \code{data.frame} contenant les variables mentionnées dans \code{formula}.
#' @param debug Logique. Si \code{TRUE}, affiche des messages de débogage détaillés (défaut: \code{FALSE}).
#'
#' @details
#' Les variables binaires (\code{0/1}) présentes dans \code{data} sont automatiquement
#' converties en facteurs avant l'identification du type de modèle. Si vous utilisez
#' un paramètre \code{debug} ou des fonctions internes (\code{.msg}, \code{.exit}, \code{.dbg}),
#' assurez-vous qu'ils soient définis dans l'environnement de travail.
#'
#' \strong{VERSION 3 (CHANGEMENT MAJEUR):} Cette fonction retourne maintenant un booléen
#' (\code{TRUE} ou \code{FALSE}) au lieu d'une chaîne de caractères (\code{"TRUE"} ou \code{"FALSE"}).
#' Ceci améliore la cohérence du code et évite les comparaisons de chaînes.
#'
#' @return
#' Un booléen logique :
#' \itemize{
#'   \item \code{TRUE} si le modèle identifié est de type ANCOVA ;
#'   \item \code{FALSE} si le modèle identifié est de type ANOVA ;
#'   \item Arrêt de la fonction avec un message d'erreur si aucun facteur catégoriel
#'         ou covariable continue n'est détecté.
#' }
#'
#' @examples
#' \dontrun{
#' # Exemple de données
#' df <- data.frame(
#'   y = rnorm(100),
#'   groupe = factor(rep(LETTERS[1:2], 50)),
#'   x = runif(100),
#'   binaire = rbinom(100, 1, 0.5)
#' )
#'
#' # Cas 1 : formule avec un facteur et une covariable continue => TRUE (ANCOVA)
#' .detect_model_type(y ~ groupe + x, data = df)
#'
#' # Cas 2 : uniquement un facteur => FALSE (ANOVA)
#' .detect_model_type(y ~ groupe, data = df)
#'
#' # Cas 3 : uniquement une variable binaire (convertie en facteur) mais pas de covariable
#' # => Ici, la fonction renvoie FALSE (ANOVA) si la variable binaire est considérée comme un facteur
#'
#' # Cas 4 : si aucune variable factorielle ni numérique n'est trouvée => arrêt
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{terms}} pour l'extraction des prédicteurs d'une formule,
#'   \item \code{\link{aov}} et \code{\link{lm}} pour l'ajustement des modèles ANOVA ou ANCOVA.
#' }
#'
#' @keywords internal
#' @export
.detect_model_type <- function(formula, data, debug=FALSE) {
  # Extraire les noms des prédicteurs du modèle
  predictors <- attr(terms(formula), "term.labels")

  # DEBUG : Afficher les informations
  if (debug) {
    cat("\n=== DEBUG .detect_model_type() ===\n")
    cat("Formule :\n")
    print(formula)
    cat("\nPrédicteurs extraits de la formule :", paste(predictors, collapse=", "), "\n")
    cat("Noms des colonnes dans data :", paste(names(data), collapse=", "), "\n")
    cat("class(data) :", class(data), "\n")
    cat("is.data.frame(data) :", is.data.frame(data), "\n")
    cat("ncol(data) :", ncol(data), "\n")
    cat("nrow(data) :", nrow(data), "\n")
  }
  # CORRECTION : Filtrer pour ne garder que les termes simples (sans interaction)
  # et qui existent dans data
  simple_predictors <- predictors[!grepl(":", predictors)]

  if (debug) {
    cat("Prédicteurs simples (sans interaction) :", paste(simple_predictors, collapse=", "), "\n")
  }

  simple_predictors <- simple_predictors[simple_predictors %in% names(data)]

  if (debug) {
    cat("Prédicteurs présents dans data :", paste(simple_predictors, collapse=", "), "\n")
    cat("=====================================\n\n")
  }

  # Si aucun prédicteur valide, erreur
  if (length(simple_predictors) == 0) {
    .exit(
      "No valid predictor found in data.",
      "Aucun prédicteur valide trouvé dans les données."
    )
  }

  # Vérifier si les prédicteurs sont binaires (0/1) et convertir en facteur si nécessaire
  binary_as_factor <- function(column) {
    if (is.numeric(column) && length(unique(column)) == 2 && all(sort(unique(column)) == c(0, 1))) {
      return(as.factor(column))
    }
    return(column)
  }

  # Appliquer la conversion sur les colonnes des prédicteurs
  data <- data.frame(lapply(data[simple_predictors], binary_as_factor))

  # Identifier les facteurs et les covariables continues
  factor_vars <- simple_predictors[sapply(data, is.factor)]
  numeric_vars <- simple_predictors[sapply(data, is.numeric)]

  # Décider du type de modèle
  # VERSION 3: Retourne maintenant un booléen (TRUE/FALSE) au lieu d'une chaîne ("TRUE"/"FALSE")
  if (length(factor_vars) > 0 && length(numeric_vars) > 0) {
    .dbg("", "Scénario détecté : ANCOVA (facteurs catégoriques et covariables continues).", debug=debug)
    return("TRUE")
  } else if (length(factor_vars) > 0) {
    .dbg("", "Scénario détecté : ANOVA (uniquement des facteurs catégoriques).", debug=debug)
    return("FALSE")
  } else {
    .exit(
      "Model not suitable: no categorical factors or continuous covariates found.",
      "Modèle non adapté : il n'y a ni facteurs catégoriques ni covariables continues."
    )
  }
}

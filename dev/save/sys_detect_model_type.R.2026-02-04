#' Detection du type de modele (ANOVA vs. ANCOVA) - VERSION 3
#'
#' Cette fonction \code{.detect_model_type} determine automatiquement, a partir d'une formule
#' R et d'un jeu de donnees, s'il s'agit d'un scenario d'ANOVA (uniquement des predicteurs
#' categoriques) ou d'ANCOVA (combinaison de predicteurs categoriques et numeriques).
#'
#' Concretement, la fonction detecte et convertit au besoin les variables binaires (0/1)
#' en facteurs, puis separe les predicteurs en deux categories (\emph{factorielles} et
#' \emph{numeriques}). Sur cette base, elle renvoie :
#' \itemize{
#'   \item \code{TRUE} si le modele correspond a une ANCOVA (facteurs + covariables continues) ;
#'   \item \code{FALSE} si le modele correspond a une ANOVA (exclusivement des facteurs) ;
#'   \item \code{exit(...)} (arret) si aucun facteur categoriel ou aucune covariable
#'         continue n'est detectee.
#' }
#'
#' @param formula Une \code{\link[stats]{formula}}. Par exemple \code{y ~ x1 + x2 + ...}.
#' @param data Un \code{data.frame} contenant les variables mentionnees dans \code{formula}.
#' @param debug Logique. Si \code{TRUE}, affiche des messages de debogage detailles (defaut: \code{FALSE}).
#'
#' @details
#' Les variables binaires (\code{0/1}) presentes dans \code{data} sont automatiquement
#' converties en facteurs avant l'identification du type de modele. Si vous utilisez
#' un parametre \code{debug} ou des fonctions internes (\code{.msg}, \code{.exit}, \code{.dbg}),
#' assurez-vous qu'ils soient definis dans l'environnement de travail.
#'
#' \strong{VERSION 3 (CHANGEMENT MAJEUR):} Cette fonction retourne maintenant un booleen
#' (\code{TRUE} ou \code{FALSE}) au lieu d'une chaine de caracteres (\code{"TRUE"} ou \code{"FALSE"}).
#' Ceci ameliore la coherence du code et evite les comparaisons de chaines.
#'
#' @return
#' Un booleen logique :
#' \itemize{
#'   \item \code{TRUE} si le modele identifie est de type ANCOVA ;
#'   \item \code{FALSE} si le modele identifie est de type ANOVA ;
#'   \item Arret de la fonction avec un message d'erreur si aucun facteur categoriel
#'         ou covariable continue n'est detecte.
#' }
#'
#' @examples
#' \dontrun{
#' # Exemple de donnees
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
#' # => Ici, la fonction renvoie FALSE (ANOVA) si la variable binaire est consideree comme un facteur
#'
#' # Cas 4 : si aucune variable factorielle ni numerique n'est trouvee => arret
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{terms}} pour l'extraction des predicteurs d'une formule,
#'   \item \code{\link{aov}} et \code{\link{lm}} pour l'ajustement des modeles ANOVA ou ANCOVA.
#' }
#'
#' @keywords internal
.detect_model_type <- function(formula, data, debug=FALSE) {
  # Extraire les noms des predicteurs du modele
  predictors <- attr(terms(formula), "term.labels")

  # DEBUG : Afficher les informations
  if (debug) {
    cat("\n=== DEBUG .detect_model_type() ===\n")
    cat("Formule :\n")
    print(formula)
    cat("\nPr\u00e9dicteurs extraits de la formule :", paste(predictors, collapse=", "), "\n")
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
    cat("Pr\u00e9dicteurs simples (sans interaction) :", paste(simple_predictors, collapse=", "), "\n")
  }

  simple_predictors <- simple_predictors[simple_predictors %in% names(data)]

  if (debug) {
    cat("Pr\u00e9dicteurs pr\u00e9sents dans data :", paste(simple_predictors, collapse=", "), "\n")
    cat("=====================================\n\n")
  }

  # Si aucun predicteur valide, erreur
  if (length(simple_predictors) == 0) {
    .exit(
      "No valid predictor found in data.",
      "Aucun pr\u00e9dicteur valide trouv\u00e9 dans les donn\u00e9es."
    )
  }

  # Verifier si les predicteurs sont binaires (0/1) et convertir en facteur si necessaire
  binary_as_factor <- function(column) {
    if (is.numeric(column) && length(unique(column)) == 2 && all(sort(unique(column)) == c(0, 1))) {
      return(as.factor(column))
    }
    return(column)
  }

  # Appliquer la conversion sur les colonnes des predicteurs
  data <- data.frame(lapply(data[simple_predictors], binary_as_factor))

  # Identifier les facteurs et les covariables continues
  factor_vars <- simple_predictors[sapply(data, is.factor)]
  numeric_vars <- simple_predictors[sapply(data, is.numeric)]

  # Decider du type de modele
  # VERSION 3: Retourne maintenant un booleen (TRUE/FALSE) au lieu d'une chaine ("TRUE"/"FALSE")
  if (length(factor_vars) > 0 && length(numeric_vars) > 0) {
    .dbg("", "Sc\u00e9nario d\u00e9tect\u00e9 : ANCOVA (facteurs cat\u00e9goriques et covariables continues).", debug=debug)
    return("TRUE")
  } else if (length(factor_vars) > 0) {
    .dbg("", "Sc\u00e9nario d\u00e9tect\u00e9 : ANOVA (uniquement des facteurs cat\u00e9goriques).", debug=debug)
    return("FALSE")
  } else {
    .exit(
      "Model not suitable: no categorical factors or continuous covariates found.",
      "Mod\u00e8le non adapt\u00e9 : il n'y a ni facteurs cat\u00e9goriques ni covariables continues."
    )
  }
}

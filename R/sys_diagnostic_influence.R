#' Diagnostic d'influence pour modèles linéaires (Cook, Leverage, DFBETAS)
#'
#' @description
#' Fonction helper pour diagnostiquer les observations influentes dans un modèle
#' linéaire (lm, ANOVA, ANCOVA). Calcule Cook's distance, leverage (hat values)
#' et DFBETAS pour identifier les points qui ont un impact disproportionné sur
#' le modèle ajusté.
#'
#' @param model Modèle lm() ajusté
#' @param data Data frame utilisé pour ajuster le modèle
#' @param alpha Seuil de significativité (défaut: 0.05) - non utilisé actuellement
#' @param verbose Logical. Afficher messages détaillés (défaut: FALSE)
#' @param k Compteur d'étapes (défaut: 0)
#' @param debug Logical. Mode debug (défaut: FALSE)
#'
#' @return Liste contenant:
#' \itemize{
#'   \item \code{cook_d}: Vecteur des distances de Cook
#'   \item \code{leverage}: Vecteur des leverages (hat values)
#'   \item \code{dfbetas}: Matrice des DFBETAS (n × p)
#'   \item \code{thresholds}: Liste des seuils (cook, leverage, dfbetas)
#'   \item \code{influential_cook}: Indices observations Cook > seuil
#'   \item \code{influential_leverage}: Indices observations leverage > seuil
#'   \item \code{influential_dfbetas}: Indices observations DFBETAS > seuil (any coef)
#'   \item \code{influential_any}: Union de tous les influents
#'   \item \code{critical}: Indices observations Cook > 1 (influence critique)
#'   \item \code{n_influential}: Nombre total d'observations influentes
#'   \item \code{n_critical}: Nombre observations critiques
#'   \item \code{max_cook}: Cook max
#'   \item \code{max_leverage}: Leverage max
#' }
#'
#' @details
#' **Méthodes de diagnostic** :
#'
#' 1. **Cook's Distance** : Mesure l'influence combinée d'une observation sur
#'    TOUS les coefficients β. Combine leverage (distance géométrique) et résidu
#'    standardisé. Seuil conservateur : 4/n ; critique si > 1.
#'
#' 2. **Leverage (hat values)** : Mesure si une observation est extrême dans
#'    l'espace des prédicteurs (covariables, facteurs). Seuil : 2p/n où p est
#'    le nombre de paramètres.
#'
#' 3. **DFBETAS** : Mesure l'impact de chaque observation sur CHAQUE coefficient
#'    individuellement. Seuil : 2/√n. Plus spécifique que Cook.
#'
#' **Interprétation** :
#' - Leverage élevé SEUL : observation extrême en X mais bien ajustée
#' - Résidu élevé SEUL : mauvais ajustement mais faible impact
#' - Cook élevé : combinaison des deux → influence forte sur β
#'
#' @references
#' - Cook, R. D. (1977). Detection of influential observations in linear regression.
#'   Technometrics, 19(1), 15-18.
#' - Fox, J. (2016). Applied Regression Analysis and Generalized Linear Models
#'   (3rd ed.). Sage.
#' - Tabachnick, B. G., & Fidell, L. S. (2013). Using Multivariate Statistics
#'   (6th ed.). Pearson.
#' - Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). Regression Diagnostics.
#'   Wiley.
#'
#' @keywords internal
#' @export
.diagnostic_influence <- function(model, data, alpha = 0.05,
                                  verbose = FALSE, k = 0, debug = FALSE) {

  # ==========================================================================
  # VALIDATION
  # ==========================================================================

  if (!inherits(model, "lm")) {
    stop("Model must be of class 'lm'")
  }

  # ==========================================================================
  # CALCULS DIAGNOSTICS
  # ==========================================================================

  n <- nrow(data)
  p <- length(coef(model))  # Nombre de paramètres (incluant intercept)

  .dbg("Computing influence diagnostics...",
       "Calcul diagnostics d'influence...",
       debug = debug)

  # Cook's distance
  cook_d <- cooks.distance(model)

  # Leverage (hat values)
  leverage <- hatvalues(model)

  # DFBETAS (impact sur chaque coefficient)
  dfbetas_vals <- dfbetas(model)

  # ==========================================================================
  # SEUILS ACADÉMIQUES
  # ==========================================================================

  # Cook : 4/n (conservateur, Fox 2016) ou 1 (critique)
  cook_threshold <- 4/n
  cook_critical <- 1

  # Leverage : 2p/n (standard) ou 3p/n (strict)
  leverage_threshold <- 2*p/n

  # DFBETAS : 2/√n (Belsley et al. 1980)
  dfbetas_threshold <- 2/sqrt(n)

  .dbg(paste0("Thresholds: Cook=", round(cook_threshold, 4),
              ", Leverage=", round(leverage_threshold, 4),
              ", DFBETAS=", round(dfbetas_threshold, 4)),
       paste0("Seuils: Cook=", round(cook_threshold, 4),
              ", Leverage=", round(leverage_threshold, 4),
              ", DFBETAS=", round(dfbetas_threshold, 4)),
       debug = debug)

  # ==========================================================================
  # DÉTECTION OBSERVATIONS INFLUENTES
  # ==========================================================================

  # Cook > seuil
  influential_cook <- which(cook_d > cook_threshold)

  # Leverage > seuil
  influential_leverage <- which(leverage > leverage_threshold)

  # DFBETAS > seuil (au moins 1 coefficient)
  influential_dfbetas <- which(apply(abs(dfbetas_vals) > dfbetas_threshold, 1, any))

  # Union (au moins 1 critère)
  influential_any <- unique(c(influential_cook, influential_leverage, influential_dfbetas))

  # Critique : Cook > 1
  critical_influence <- which(cook_d > cook_critical)

  # Max
  max_cook <- max(cook_d, na.rm = TRUE)
  max_leverage <- max(leverage, na.rm = TRUE)

  .dbg(paste0("Influential: ", length(influential_any),
              " (Cook: ", length(influential_cook),
              ", Leverage: ", length(influential_leverage),
              ", DFBETAS: ", length(influential_dfbetas), ")"),
       paste0("Influents: ", length(influential_any),
              " (Cook: ", length(influential_cook),
              ", Leverage: ", length(influential_leverage),
              ", DFBETAS: ", length(influential_dfbetas), ")"),
       debug = debug)

  # ==========================================================================
  # RETOUR
  # ==========================================================================

  results <- list(
    # Valeurs brutes
    cook_d = cook_d,
    leverage = leverage,
    dfbetas = dfbetas_vals,

    # Seuils
    thresholds = list(
      cook = cook_threshold,
      cook_critical = cook_critical,
      leverage = leverage_threshold,
      dfbetas = dfbetas_threshold
    ),

    # Détection par critère
    influential_cook = influential_cook,
    influential_leverage = influential_leverage,
    influential_dfbetas = influential_dfbetas,

    # Détection combinée
    influential_any = influential_any,
    critical = critical_influence,

    # Compteurs
    n_influential = length(influential_any),
    n_critical = length(critical_influence),

    # Max
    max_cook = max_cook,
    max_leverage = max_leverage
  )

  return(results)
}

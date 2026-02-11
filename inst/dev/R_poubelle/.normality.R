#' Test de normalité adaptatif avec stratégies multiples
#'
#' Évalue la normalité d'une variable numérique avec sélection automatique du test
#' approprié selon la taille des échantillons et l'hétérogénéité des groupes.
#'
#' @param x Vecteur numérique. Variable dont la normalité doit être testée.
#' @param g Vecteur ou facteur de regroupement. Si \code{NULL}, toutes les données
#'   sont traitées comme un seul groupe (défaut : \code{NULL}).
#' @param alpha Numérique. Seuil de signification pour évaluer la normalité
#'   (défaut : \code{0.05}).
#' @param tolerance Chaîne de caractères. Niveau de tolérance pour les critères
#'   de Skewness/Kurtosis :
#'   \itemize{
#'     \item \code{"basic"} : |Skewness| ≤ 1 et |Kurtosis| ≤ 1.5 (Kline, 2011)
#'     \item \code{"extrem"} : |Skewness| ≤ 2 et |Kurtosis| ≤ 7 (Blanca et al., 2017)
#'   }
#'   Défaut : \code{"basic"}.
#' @param paired Logique. Si \code{TRUE}, teste la normalité des différences
#'   pour données appariées (défaut : \code{FALSE}).
#' @param debug Logique. Si \code{TRUE}, affiche des messages de débogage
#'   (défaut : \code{FALSE}).
#' @param verbose Logique. Si \code{TRUE}, affiche des messages détaillés
#'   (défaut : \code{FALSE}).
#' @param k Entier. Compteur de messages (défaut : \code{0}).
#' @param cpt Chaîne de caractères. Mode de comptage des messages
#'   (défaut : \code{"on"}).
#'
#' @return Une liste contenant :
#'   \itemize{
#'     \item \code{check_normality} : Logique. \code{TRUE} si les données sont
#'       considérées normales, \code{FALSE} sinon.
#'     \item \code{k} : Entier. Compteur de messages mis à jour.
#'   }
#'
#' @details
#' \strong{Stratégie de sélection des tests :}
#' \enumerate{
#'   \item \strong{Groupes homogènes (n ≤ 50)} : Test de Shapiro-Wilk
#'   \item \strong{Groupes homogènes (50 < n ≤ 500)} : Test de Jarque-Bera modifié
#'   \item \strong{Grands échantillons (n > 500) ou hétérogénéité forte} :
#'         Critères heuristiques basés sur Skewness et Kurtosis
#' }
#'
#' \strong{Gestion de l'hétérogénéité :}
#' Lorsque les tailles de groupes traversent les seuils critiques (50 ou 500),
#' une stratégie uniforme (heuristique) est appliquée pour éviter les biais
#' liés aux différences de puissance statistique entre tests.
#'
#' \strong{Correction pour comparaisons multiples :}
#' \itemize{
#'   \item 2-10 groupes : Correction de Sidak
#'   \item > 10 groupes : Correction FDR (Benjamini-Hochberg)
#'   \item Heuristique : Pas de correction (valeurs non probabilistes)
#' }
#'
#' @references
#' Blanca, M. J., Alarcón, R., Arnau, J., Bono, R., & Bendayan, R. (2017).
#' Non-normal data: Is ANOVA still a valid option? \emph{Psicothema}, 29(4), 552-557.
#' \url{https://diposit.ub.edu/dspace/bitstream/2445/122126/1/671797.pdf}
#'
#' Chaffin, W. W., & Rhiel, S. G. (1993). The effect of skewness and kurtosis
#' on the one-sample T test and the impact of knowledge of the population
#' standard deviation. \emph{Journal of Statistical Computation and Simulation},
#' 46(1-2), 79-90.
#'
#' Glinskiy, V. V., Zhukova, M. A., & Tsybatov, A. E. (2024).
#' Modifications to the Jarque-Bera Test. \emph{Mathematics}, 12(16), 2523.
#' \doi{10.3390/math12162523}
#'
#' Kline, R. B. (2011). \emph{Principles and practice of structural equation
#' modeling} (4th ed.). Guilford Press.
#'
#' Korkmaz, S., & Demir, Y. (2023). Investigation of some univariate normality
#' tests in terms of type-I errors and test power. \emph{Journal of Scientific
#' Reports-A}, (52), 376-395.
#'
#' @seealso
#' \code{\link{shapiro.test}}, \code{\link{jb.norm.test}},
#' \code{\link[agricolae]{skewness}}, \code{\link[agricolae]{kurtosis}}
#'
#' @examples
#' # Exemple 1 : Un seul groupe (petite taille)
#' x1 <- rnorm(30)
#' .normality(x1)  # Utilise Shapiro-Wilk
#'
#' # Exemple 2 : Plusieurs groupes homogènes
#' x2 <- c(rnorm(40, mean = 5), rnorm(40, mean = 6), rnorm(40, mean = 7))
#' g2 <- factor(rep(c("A", "B", "C"), each = 40))
#' .normality(x2, g2)  # Utilise Shapiro-Wilk avec correction Sidak
#'
#' # Exemple 3 : Grand échantillon
#' x3 <- rnorm(600)
#' .normality(x3)  # Utilise critères Skewness/Kurtosis
#'
#' # Exemple 4 : Groupes hétérogènes (force heuristique)
#' x4 <- c(rnorm(30), rnorm(200))
#' g4 <- factor(rep(c("Small", "Large"), c(30, 200)))
#' .normality(x4, g4)  # Stratégie uniforme appliquée
#'
#' # Exemple 5 : Tolérance extrême pour données réelles
#' x5 <- rgamma(100, shape = 2, rate = 1)  # Légèrement asymétrique
#' .normality(x5, tolerance = "extrem")
#'
#' @keywords internal
#' @export
.normality <- function(x, g = NULL, alpha=0.05, tolerance="basic", paired=FALSE,
                       debug = FALSE, verbose=FALSE, k=0, cpt="on") {

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 1 : INITIALISATION ET VALIDATION
  # ═══════════════════════════════════════════════════════════════════════════
  # Korkmaz, S., & Demir, Y. (2023). Investigation of some univariate normality tests
  # in terms of type-I errors and test power. Journal of Scientific Reports-A, (52), 376–395.
  # https://dergipark.org.tr/en/download/article-file/2847569

  .dbg("Start of normality assessment.",
       "Début de l'évaluation de la normalité.", debug = debug)

  if (is.null(g)) {
    g <- factor(rep("A", length(x)))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 2 : ANALYSE DES TAILLES DE GROUPES
  # ═══════════════════════════════════════════════════════════════════════════
  # Déterminer la taille du plus petit groupe pour uniformiser le choix du test
  group_sizes <- table(g)
  min_size <- min(group_sizes)
  max_size <- max(group_sizes)

  if (length(unique(g))==1) {
    multi <- FALSE
    # Limites pour un seul groupe (contrôles des résidus d'une ANOVA)
  } else {
    multi <- TRUE
  }

  # Détecter hétérogénéité forte des tailles de groupes
  # Si hétérogénéité, la normalité est limitée à des contrôles de skewness et kurtosis
  if ((min_size<=50 & max_size>50) | (min_size<=500 & max_size>500)) {
    reajust <- TRUE
  } else {
    reajust <- FALSE
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 3 : FONCTION DE TEST PAR GROUPE
  # ═══════════════════════════════════════════════════════════════════════════
  subnormality <- function(vector, multi=multi) {

    n <- length(vector)  # Taille du groupe ACTUEL

    #================================
    # Scénario de base
    #================================
    if (tolerance=="basic") {
      if ((n <= 50) & (reajust==FALSE)) {
        return(shapiro.test(as.numeric(vector))$p.value)
      } else if ((n <= 500) & (reajust==FALSE)) {
        # Ce groupe spécifique dépasse la limite de Shapiro
        # Il présente une bonne puissance au-delà de 300
        # Glinskiy, V. V., Zhukova, M. A., & Tsybatov, A. E. (2024). Modifications to the Jarque–Bera Test. Mathematics, 12(16), 2523. https://doi.org/10.3390/math12162523
        return(jb.norm.test(vector)$p.value)
      } else {
        # |Skewness| < 1 et |Kurtosis| < 1.5 sont les seuils recommandés pour considérer une distribution comme "approximativement normale", selon Kline (2011).
        # Kline, R. B. (2011). Principles and practice of structural equation modeling (4th ed.). New York, NY: Guilford Press.
        # Blanca, M.J., Alarcón, R., Arnau, J., Bono, R., & Bendayan, R. (2017). Non‑normal data: Is ANOVA still a valid option? Psicothema, 29(4), 552‑557. DOI 10.7334/psicothema2016.383
        mysk <- abs(skewness(vector))
        myku <- abs(kurtosis(vector))
        if ((mysk<=1) & (myku<=1.5)) {
          return(1)
        } else {
          return(0)
        }
      }
    } else if (tolerance=="extrem") {
      mysk <- abs(skewness(vector))
      myku <- abs(kurtosis(vector))
      if (multi==TRUE) { # Contrôle des groupes Chaffin et al.
        # Chaffin, W. W., & Rhiel, S. G. (1993). The effect of skewness and kurtosis on the one-sample T test and the impact of knowledge of the population standard deviation. Journal of Statistical Computation and Simulation, 46(1-2), 79-90.
        if ((mysk<=1) & (myku<=4.5) & (n>=20)) {
          return(1)
        } else if ((mysk<=1.5) & (myku<=5) & (n>=30)) {
          return(1)
        } else if ((mysk<=2) & (myku<=6.5) & (n>=50)) {
          return(1)
        } else {
          return(0)
        }
      } else if ((multi==FALSE) & (mysk<=2) & (myku<=7) & (n>=15)) { # Contrôle des résidus Blanca et al.
        return(1)
      } else {
        return(0)
      }
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 4 : APPLICATION AUX GROUPES
  # ═══════════════════════════════════════════════════════════════════════════
  pvals <- by(x, g, subnormality, multi=multi)

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 5 : CORRECTION POUR COMPARAISONS MULTIPLES
  # ═══════════════════════════════════════════════════════════════════════════
  if ((reajust==TRUE) | (min_size>500) | (tolerance=="extrem")) {
    # Aucune correction des valeurs p qui n'en sont pas.
    check_normality <- min(pvals) >= alpha
  } else {
    # Sidak devient trop conservatif au-delà de 10 groupes...
    # Bender, R., & Lange, S. (2001). Adjusting for multiple testing—when and how?  DOI: 10.1016/S0895-4356(00)00314-0
    # Méthode FDR la correction par FDR (False Discovery Rate) – notamment via la méthode Benjamini-Hochberg (BH)
    # Murray, E. J., Berrie, L., & Matthews, J. N. S. (2021). Understanding multiplicity in clinical trials: the impact of multiple endpoints in the interpretation of treatment effects.
    if (length(unique(g))>1 && length(unique(g))<=10) {
      pval <- 1-(1-alpha)^(1/length(unique(g)))
      check_normality <- min(pvals) >= pval
    } else if (length(unique(g))>10) {
      pvals_fdr <- p.adjust(pvals, method = "BH")
      check_normality <- min(pvals_fdr) >= alpha
    } else {
      # Pas de correction (1 groupe)
      check_normality <- min(pvals) >= alpha
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 6 : MESSAGES VERBOSE
  # ═══════════════════════════════════════════════════════════════════════════
  if (verbose==TRUE) {

    if (paired==FALSE) {

      if ((multi==TRUE) & (tolerance=="basic")) {
        # Affichage du bilan pour plusieurs groupes
        k <- .vbse(
          paste0("Normality check - Test selection based on group size.\n\t",
                 "Shapiro-Wilk (≤50), Jarque-Bera (≤500), or skewness/kurtosis (>500).\n\t",
                 if ((min_size <= 50 && max_size > 50) || (min_size <= 500 && max_size > 500)) {
                   "Note: Due to data imbalance, all groups are checked using skewness/kurtosis.\n\t"
                 } else {
                   if (max_size > 500) {
                     "Note: Due to large group sizes, all groups are checked using skewness/kurtosis.\n\t"
                   } else {
                     paste0("Note: p-values are compared to a Sidak-corrected alpha of ", .format_pval(pval), ".\n\t")
                   }
                 },
                 if (check_normality==TRUE) {
                   paste0("==> All groups are normal.",
                          if ((min_size <= 50 && max_size > 50) || (max_size > 500)) {
                            "\n\t...with |skewness| ≤ 1 and |kurtosis| ≤ 1.5."
                          } else if (min(pvals) != 1) {
                            paste0(" (min p-value: ", .format_pval(min(pvals)), ")")
                          })
                 } else {
                   paste0("==> At least one group is non-normal.",
                          if ((min_size <= 50 && max_size > 50) || (max_size > 500)) {
                            "\n\t...with |skewness| > 1 or |kurtosis| > 1.5."
                          } else if (min(pvals) != 0) {
                            paste0(" (min p-value: ", .format_pval(min(pvals)), ")")
                          })
                 }),
          paste0("Contrôle de normalité - Sélection du test selon la taille des groupes.\n\t",
                 "Shapiro-Wilk (≤50), Jarque-Bera (≤500), ou skewness/kurtosis (>500).\n\t",
                 if ((min_size <= 50 && max_size > 50) || (min_size <= 500 && max_size > 500)) {
                   "Note : du fait du déséquilibre des données, les groupes sont tous contrôlés par skewness/kurtosis.\n\t"
                 } else {
                   if (max_size > 500) {
                     "Note : du fait de la taille importante des groupes, ceux-ci sont tous contrôlés par skewness/kurtosis.\n\t"
                   } else {
                     paste0("Note : les p-value sont comparées à un alpha corrigé (Sidak) de ", .format_pval(pval), ".\n\t")
                   }
                 },
                 if (check_normality==TRUE) {
                   paste0("==> Tous les groupes sont normaux.",
                          if ((min_size <= 50 && max_size > 50) || (max_size > 500)) {
                            "\n\t...avec un |skewness| ≤ 1 et un |kurtosis| ≤ 1.5."
                          } else if (min(pvals) != 1) {
                            paste0(" (p-value min : ", .format_pval(min(pvals)), ")")
                          })
                 } else {
                   paste0("==> Au moins un des groupes est non-normal.",
                          if ((min_size <= 50 && max_size > 50) || (max_size > 500)) {
                            "\n\t...avec un |skewness| > 1 ou un |kurtosis| > 1.5."
                          } else if (min(pvals) != 0) {
                            paste0(" (p-value min : ", .format_pval(min(pvals)), ")")
                          })
                 }),
          verbose = verbose, k = k, cpt = cpt)

      } else if ((multi==FALSE) & (tolerance=="basic")) {
        # Contrôle des résidus
        n_g <- length(x)
        k <- .vbse(
          paste0("Academic check of residual normality.\n\t",
                 "Shapiro-Wilk (≤50), Jarque-Bera (≤500), or skewness/kurtosis (>500).\n\t",
                 if (check_normality==TRUE) {
                   paste0("==> Residuals are normal.",
                          if (n_g > 500) {
                            "\n\t...with |skewness| ≤ 1 and |kurtosis| ≤ 1.5.\n\t\tWe maintain the choice of performing ANOVA."
                          } else if (min(pvals) != 1) {
                            paste0(" (p-value: ", .format_pval(min(pvals)), ")")
                          })
                 } else {
                   paste0("==> Residuals are non-normal.",
                          if (n_g > 500) {
                            "\n\t...with |skewness| > 1 or |kurtosis| > 1.5."
                          } else if (min(pvals) != 0) {
                            paste0(" (p-value: ", .format_pval(min(pvals)), ")")
                          })
                 }),
          paste0("Contrôle ACADEMIQUE de la normalité des résidus.\n\t",
                 "Shapiro-Wilk (≤50), Jarque-Bera (≤500), ou skewness/kurtosis (>500).\n\t",
                 if (check_normality==TRUE) {
                   paste0("==> Les résidus sont normaux.",
                          if (n_g > 500) {
                            "\n\t...avec un |skewness| ≤ 1 et un |kurtosis| ≤ 1.5.\n\t\tOn reste sur le choix de faire une ANOVA."
                          } else if (min(pvals) != 1) {
                            paste0(" (p-value : ", .format_pval(min(pvals)), ")")
                          })
                 } else {
                   paste0("==> Les résidus sont non-normaux.",
                          if (n_g > 500) {
                            "\n\t...avec un |skewness| > 1 ou un |kurtosis| > 1.5."
                          } else if (min(pvals) != 0) {
                            paste0(" (p-value : ", .format_pval(min(pvals)), ")")
                          })
                 }),
          verbose = verbose, k = k, cpt = cpt)

      } else if ((multi==TRUE) & (tolerance=="extrem")) {
        # Affichage du bilan avec tolérance extrême
        k <- .vbse(
          paste0("Normality check with variable tolerance based on group size.\n\t",
                 "Skewness/Kurtosis ≤1/4.5 (n≥20); ≤1.5/5 (n≥30); ≤2/6.5 (n≥50) (Chaffin et al. 1993).\n\t",
                 if (check_normality==TRUE) {
                   "==> Groups show acceptable normality."
                 } else {
                   "==> At least one group shows extreme non-normality."
                 }),
          paste0("Contrôle de la normalité avec tolérance variable selon la taille des groupes.\n\t",
                 "Skewness/Kurtosis ≤1/4.5 (n≥20) ; ≤1.5/5 (n≥30) ; ≤2/6.5 (n≥50) (Chaffin et al. 1993).\n\t",
                 if (check_normality==TRUE) {
                   "==> Les groupes présentent une normalité acceptable."
                 } else {
                   "==> Au moins un des groupes présente une non-normalité trop extrême."
                 }),
          verbose = verbose, k = k, cpt = cpt)

      } else if ((multi==FALSE) & (tolerance=="extrem")) {
        k <- .vbse(
          paste0("Residual normality check with tolerance if n>15.\n\t",
                 "Skewness ≤2 and Kurtosis ≤7 (Blanca et al. 2018).\n\t",
                 if (check_normality==TRUE) {
                   "==> Residuals show acceptable normality."
                 } else {
                   "==> Residuals show extreme non-normality."
                 }),
          paste0("Contrôle de la normalité des résidus avec tolérance si n>15.\n\t",
                 "Skewness ≤2 et Kurtosis ≤7 (Blanca et al. 2018).\n\t",
                 if (check_normality==TRUE) {
                   "==> Les résidus présentent une normalité acceptable."
                 } else {
                   "==> Les résidus présentent une non-normalité trop extrême."
                 }),
          verbose = verbose, k = k, cpt = cpt)
      }

    } else if (paired==TRUE) {
      # Messages pour données appariées
      if (check_normality == FALSE) {
        k <- .vbse(
          paste0("Normality check of paired differences:\n\tShapiro–Wilk (≤50), Jarque–Bera (≤500), or skewness/kurtosis (>500).\n\t",
                 "The differences are not normally distributed",
                 if (min(pvals) != 0) {
                   paste0(" (p-value: ", .format_pval(min(pvals)), ")")
                 },
                 "."),
          paste0("Contrôle de la normalité des différences appariées :\n\tTest de Shapiro-Wilk (≤50), Jarque-Bera (≤500) ou skewness/kurtosis (>500).\n\t",
                 "Les différences ne sont pas normales",
                 if (min(pvals) != 0) {
                   paste0(" (p-value : ", .format_pval(min(pvals)), ")")
                 },
                 "."),
          verbose = verbose, k = k, cpt=cpt)
      } else {
        k <- .vbse(
          paste0("Normality check of paired differences:\n\tShapiro–Wilk (≤50), Jarque–Bera (≤500), or skewness/kurtosis (>500).\n\t",
                 "The differences are normally distributed",
                 if (min(pvals) != 1) {
                   paste0(" (p-value: ", .format_pval(min(pvals)), ")")
                 },
                 "."),
          paste0("Contrôle de la normalité des différences appariées :\n\tTest de Shapiro-Wilk (≤50), Jarque-Bera (≤500) ou skewness/kurtosis (>500).\n\t",
                 "Les différences sont normales",
                 if (min(pvals) != 1) {
                   paste0(" (p-value : ", .format_pval(min(pvals)), ")")
                 },
                 "."),
          verbose = verbose, k = k, cpt=cpt)
      }
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # SECTION 7 : RETOUR
  # ═══════════════════════════════════════════════════════════════════════════
  normal <- list()
  normal[[1]] <- check_normality
  normal[[2]] <- k
  return(normal)
}

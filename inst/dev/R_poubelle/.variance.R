#' Test d'homogénéité des variances
#'
#' @param x Numeric vector of data values
#' @param g Factor or grouping vector
#' @param check_normality Logical. TRUE if data are normal, FALSE otherwise
#' @param alpha Numeric. Significance threshold (default 0.05)
#' @param paired Logical. TRUE for paired designs, FALSE otherwise
#' @param debug Logical. Debug mode
#' @param verbose Logical. Display messages
#' @param k Integer. Message counter
#' @param code Logical. Display equivalent R code
#'
#' @return A list containing:
#'   \itemize{
#'     \item check_variance: Logical indicating if variances are homogeneous
#'     \item k: Updated message counter
#'   }
#'
#' @importFrom car leveneTest
#' @importFrom onewaytests bf.test
#'
#' @keywords internal
.variance <- function(x, g, check_normality=TRUE, alpha=0.05,
                      paired=FALSE, debug=FALSE, verbose=FALSE, k=0, code=FALSE) {
  ####################################
  # Initialisation
  ####################################
  n_g <- length(unique(g))
  check_variance <- TRUE  # Par défaut

  if (paired==TRUE) {
    #========================================
    # CAS APPARIÉ : Pas de test nécessaire
    #========================================
    k <- .vbse(
      "No variance homogeneity check needed for paired designs.",
      "Il n'est pas nécessaire de contrôler l'homogénéité des variances pour un plan apparié.",
      verbose = verbose, k = k, cpt="on"
    )
    check_variance <- TRUE

  } else if (paired==FALSE) {

    if (n_g == 2) {
      #========================================
      # 2 GROUPES : Test de Fisher-Snedecor
      #========================================
      if (code==TRUE){
        cat("# Test de Fisher-Snedecor\n")
        cat("var.test(x~g)\n")
      }

      fm <- formula(x~g)
      pvals <- var.test(fm)$p.value

      if (pvals > alpha) {
        # VARIANCES HOMOGÈNES
        check_variance <- TRUE
        ang <- paste0("Fisher-Snedecor test [var.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
        fr <- paste0("Test de Fisher-Snedecor [var.test()] - Variances homogènes (p = ", .format_pval(pvals), ").")
        k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="on")
      } else {
        # VARIANCES HÉTÉROGÈNES
        check_variance <- FALSE
        ang <- paste0("Fisher-Snedecor test [var.test()] - Heterogeneous variances (p = ", .format_pval(pvals), ").")
        fr <- paste0("Test de Fisher-Snedecor [var.test()] - Variances hétérogènes (p = ", .format_pval(pvals), ").")
        k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="on")
      }

    } else if (n_g > 2) {

      if (check_normality==TRUE) {
        #========================================
        # DONNÉES NORMALES : Test de Bartlett
        #========================================
        if (code==TRUE){
          cat("# Test de Bartlett\n")
          cat("bartlett.test(x, g)\n")
        }

        pvals <- bartlett.test(x, g)$p.value

        if (pvals > alpha) {
          # VARIANCES HOMOGÈNES
          check_variance <- TRUE
          ang <- paste0("Bartlett test [bartlett.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances homogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="on")
        } else {
          # VARIANCES HÉTÉROGÈNES
          check_variance <- FALSE
          ang <- paste0("Bartlett test [bartlett.test()] - Heterogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances hétérogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="on")
        }

      } else if (check_normality==FALSE) {
        #========================================
        # DONNÉES NON NORMALES : Levene (puis éventuellement Brown-Forsyth)
        #========================================
        if (code==TRUE){
          cat("# Tests de Levene et Brown-Forsyth\n")
          cat("library(car)\n")
          cat("leveneTest(x, g)\n")
          cat("library(onewaytests)\n")
          cat("bf.test(x~g, data=data.frame(x, 'g'=factor(g)))\n")
        }

        # Test de Levene (basé sur la moyenne)
        pvals <- suppressWarnings(car::leveneTest(x, g))[1, 3]

        k <- .vbse(
          "Academic check of the homogeneity of group variances.",
          "Contrôle académique de l'homogénéité de la variance des groupes.",
          verbose = verbose, k = k, cpt="on"
        )

        if (pvals <= alpha) {
          # VARIANCE HÉTÉROGÈNE détectée par Levene
          ang <- paste0("Levene test detects heteroscedasticity (heterogeneous variances).\n\t",
                        "p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Le test de Levene met en évidence une hétéroscédasticité (variances hétérogènes).\n\t",
                       "p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir été faussé par une asymétrie des données ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="off")

          #========================================
          # Test de Brown-Forsyth (plus robuste, basé sur la médiane)
          # UNIQUEMENT si Levene a détecté une hétérogénéité
          #========================================
          pvals2 <- suppressWarnings(bf.test(x~g, data=data.frame(x, "g"=factor(g)), verbose=FALSE))$p.value

          if (pvals2 <= alpha) {
            # VARIANCE HÉTÉROGÈNE confirmée par Brown-Forsyth
            check_variance <- FALSE
            ang <- paste0("Brown-Forsyth test confirms heteroscedasticity.\n\t",
                          "p-value: ", .format_pval(pvals2), ".")
            fr <- paste0("Le test de Brown-Forsyth confirme l'hétéroscédasticité.\n\t",
                         "p-value : ", .format_pval(pvals2), ".")
            k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="off")
          } else {
            # Brown-Forsyth ne confirme PAS → Levene probablement faussé
            check_variance <- TRUE
            ang <- paste0("Brown-Forsyth test does NOT confirm heteroscedasticity.\n\t",
                          "p-value: ", .format_pval(pvals2), ".\n\t",
                          "→ Variances considered homogeneous (Levene likely affected by asymmetry/outliers).")
            fr <- paste0("Le test de Brown-Forsyth ne confirme PAS l'hétéroscédasticité.\n\t",
                         "p-value : ", .format_pval(pvals2), ".\n\t",
                         "→ Variances considérées comme homogènes (Levene probablement faussé par asymétrie/outliers).")
            k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="off")
          }

        } else {
          # VARIANCE HOMOGÈNE selon Levene → PAS de test Brown-Forsyth
          check_variance <- TRUE
          ang <- paste0("Levene test shows homogeneous variances (homoscedasticity).\n\t",
                        "p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Le test de Levene montre que la variance est homogène (homoscédasticité).\n\t",
                       "p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir été faussé par une asymétrie des données ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, k = k, cpt="off")
        }

      } # Fin check_normality==FALSE

    } # Fin n_g > 2

  } # Fin paired==FALSE

  # Retour de la liste
  variance <- list()
  variance[[1]] <- check_variance
  variance[[2]] <- k
  return(variance)

} # Fin de la fonction

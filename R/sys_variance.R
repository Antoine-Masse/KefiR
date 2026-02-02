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
#' @param assumption_label Character. Label for assumption numbering (e.g., "3/3", "2/3", NULL for no label)
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
                      paired=FALSE, debug=FALSE, verbose=FALSE, k=0, code=FALSE, cpt="on",
                      assumption_label=NULL) {
  ####################################
  # Initialisation
  ####################################
  n_g <- length(unique(g))
  check_variance <- TRUE  # Par défaut

  if (paired==TRUE) {
    #========================================
    # CAS APPARIÉ : Pas de test nécessaire (plan within pur)
    #========================================
    # Note technique: Pour un plan MIXTE (within + between), l'homogénéité des variances
    # devrait être testée pour les facteurs between. Cependant, cette fonction est appelée
    # dans un contexte où le test de variance n'est pas nécessaire pour les comparaisons
    # intra-sujet (within). Si un facteur between est présent, le test devrait être géré
    # séparément dans .multi_factor_analysis.R avant appel à cette fonction.

    # NE PAS AFFICHER ce message pour mesures répétées pures (cpt="off" si pas verbose)
    # Le message ne sera affiché QUE si verbose=TRUE ET que l'utilisateur veut des détails
    # Pour éviter confusion, on ne l'affiche PAS du tout (cpt="off" désactive numérotation)
    # Message supprimé complètement pour éviter confusion utilisateur
    # k <- .vbse(..., verbose = FALSE, k = k, cpt="off")  # Désactivé

    check_variance <- TRUE

  } else if (paired==FALSE) {

    if (n_g == 2) {
      #========================================
      # 2 GROUPES : Test de Fisher-Snedecor
      #========================================
      if (isTRUE(code)){
        cat("# Test de Fisher-Snedecor\n")
        cat("var.test(x~g)\n")
      }

      fm <- formula(x~g)
      pvals <- var.test(fm)$p.value

      if (pvals > alpha) {
        # VARIANCES HOMOGÈNES
        check_variance <- TRUE
        if (check_normality == FALSE) {
          # Cas non-normal + variances homogènes : passage vers contrôle tolérance
          # pour tester éventualité retour vers paramétrique
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Homogeneous variances (p = ", .format_pval(pvals), ").\n\t",
                       "--> Towards a less stringent normality check for...\n\t",
                       "    ...testing the possibility of returning to a parametric comparison test.")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances homogènes (p = ", .format_pval(pvals), ").\n\t",
                      "--> Vers un contrôle de la normalité moins exigeant pour...\n\t",
                      "    ...tester éventualité d'un retour à un test de comparaison paramétrique.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        } else {
          # Cas normal : pas de message directionnel
          ang <- paste0("Fisher-Snedecor test [var.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Fisher-Snedecor [var.test()] - Variances homogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        }
      } else {
        # VARIANCES HÉTÉROGÈNES
        check_variance <- FALSE
        if (check_normality == FALSE) {
          # Cas non-normal : indiquer passage vers test non-paramétrique
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Heterogeneous variances (p = ", .format_pval(pvals), ").\n\t",
                       "--> Towards a non-parametric comparison test (Wilcoxon-Mann-Whitney)")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances hétérogènes (p = ", .format_pval(pvals), ").\n\t",
                      "--> Vers un test de comparaison non paramétrique (Wilcoxon-Mann-Whitney)")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        } else {
          # Cas normal : pas de message directionnel
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Heterogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances hétérogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        }
      }

    } else if (n_g > 2) {

      if (check_normality==TRUE) {
        #========================================
        # DONNÉES NORMALES : Test de Bartlett
        #========================================

        # Vérification précoce : au moins 2 observations par groupe
        group_counts <- table(g)
        groups_with_lt2 <- names(group_counts)[group_counts < 2]

        if (length(groups_with_lt2) > 0) {
          k <- .vbse(
            paste0("WARNING: Cannot perform Bartlett test - insufficient observations.\n\t",
                   "Groups with < 2 observations: ", paste(groups_with_lt2, collapse = ", "), "\n\t",
                   "Bartlett test requires at least 2 observations per group.\n\t",
                   "==> Variance homogeneity assumed (test skipped)."),
            paste0("ATTENTION : Impossible d'effectuer le test de Bartlett - observations insuffisantes.\n\t",
                   "Groupes avec < 2 observations : ", paste(groups_with_lt2, collapse = ", "), "\n\t",
                   "Le test de Bartlett nécessite au moins 2 observations par groupe.\n\t",
                   "==> Homogénéité des variances assumée (test ignoré)."),
            verbose = verbose, code = code, k = k, cpt = cpt
          )
          check_variance <- TRUE
          return(list(check_variance, k))
        }

        if (isTRUE(code)){
          cat("# Test de Bartlett\n")
          cat("bartlett.test(x, g)\n")
        }

        pvals <- bartlett.test(x, g)$p.value

        # Message de contrôle académique avec numérotation ASSOMPTION (si fournie)
        if (!is.null(assumption_label)) {
          k <- .vbse(
            paste0("ASSUMPTION ", assumption_label, ": Academic check of group variance homogeneity."),
            paste0("ASSOMPTION ", assumption_label, " : Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes."),
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        } else {
          k <- .vbse(
            "Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes.",
            "Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes.",
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        }

        if (pvals > alpha) {
          # VARIANCES HOMOGÈNES
          check_variance <- TRUE
          ang <- paste0("Bartlett test [bartlett.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances homogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
        } else {
          # VARIANCES HÉTÉROGÈNES
          check_variance <- FALSE
          ang <- paste0("Bartlett test [bartlett.test()] - Heterogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances hétérogènes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
        }

      } else if (check_normality==FALSE) {
        #========================================
        # DONNÉES NON NORMALES : Levene (puis éventuellement Brown-Forsyth)
        #========================================

        # Vérification précoce : au moins 2 observations par groupe
        group_counts <- table(g)
        groups_with_lt2 <- names(group_counts)[group_counts < 2]

        if (length(groups_with_lt2) > 0) {
          k <- .vbse(
            paste0("WARNING: Cannot perform Levene test - insufficient observations.\n\t",
                   "Groups with < 2 observations: ", paste(groups_with_lt2, collapse = ", "), "\n\t",
                   "Levene test requires at least 2 observations per group.\n\t",
                   "==> Variance homogeneity assumed (test skipped)."),
            paste0("ATTENTION : Impossible d'effectuer le test de Levene - observations insuffisantes.\n\t",
                   "Groupes avec < 2 observations : ", paste(groups_with_lt2, collapse = ", "), "\n\t",
                   "Le test de Levene nécessite au moins 2 observations par groupe.\n\t",
                   "==> Homogénéité des variances assumée (test ignoré)."),
            verbose = verbose, code = code, k = k, cpt = cpt
          )
          check_variance <- TRUE
          return(list(check_variance, k))
        }

        if (isTRUE(code)){
          cat("# Tests de Levene et Brown-Forsyth\n")
          cat("library(car)\n")
          cat("leveneTest(x, g)\n")
          cat("library(onewaytests)\n")
          cat("bf.test(x~g, data=data.frame(x, 'g'=factor(g)))\n")
        }

        # Test de Levene (basé sur la moyenne)
        pvals <- suppressWarnings(car::leveneTest(x, g))[1, 3]

        # Message de contrôle académique avec numérotation ASSOMPTION (si fournie)
        if (!is.null(assumption_label)) {
          k <- .vbse(
            paste0("ASSUMPTION ", assumption_label, ": Academic check of group variance homogeneity."),
            paste0("ASSOMPTION ", assumption_label, " : Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes."),
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        } else {
          k <- .vbse(
            "Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes.",
            "Contrôle ACADÉMIQUE de l'homogénéité de la variance des groupes.",
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        }

        if (pvals <= alpha) {
          # VARIANCE HÉTÉROGÈNE détectée par Levene
          ang <- paste0("Levene test detects heteroscedasticity (heterogeneous variances).\n\t",
                        "\t==> p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Le test de Levene met en évidence une hétéroscédasticité (variances hétérogènes).\n\t",
                       "\t==> p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir été faussé par une asymétrie des données ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

          #========================================
          # Test de Brown-Forsyth (plus robuste, basé sur la médiane)
          # UNIQUEMENT si Levene a détecté une hétérogénéité
          #========================================
          pvals2 <- suppressWarnings(bf.test(x~g, data=data.frame(x, "g"=factor(g)), verbose=FALSE))$p.value

          if (pvals2 <= alpha) {
            # VARIANCE HÉTÉROGÈNE confirmée par Brown-Forsyth
            check_variance <- FALSE
            ang <- paste0("Brown-Forsyth test confirms heteroscedasticity.\n\t",
                          "\t==> p-value: ", .format_pval(pvals2), ".")
            fr <- paste0("Le test de Brown-Forsyth confirme l'hétéroscédasticité.\n\t",
                         "\t==> p-value : ", .format_pval(pvals2), ".")
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

            # Message de redirection vers non-paramétrique
            ang <- "--> Moving towards non-parametric analysis."
            fr <- "--> On part vers une analyse non-paramétrique."
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
          } else {
            # Brown-Forsyth ne confirme PAS → Levene probablement faussé
            check_variance <- TRUE
            ang <- paste0("Brown-Forsyth test does NOT confirm heteroscedasticity.\n\t",
                          "\t==> p-value: ", .format_pval(pvals2), ".\n\t",
                          "==> Variances considered homogeneous (Levene likely affected by asymmetry/outliers).")
            fr <- paste0("Le test de Brown-Forsyth ne confirme PAS l'hétéroscédasticité.\n\t",
                         "\t==> p-value : ", .format_pval(pvals2), ".\n\t",
                         "==> Variances considérées comme homogènes (Levene probablement faussé par asymétrie/outliers).")
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
          }

        } else {
          # VARIANCE HOMOGÈNE selon Levene → PAS de test Brown-Forsyth
          check_variance <- TRUE
          ang <- paste0("Levene test [car::leveneTest()] shows homogeneous variances (homoscedasticity).\n\t",
                        "p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Test de Levene [car::leveneTest()] - Variance homogène (homoscédasticité).\n\t",
                       "p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir été faussé par une asymétrie des données ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

          # Message de transition vers contrôle tolérant de normalité
          ang <- "--> Variance homogeneity invites a more tolerant normality check."
          fr <- "--> L'homogénéité de la variance invite à un contrôle plus tolérant de la normalité."
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
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

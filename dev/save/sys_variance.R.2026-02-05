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
  check_variance <- TRUE  # Par d\u00e9faut

  if (paired==TRUE) {
    #========================================
    # CAS APPARI\u00c9 : Pas de test n\u00e9cessaire (plan within pur)
    #========================================
    # Note technique: Pour un plan MIXTE (within + between), l'homog\u00e9n\u00e9it\u00e9 des variances
    # devrait \u00eatre test\u00e9e pour les facteurs between. Cependant, cette fonction est appel\u00e9e
    # dans un contexte o\u00f9 le test de variance n'est pas n\u00e9cessaire pour les comparaisons
    # intra-sujet (within). Si un facteur between est pr\u00e9sent, le test devrait \u00eatre g\u00e9r\u00e9
    # s\u00e9par\u00e9ment dans .multi_factor_analysis.R avant appel \u00e0 cette fonction.

    # NE PAS AFFICHER ce message pour mesures r\u00e9p\u00e9t\u00e9es pures (cpt="off" si pas verbose)
    # Le message ne sera affich\u00e9 QUE si verbose=TRUE ET que l'utilisateur veut des d\u00e9tails
    # Pour \u00e9viter confusion, on ne l'affiche PAS du tout (cpt="off" d\u00e9sactive num\u00e9rotation)
    # Message supprim\u00e9 compl\u00e8tement pour \u00e9viter confusion utilisateur
    # k <- .vbse(..., verbose = FALSE, k = k, cpt="off")  # D\u00e9sactiv\u00e9

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
        # VARIANCES HOMOG\u00c8NES
        check_variance <- TRUE
        if (check_normality == FALSE) {
          # Cas non-normal + variances homog\u00e8nes : passage vers contr\u00f4le tol\u00e9rance
          # pour tester \u00e9ventualit\u00e9 retour vers param\u00e9trique
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Homogeneous variances (p = ", .format_pval(pvals), ").\n\t",
                       "--> Towards a less stringent normality check for...\n\t",
                       "    ...testing the possibility of returning to a parametric comparison test.")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances homog\u00e8nes (p = ", .format_pval(pvals), ").\n\t",
                      "--> Vers un contr\u00f4le de la normalit\u00e9 moins exigeant pour...\n\t",
                      "    ...tester \u00e9ventualit\u00e9 d'un retour \u00e0 un test de comparaison param\u00e9trique.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        } else {
          # Cas normal : pas de message directionnel
          ang <- paste0("Fisher-Snedecor test [var.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Fisher-Snedecor [var.test()] - Variances homog\u00e8nes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        }
      } else {
        # VARIANCES H\u00c9T\u00c9ROG\u00c8NES
        check_variance <- FALSE
        if (check_normality == FALSE) {
          # Cas non-normal : indiquer passage vers test non-param\u00e9trique
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Heterogeneous variances (p = ", .format_pval(pvals), ").\n\t",
                       "--> Towards a non-parametric comparison test (Wilcoxon-Mann-Whitney)")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances h\u00e9t\u00e9rog\u00e8nes (p = ", .format_pval(pvals), ").\n\t",
                      "--> Vers un test de comparaison non param\u00e9trique (Wilcoxon-Mann-Whitney)")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        } else {
          # Cas normal : pas de message directionnel
          ang <- paste0("Fisher-Snedecor test [var.test()]\n\t",
                       "==> Heterogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Fisher-Snedecor [var.test()]\n\t",
                      "==> Variances h\u00e9t\u00e9rog\u00e8nes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="on")
        }
      }

    } else if (n_g > 2) {

      if (check_normality==TRUE) {
        #========================================
        # DONN\u00c9ES NORMALES : Test de Bartlett
        #========================================

        # V\u00e9rification pr\u00e9coce : au moins 2 observations par groupe
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
                   "Le test de Bartlett n\u00e9cessite au moins 2 observations par groupe.\n\t",
                   "==> Homog\u00e9n\u00e9it\u00e9 des variances assum\u00e9e (test ignor\u00e9)."),
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

        # Message de contr\u00f4le acad\u00e9mique avec num\u00e9rotation ASSOMPTION (si fournie)
        if (!is.null(assumption_label)) {
          k <- .vbse(
            paste0("ASSUMPTION ", assumption_label, ": Academic check of group variance homogeneity."),
            paste0("ASSOMPTION ", assumption_label, " : Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes."),
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        } else {
          k <- .vbse(
            "Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes.",
            "Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes.",
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        }

        if (pvals > alpha) {
          # VARIANCES HOMOG\u00c8NES
          check_variance <- TRUE
          ang <- paste0("Bartlett test [bartlett.test()] - Homogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances homog\u00e8nes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
        } else {
          # VARIANCES H\u00c9T\u00c9ROG\u00c8NES
          check_variance <- FALSE
          ang <- paste0("Bartlett test [bartlett.test()] - Heterogeneous variances (p = ", .format_pval(pvals), ").")
          fr <- paste0("Test de Bartlett [bartlett.test()] - Variances h\u00e9t\u00e9rog\u00e8nes (p = ", .format_pval(pvals), ").")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
        }

      } else if (check_normality==FALSE) {
        #========================================
        # DONN\u00c9ES NON NORMALES : Levene (puis \u00e9ventuellement Brown-Forsyth)
        #========================================

        # V\u00e9rification pr\u00e9coce : au moins 2 observations par groupe
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
                   "Le test de Levene n\u00e9cessite au moins 2 observations par groupe.\n\t",
                   "==> Homog\u00e9n\u00e9it\u00e9 des variances assum\u00e9e (test ignor\u00e9)."),
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

        # Test de Levene (bas\u00e9 sur la moyenne)
        pvals <- suppressWarnings(car::leveneTest(x, g))[1, 3]

        # Message de contr\u00f4le acad\u00e9mique avec num\u00e9rotation ASSOMPTION (si fournie)
        if (!is.null(assumption_label)) {
          k <- .vbse(
            paste0("ASSUMPTION ", assumption_label, ": Academic check of group variance homogeneity."),
            paste0("ASSOMPTION ", assumption_label, " : Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes."),
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        } else {
          k <- .vbse(
            "Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes.",
            "Contr\u00f4le ACAD\u00c9MIQUE de l'homog\u00e9n\u00e9it\u00e9 de la variance des groupes.",
            verbose = verbose, code = code, k = k, cpt=cpt
          )
        }

        if (pvals <= alpha) {
          # VARIANCE H\u00c9T\u00c9ROG\u00c8NE d\u00e9tect\u00e9e par Levene
          ang <- paste0("Levene test detects heteroscedasticity (heterogeneous variances).\n\t",
                        "\t==> p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Le test de Levene met en \u00e9vidence une h\u00e9t\u00e9rosc\u00e9dasticit\u00e9 (variances h\u00e9t\u00e9rog\u00e8nes).\n\t",
                       "\t==> p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir \u00e9t\u00e9 fauss\u00e9 par une asym\u00e9trie des donn\u00e9es ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

          #========================================
          # Test de Brown-Forsyth (plus robuste, bas\u00e9 sur la m\u00e9diane)
          # UNIQUEMENT si Levene a d\u00e9tect\u00e9 une h\u00e9t\u00e9rog\u00e9n\u00e9it\u00e9
          #========================================
          pvals2 <- suppressWarnings(bf.test(x~g, data=data.frame(x, "g"=factor(g)), verbose=FALSE))$p.value

          if (pvals2 <= alpha) {
            # VARIANCE H\u00c9T\u00c9ROG\u00c8NE confirm\u00e9e par Brown-Forsyth
            check_variance <- FALSE
            ang <- paste0("Brown-Forsyth test confirms heteroscedasticity.\n\t",
                          "\t==> p-value: ", .format_pval(pvals2), ".")
            fr <- paste0("Le test de Brown-Forsyth confirme l'h\u00e9t\u00e9rosc\u00e9dasticit\u00e9.\n\t",
                         "\t==> p-value : ", .format_pval(pvals2), ".")
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

            # Message de redirection vers non-param\u00e9trique
            ang <- "--> Moving towards non-parametric analysis."
            fr <- "--> On part vers une analyse non-param\u00e9trique."
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
          } else {
            # Brown-Forsyth ne confirme PAS \u2192 Levene probablement fauss\u00e9
            check_variance <- TRUE
            ang <- paste0("Brown-Forsyth test does NOT confirm heteroscedasticity.\n\t",
                          "\t==> p-value: ", .format_pval(pvals2), ".\n\t",
                          "==> Variances considered homogeneous (Levene likely affected by asymmetry/outliers).")
            fr <- paste0("Le test de Brown-Forsyth ne confirme PAS l'h\u00e9t\u00e9rosc\u00e9dasticit\u00e9.\n\t",
                         "\t==> p-value : ", .format_pval(pvals2), ".\n\t",
                         "==> Variances consid\u00e9r\u00e9es comme homog\u00e8nes (Levene probablement fauss\u00e9 par asym\u00e9trie/outliers).")
            k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")
          }

        } else {
          # VARIANCE HOMOG\u00c8NE selon Levene \u2192 PAS de test Brown-Forsyth
          check_variance <- TRUE
          ang <- paste0("Levene test [car::leveneTest()] shows homogeneous variances (homoscedasticity).\n\t",
                        "p-value: ", .format_pval(pvals), ".\n\t",
                        "Note: Levene test may be affected by data asymmetry or outliers.")
          fr <- paste0("Test de Levene [car::leveneTest()] - Variance homog\u00e8ne (homosc\u00e9dasticit\u00e9).\n\t",
                       "p-value : ", .format_pval(pvals), ".\n\t",
                       "Note : Le test de Levene peut avoir \u00e9t\u00e9 fauss\u00e9 par une asym\u00e9trie des donn\u00e9es ou des outliers.")
          k <- .vbse(ang, fr, verbose = verbose, code = code, k = k, cpt="off")

          # Message de transition vers contr\u00f4le tol\u00e9rant de normalit\u00e9
          ang <- "--> Variance homogeneity invites a more tolerant normality check."
          fr <- "--> L'homog\u00e9n\u00e9it\u00e9 de la variance invite \u00e0 un contr\u00f4le plus tol\u00e9rant de la normalit\u00e9."
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

#' Post-hoc analyses for MANOVA
#'
#' Performs appropriate follow-up analyses after a significant MANOVA result.
#' Includes discriminant analysis to identify which dependent variables
#' separate the groups, and protected univariate ANOVAs with Bonferroni correction.
#'
#' @param x Matrix of dependent variables (p columns = p DVs)
#' @param g Factor vector of group memberships
#' @param manova_result The MANOVA result object from .manova_analysis()
#' @param alpha Significance level (default = 0.05)
#' @param method Method for post-hoc: "discriminant" (default), "protected_anova", or "both"
#' @param verbose Logical, print detailed messages
#' @param code Logical. If TRUE, suppresses verbose messages for code-only mode (default: FALSE)
#' @param debug Logical, print debug messages
#' @param k Message counter
#'
#' @return List containing:
#'   \itemize{
#'     \item discriminant_analysis: Results from discriminant analysis (if requested)
#'     \item protected_anovas: Results from protected univariate ANOVAs (if requested)
#'     \item interpretation: Text interpretation of results
#'   }
#'
#' @details
#' **APPROPRIATE POST-HOC METHODS FOR MANOVA**
#'
#' After a significant MANOVA, standard univariate post-hocs (Tukey, etc.) are
#' inappropriate because they ignore the multivariate structure and correlations
#' between dependent variables.
#'
#' **Recommended approaches** (Huberty & Olejnik, 2006):
#'
#' 1. **Discriminant Analysis** (Primary recommendation)
#'    - Identifies which linear combination of DVs best separates groups
#'    - Shows which DVs contribute most to group separation
#'    - Provides canonical discriminant functions
#'    - Reference: Huberty & Olejnik (2006). Applied MANOVA and Discriminant Analysis, 2nd ed.
#'
#' 2. **Protected Univariate ANOVAs**
#'    - Perform separate ANOVA on each DV
#'    - Apply Bonferroni correction: alpha_adjusted = alpha / p (p = number of DVs)
#'    - Only valid if MANOVA was significant (protection against Type I error inflation)
#'    - Reference: Stevens, J. P. (2009). Applied Multivariate Statistics for the Social Sciences, 5th ed.
#'
#' 3. **Roy-Bargmann Stepdown Analysis** (Advanced, not implemented yet)
#'    - Tests DVs sequentially in order of theoretical importance
#'    - Controls for previously tested DVs
#'    - Reference: Roy, S. N., & Bargmann, R. E. (1958). Tests of multiple independence
#'
#' **NOT recommended:**
#' - Pairwise multivariate comparisons (complex, low power)
#' - Univariate post-hocs without MANOVA protection (inflates Type I error)
#'
#' @note
#' This function should only be called after a significant MANOVA result.
#' The MANOVA provides protection against Type I error inflation when
#' examining individual dependent variables.
#'
#' @references
#' Huberty, C. J., & Olejnik, S. (2006). Applied MANOVA and Discriminant Analysis (2nd ed.). Wiley.
#' Stevens, J. P. (2009). Applied Multivariate Statistics for the Social Sciences (5th ed.). Routledge.
#' Tabachnick, B. G., & Fidell, L. S. (2019). Using Multivariate Statistics (7th ed.). Pearson.
#'
#' @examples
#' \dontrun{
#' # After significant MANOVA
#' manova_res <- .manova_analysis(cbind(A,B,C) ~ Group, data=data)
#' if (manova_res$test_statistics["Wilks", "p.value"] < 0.05) {
#'   posthoc <- .posthoc_MANOVA(
#'     x = manova_res$x,
#'     g = manova_res$g,
#'     manova_result = manova_res,
#'     method = "both",
#'     verbose = TRUE
#'   )
#' }
#' }
#' @keywords internal
.posthoc_MANOVA <- function(x = NULL,
                            g = NULL,
                            manova_result = NULL,
                            alpha = 0.05,
                            method = "both",
                            verbose = FALSE,
                            code = FALSE,
                            debug = FALSE,
                            k = NULL) {

  # =============================================================================
  # VALIDATION
  # =============================================================================

  if (is.null(x) || is.null(g)) {
    .exit("x and g are required for MANOVA post-hoc analysis",
          "x et g sont requis pour les post-hocs MANOVA")
  }

  if (!is.matrix(x)) {
    .exit("x must be a matrix for MANOVA post-hoc (multiple dependent variables)",
          "x doit \u00eatre une matrice pour post-hocs MANOVA (plusieurs variables d\u00e9pendantes)")
  }

  if (!is.factor(g)) {
    g <- as.factor(g)
  }

  # Vérifier que MANOVA était significative
  if (!is.null(manova_result) && !is.null(manova_result$test_statistics)) {
    # Test retenu par la MANOVA ; avec plusieurs facteurs, au moins un terme significatif suffit
    effects <- manova_result$effects
    manova_p <- if (!is.null(effects)) {
      min(if (isFALSE(manova_result$check_covariance_homogeneity) ||
              isFALSE(manova_result$check_multivariate_normality)) effects$p.Pillai else effects$p.Wilks)
    } else {
      manova_result$global_pvalue
    }

    if (isTRUE(manova_p >= alpha)) {
      k <- .vbse(
        paste0("WARNING: MANOVA was not significant (p = ", .format_pval(manova_p),
               ").\n\tThe post-hocs below are descriptive."),
        paste0("ATTENTION : la MANOVA n'est pas significative (p = ", .format_pval(manova_p),
               ").\n\tLes post-hocs ci-dessous sont descriptifs."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }
  }

  # =============================================================================
  # PARAMÈTRES
  # =============================================================================

  n <- nrow(x)
  p <- ncol(x)  # Nombre de variables dépendantes
  n_groups <- length(levels(g))

  k <- .vbse(
    paste0("POST-HOC MANOVA ANALYSIS\n",
           "\t", n, " observations, ", p, " DVs, ", n_groups, " groups\n",
           "\tMethod: ", method),
    paste0("ANALYSE POST-HOC MANOVA\n",
           "\t", n, " observations, ", p, " VDs, ", n_groups, " groupes\n",
           "\tM\u00e9thode : ", method),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  result <- list()

  # =============================================================================
  # 1. ANALYSE DISCRIMINANTE
  # =============================================================================

  if (method %in% c("discriminant", "both")) {

    .dbg("Performing discriminant analysis...",
         "Analyse discriminante en cours...", debug = debug)

    # Vérifier si package MASS disponible
    if (requireNamespace("MASS", quietly = TRUE)) {

      # Créer data frame pour lda()
      df <- as.data.frame(x)
      colnames(df) <- paste0("DV", 1:p)
      df$Group <- g

      # Analyse discriminante linéaire
      lda_result <- tryCatch({
        MASS::lda(Group ~ ., data = df)
      }, error = function(e) {
        k <<- .vbse(
          paste0("ERROR in discriminant analysis: ", e$message),
          paste0("ERREUR dans l'analyse discriminante: ", e$message),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        return(NULL)
      })

      if (!is.null(lda_result)) {

        # Proportions de variance expliquée par chaque fonction discriminante
        prop_trace <- lda_result$svd^2 / sum(lda_result$svd^2)

        # Coefficients standardisés (structure matrix)
        # Corrélations entre DVs et fonctions discriminantes
        lda_predict <- predict(lda_result)
        # Matrice de structure INTRA-groupes, la d\u00e9finition usuelle (corr\u00e9lations sur les
        # donn\u00e9es brutes : gonfl\u00e9es par l'\u00e9cart entre groupes)
        # R\u00e9f\u00e9rence: Huberty & Olejnik (2006). Applied MANOVA and Discriminant Analysis. Wiley.
        correlations <- cor(residuals(lm(x ~ g)), residuals(lm(lda_predict$x ~ g)))

        # Message d'en-tête compact
        n_functions <- length(lda_result$svd)
        k <- .vbse(
          paste0("DISCRIMINANT ANALYSIS:\n",
                 "\tNumber of discriminant functions: ", n_functions),
          paste0("ANALYSE DISCRIMINANTE :\n",
                 "\tNombre de fonctions discriminantes : ", n_functions),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # Afficher proportion de variance pour chaque fonction (max 3)
        if (n_functions > 0) {
          prop_messages_en <- character()
          prop_messages_fr <- character()
          for (i in 1:min(n_functions, 3)) {
            prop_messages_en <- c(prop_messages_en,
                                  paste0("\tLD", i, ": ", round(prop_trace[i] * 100), "%"))
            prop_messages_fr <- c(prop_messages_fr,
                                  paste0("\tFD", i, " : ", round(prop_trace[i] * 100), "%"))
          }

          k <- .vbse(
            paste0("Proportion of between-group variance explained:\n",
                   paste(prop_messages_en, collapse = "\n")),
            paste0("Proportion de variance inter-groupes expliqu\u00e9e :\n",
                   paste(prop_messages_fr, collapse = "\n")),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }

        # Identifier les DVs les plus importantes pour chaque fonction (max 2)
        for (i in 1:min(ncol(correlations), 2)) {
          max_idx <- which.max(abs(correlations[, i]))
          max_cor <- correlations[max_idx, i]
          dv_name <- if (!is.null(colnames(x))) colnames(x)[max_idx] else paste0("DV", max_idx)

          k <- .vbse(
            paste0("\tLD", i, ": Main variable = ", dv_name, " (r = ", round(max_cor, 3), ")"),
            paste0("\tFD", i, " : Variable la plus importante = ", dv_name, " (r = ", round(max_cor, 3), ")"),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }

        # Stocker SEULEMENT les résumés (pas les objets complets)
        # Pour éviter un retour trop verbeux
        result$discriminant_analysis <- list(
          n_functions = n_functions,
          prop_variance = prop_trace,
          structure_matrix = correlations
          # NOTE: lda et predictions disponibles en relançant MASS::lda() si besoin
        )
      }

    } else {
      k <- .vbse(
        "Package 'MASS' not available. Discriminant analysis skipped.\n\tInstall with: install.packages('MASS')",
        "Package 'MASS' non disponible. Analyse discriminante ignor\u00e9e.\n\tInstaller avec: install.packages('MASS')",
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }
  }

  # =============================================================================
  # 2. ANOVAS UNIVARIÉES PROTÉGÉES (Bonferroni)
  # =============================================================================

  if (method %in% c("protected_anova", "both")) {

    .dbg("Performing protected univariate ANOVAs...",
         "ANOVAs univari\u00e9es prot\u00e9g\u00e9es en cours...", debug = debug)

    # Ajustement Bonferroni
    alpha_adjusted <- alpha / p

    k <- .vbse(
      paste0("PROTECTED UNIVARIATE ANOVAs (Bonferroni correction):\n",
             "\tOriginal alpha: ", alpha, "\n",
             "\tAdjusted alpha: ", round(alpha_adjusted, 3), " (alpha / ", p, " DVs)"),
      paste0("ANOVAs UNIVARI\u00c9ES PROT\u00c9G\u00c9ES (correction Bonferroni) :\n",
             "\tAlpha original : ", alpha, "\n",
             "\tAlpha ajust\u00e9 : ", round(alpha_adjusted, 3), " (alpha / ", p, " VDs)"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # Test de chaque VD coh\u00e9rent avec les assomptions contr\u00f4l\u00e9es par la MANOVA
    dv_route <- if (isTRUE(manova_result$robust)) {
      "kruskal"
    } else if (isFALSE(manova_result$check_covariance_homogeneity)) {
      "welch"
    } else {
      "anova"
    }
    route_label <- switch(dv_route, kruskal = "Kruskal-Wallis [kruskal.test()]",
                          welch = "Welch [oneway.test()]", anova = "ANOVA [aov()]")
    k <- .vbse(paste0("\tTest used for each variable: ", route_label),
               paste0("\tTest utilis\u00e9 pour chaque variable : ", route_label),
               verbose = verbose, code = code, k = k, cpt = "off")

    anova_results <- list()
    sig_dvs <- character()
    nonsig_dvs <- character()

    for (i in 1:p) {
      dv_name <- if (!is.null(colnames(x))) colnames(x)[i] else paste0("DV", i)

      # ANOVA pour cette DV
      dv_data <- x[, i]
      if (dv_route == "kruskal") {
        kw <- kruskal.test(dv_data ~ g)
        f_stat <- unname(kw$statistic); df1 <- unname(kw$parameter); df2 <- NA; p_value <- kw$p.value
        stat_label <- paste0("H(", df1, ")")
      } else if (dv_route == "welch") {
        welch <- oneway.test(dv_data ~ g, var.equal = FALSE)
        f_stat <- unname(welch$statistic); df1 <- unname(welch$parameter[1]); df2 <- unname(welch$parameter[2])
        p_value <- welch$p.value
        stat_label <- paste0("F(", df1, ",", round(df2, 1), ")")
      } else {
        aov_summary <- summary(aov(dv_data ~ g))
        f_stat <- aov_summary[[1]]["g", "F value"]
        p_value <- aov_summary[[1]]["g", "Pr(>F)"]
        df1 <- aov_summary[[1]]["g", "Df"]
        df2 <- aov_summary[[1]]["Residuals", "Df"]
        stat_label <- paste0("F(", df1, ",", df2, ")")
      }

      # Test de significativité avec alpha ajusté
      is_sig <- p_value < alpha_adjusted

      # Collecter pour synthèse
      dv_result <- paste0(dv_name, ": ", stat_label, "=", round(f_stat, 2),
                          ", p=", .format_pval(p_value))

      # Comparaison des groupes pour cette VD, au seuil corrig\u00e9 de Bonferroni
      dv_pairs <- if (dv_route == "kruskal") {
        suppressWarnings(pairwise.wilcox.test(dv_data, g, p.adjust.method = "holm", exact = FALSE))
      } else {
        pairwise.t.test(dv_data, g, p.adjust.method = "holm", pool.sd = (dv_route == "anova"))
      }
      dv_letters <- .letters_from_pairs(dv_pairs$p.value, tapply(dv_data, g, mean), alpha_adjusted)
      colnames(dv_letters)[2] <- if (dv_route == "kruskal") "Wilcoxon_Holm" else "Student_Holm"

      if (is_sig) {
        sig_dvs <- c(sig_dvs, dv_result)
      } else {
        nonsig_dvs <- c(nonsig_dvs, dv_result)
      }

      # Stocker SEULEMENT les statistiques (pas le modèle complet)
      # Pour éviter un retour trop verbeux
      anova_results[[dv_name]] <- list(
        F = f_stat,
        df1 = df1,
        df2 = df2,
        p_value = p_value,
        significant = is_sig,
        alpha_adjusted = alpha_adjusted,
        test = dv_route,
        groups = dv_letters
        # NOTE: model aov disponible en relançant aov() si besoin
      )
    }

    result$protected_anovas <- anova_results

    # Plusieurs facteurs : le test de chaque terme pour chaque VD (Type II), compl\u00e9mentaire
    # de la comparaison des groupes crois\u00e9s ci-dessus
    if (!is.null(manova_result$term_labels) && !is.null(manova_result$g_terms)) {
      term_formula <- paste("dv ~", paste(manova_result$term_labels, collapse = " + "))
      result$effects_by_dv <- lapply(stats::setNames(seq_len(p), colnames(x)), function(i_dv) {
        dv_data <- data.frame(dv = x[, i_dv], manova_result$g_terms)
        .effects_from_table(suppressMessages(car::Anova(lm(as.formula(term_formula), data = dv_data), type = "II")))
      })
      k <- .vbse("Test of each term for each variable, Type II [car::Anova()]:",
                 "Test de chaque terme pour chaque variable, Type II [car::Anova()] :",
                 verbose = verbose, code = code, k = k, cpt = "off")
      if (verbose) {
        for (dv_name in names(result$effects_by_dv)) {
          cat("\t", dv_name, "\n", sep = "")
          print(result$effects_by_dv[[dv_name]], row.names = FALSE, digits = 4)
        }
        cat("\n")
      }
    }

    # Synthèse des résultats
    n_sig <- length(sig_dvs)

    # Afficher chaque VD : une absence d'effet est aussi un r\u00e9sultat
    for (dv_info in sig_dvs) {
      k <- .vbse(
        paste0("\t", dv_info, " *SIGNIFICANT*"),
        paste0("\t", dv_info, " *SIGNIFICATIF*"),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
    for (dv_info in nonsig_dvs) {
      k <- .vbse(paste0("\t", dv_info), paste0("\t", dv_info),
                 verbose = verbose, code = code, k = k, cpt = "off")
    }

    # Message de résumé final (cpt = "off" car c'est une suite de l'étape précédente)
    k <- .vbse(
      paste0("Summary: ", n_sig, "/", p, " dependent variables show significant differences\n",
             "\t(after Bonferroni correction)"),
      paste0("R\u00e9sum\u00e9 : ", n_sig, "/", p, " variables d\u00e9pendantes montrent des diff\u00e9rences significatives\n",
             "\t(apr\u00e8s correction Bonferroni)"),
      verbose = verbose, code = code, k = k, cpt = "off"
    )
  }

  # =============================================================================
  # 2b. COMPARAISONS MULTIVARI\u00c9ES ENTRE GROUPES (Pillai par paire, correction de Holm)
  # =============================================================================
  # R\u00e9f\u00e9rence: Tabachnick & Fidell (2019). Using Multivariate Statistics (7th ed.). Chapter 7.

  if (n_groups > 2) {
    pair_levels <- levels(g)
    pair_p <- matrix(NA_real_, n_groups, n_groups, dimnames = list(pair_levels, pair_levels))
    for (i in seq_len(n_groups - 1)) {
      for (j in (i + 1):n_groups) {
        keep <- g %in% pair_levels[c(i, j)]
        pair_x <- x[keep, , drop = FALSE]
        pair_g <- droplevels(g[keep])
        pair_p[j, i] <- tryCatch(summary(manova(pair_x ~ pair_g), test = "Pillai")$stats[1, "Pr(>F)"],
                                 error = function(e) NA_real_)
      }
    }
    lower <- lower.tri(pair_p)
    pair_p[lower] <- p.adjust(pair_p[lower], method = "holm")
    pair_p <- pair_p[-1, -n_groups, drop = FALSE]
    position <- if (exists("lda_predict")) tapply(lda_predict$x[, 1], g, mean) else tapply(x[, 1], g, mean)
    result$pairwise_multivariate <- pair_p
    result$groups_multivariate <- .letters_from_pairs(pair_p, position, alpha)
    colnames(result$groups_multivariate)[2] <- "Pillai_Holm"
    k <- .vbse(
      "MULTIVARIATE COMPARISONS BETWEEN GROUPS (Pillai's trace, Holm correction):",
      "COMPARAISONS MULTIVARI\u00c9ES ENTRE GROUPES (trace de Pillai, correction de Holm) :",
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    if (verbose) {
      print(result$groups_multivariate, row.names = FALSE)
      cat("\n")
    }
  }

  # =============================================================================
  # 3. RETOUR - Nettoyage des champs internes
  # =============================================================================

  # Conformément à bp.log section 7.4.6 : Supprimer champs internes
  result$k <- NULL         # Compteur .vbse() (usage interne uniquement)
  result$verbose <- NULL   # Paramètres internes
  result$debug <- NULL     # Paramètres internes

  # =============================================================================
  # 4. NORMALISATION AU FORMAT .posthoc() STANDARD
  # =============================================================================
  # Ajouter $groups et $p.value pour compatibilité avec le standard KefiR
  # Référence: Format standard défini dans .posthoc.R

  if (!is.null(result$protected_anovas) && length(result$protected_anovas) > 0) {
    # Créer dataframe $groups avec lettres de significativité
    dv_names <- names(result$protected_anovas)
    sig_status <- sapply(result$protected_anovas, function(x) x$significant)

    # Lettres : 'a' si non significatif, 'b' si significatif (simplifié)
    # Note: Pour MANOVA, les "groupes" sont en fait les DVs
    groups_df <- data.frame(
      dependent_variable = dv_names,
      significance = ifelse(sig_status, "significant", "ns"),
      p_value = sapply(result$protected_anovas, function(x) x$p_value),
      stringsAsFactors = FALSE
    )

    result$groups <- groups_df

    # Créer matrice $p.value (simplifié : p-values des ANOVAs protégées)
    pval_vec <- sapply(result$protected_anovas, function(x) x$p_value)
    names(pval_vec) <- dv_names
    result$p.value <- pval_vec
  } else {
    # Si pas d'ANOVAs protégées, format minimal
    result$groups <- data.frame(
      dependent_variable = character(0),
      significance = character(0),
      stringsAsFactors = FALSE
    )
    result$p.value <- numeric(0)
  }

  # Ajouter classe "posthoc" pour cohérence
  class(result) <- c("posthoc_manova", "posthoc", "list")

  return(result)
}

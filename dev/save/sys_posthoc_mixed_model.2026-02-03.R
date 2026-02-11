################################################################################
#
#                  .posthoc_mixed_model() - Post-hocs pour modèles mixtes
#
################################################################################

#' Post-hoc comparisons for Linear Mixed-Effects Models
#'
#' @description
#' Performs pairwise comparisons on a fitted linear mixed-effects model (lmer)
#' using estimated marginal means (emmeans) with appropriate adjustments.
#' Only executes if at least one fixed effect is significant.
#'
#' @param mixed_model Fitted lmer model object from lme4/lmerTest
#' @param significant_effects Character vector of significant effect names
#'   (from ANOVA Type III). If NULL or length 0, no post-hocs are performed.
#' @param alpha Significance level (default: 0.05)
#' @param method Adjustment method: "tukey" (default), "bonferroni", "holm", "none"
#' @param conf.level Confidence level for intervals (default: 0.95)
#' @param verbose Logical. Print detailed output (default: TRUE)
#' @param debug Logical. Print debug messages (default: FALSE)
#' @param code Logical. Print R code for reproduction (default: FALSE)
#' @param k Message counter for .vbse()
#'
#' @return List containing:
#'   \item{posthoc_results}{List of emmeans results for each significant effect}
#'   \item{compact_letters}{List of compact letter displays (if applicable)}
#'   \item{emmeans_summary}{Summary tables for each effect}
#'   \item{k}{Updated message counter}
#'
#' @details
#' This function uses the emmeans package to compute estimated marginal means
#' (EMMs) and pairwise comparisons for significant fixed effects in a mixed model.
#'
#' **Key features:**
#' - Only performs post-hocs if at least one effect is significant
#' - Uses Satterthwaite approximation for degrees of freedom
#' - Supports multiple comparison adjustment methods
#' - Generates compact letter displays for visualization
#'
#' **Academic rationale:**
#' For mixed models, post-hocs should be performed on estimated marginal means
#' (EMMs) rather than raw means, as EMMs account for the random effects structure
#' and unbalanced designs.
#'
#' @references
#' Lenth, R. V. (2021). emmeans: Estimated Marginal Means, aka Least-Squares Means.
#' R package version 1.7.0. https://CRAN.R-project.org/package=emmeans
#'
#' Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population marginal means
#' in the linear model: An alternative to least squares means. The American
#' Statistician, 34(4), 216-221. https://doi.org/10.1080/00031305.1980.10483031
#'
#' @keywords internal
#' @export
.posthoc_mixed_model <- function(mixed_model,
                                  significant_effects = NULL,
                                  alpha = 0.05,
                                  method = "tukey",
                                  conf.level = 0.95,
                                  verbose = TRUE,
                                  debug = FALSE,
                                  code = FALSE,
                                  k = NULL) {

  # Initialiser k
  if (is.null(k)) k <- 0

  #=============================================================================
  # BLOC 1 : VÉRIFICATIONS PRÉALABLES
  #=============================================================================

  # Vérifier package emmeans
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    k <- .vbse(
      paste0("Note: Package 'emmeans' not available. Post-hoc tests for mixed models skipped.\n",
             "\tInstall with: install.packages('emmeans')"),
      paste0("Note : Package 'emmeans' non disponible. Tests post-hoc pour modèles mixtes ignorés.\n",
             "\tInstaller avec : install.packages('emmeans')"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    return(list(posthoc_results = NULL, k = k))
  }

  # Vérifier que le modèle est valide
  if (is.null(mixed_model) || !inherits(mixed_model, "lmerMod")) {
    k <- .vbse(
      "Error: Invalid mixed model object. Post-hocs skipped.",
      "Erreur : Objet modèle mixte invalide. Post-hocs ignorés.",
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    return(list(posthoc_results = NULL, k = k))
  }

  # Vérifier si au moins un effet significatif
  if (is.null(significant_effects) || length(significant_effects) == 0) {
    k <- .vbse(
      "No significant fixed effects detected. Post-hoc tests not performed for mixed model.",
      "Aucun effet fixe significatif détecté. Tests post-hoc non effectués pour modèle mixte.",
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    return(list(posthoc_results = NULL, k = k))
  }

  #=============================================================================
  # BLOC 2 : HEADER POST-HOCS
  #=============================================================================

  k <- .vbse(
    paste0("Post-hoc pairwise comparisons for mixed model (emmeans method):\n",
           "\tSignificant effect(s): ", paste(significant_effects, collapse = ", "), "\n",
           "\tAdjustment method: ", method, "\n",
           "\tConfidence level: ", conf.level),
    paste0("Comparaisons post-hoc par paires pour modèle mixte (méthode emmeans) :\n",
           "\tEffet(s) significatif(s) : ", paste(significant_effects, collapse = ", "), "\n",
           "\tMéthode d'ajustement : ", method, "\n",
           "\tNiveau de confiance : ", conf.level),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  #=============================================================================
  # BLOC 3 : MODE CODE
  #=============================================================================

  if (isTRUE(code)) {
    cat("# ============================================================\n")
    cat("# Comparaisons post-hoc pour modèle mixte (emmeans)\n")
    cat("# ============================================================\n\n")
    cat("library(emmeans)\n")
    cat("library(lmerTest)\n\n")
    cat("# Ajuster le modèle mixte (déjà fait)\n")
    cat("# mixed_model <- lmer(formule, data = data)\n\n")
    cat("# Effectuer les comparaisons post-hoc pour chaque effet significatif\n")
    for (effect in significant_effects) {
      cat(paste0("emm_", effect, " <- emmeans(mixed_model, specs = ~ ", effect, ")\n"))
      cat(paste0("pairs_", effect, " <- pairs(emm_", effect, ", adjust = '", method, "')\n"))
      cat(paste0("summary(pairs_", effect, ")\n\n"))
    }
    cat("# ============================================================\n\n")
  }

  #=============================================================================
  # BLOC 4 : CALCUL POST-HOCS POUR CHAQUE EFFET SIGNIFICATIF
  #=============================================================================

  posthoc_results <- list()
  emmeans_summary <- list()
  compact_letters <- list()

  for (effect in significant_effects) {

    .dbg(paste0("Computing post-hocs for effect: ", effect),
         paste0("Calcul post-hocs pour effet : ", effect),
         debug = debug)

    # Calculer EMMs (suppressMessages pour masquer avertissements emmeans en anglais)
    emm_result <- tryCatch({
      suppressMessages(
        emmeans::emmeans(mixed_model, specs = as.formula(paste0("~ ", effect)))
      )
    }, error = function(e) {
      k <<- .vbse(
        paste0("Warning: Could not compute emmeans for effect '", effect, "': ", e$message),
        paste0("Attention : Impossible de calculer emmeans pour effet '", effect, "' : ", e$message),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      return(NULL)
    })

    if (is.null(emm_result)) next

    # Comparaisons par paires
    # Note: pairs() est une méthode S3, il faut charger emmeans et appeler sans namespace
    pairs_result <- tryCatch({
      if (requireNamespace("emmeans", quietly = TRUE)) {
        # Utiliser contrast() qui est exportée, plus robuste que pairs()
        emmeans::contrast(emm_result, method = "pairwise", adjust = method)
      } else {
        stop("emmeans package required")
      }
    }, error = function(e) {
      k <<- .vbse(
        paste0("Warning: Could not compute pairwise comparisons for '", effect, "': ", e$message),
        paste0("Attention : Impossible de calculer comparaisons par paires pour '", effect, "' : ", e$message),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      return(NULL)
    })

    if (is.null(pairs_result)) next

    # Compact letter display
    cld_result <- tryCatch({
      if (requireNamespace("multcomp", quietly = TRUE)) {
        multcomp::cld(emm_result, adjust = method, Letters = letters, alpha = alpha)
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })

    # Stocker résultats
    posthoc_results[[effect]] <- list(
      emmeans = emm_result,
      pairs = pairs_result,
      cld = cld_result
    )

    emmeans_summary[[effect]] <- summary(pairs_result)

    # Afficher résultats
    pairs_summary <- summary(pairs_result)
    k <- .vbse(
      paste0("\tEffect: ", effect),
      paste0("\tEffet : ", effect),
      verbose = verbose, code = code, k = k, cpt = "off"
    )

    for (i in 1:nrow(pairs_summary)) {
      comparison <- as.character(pairs_summary[i, "contrast"])
      estimate <- pairs_summary[i, "estimate"]
      se <- pairs_summary[i, "SE"]
      p_val <- pairs_summary[i, "p.value"]

      sig_mark <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"

      k <- .vbse(
        paste0("\t  ", comparison, ": diff = ", round(estimate, 3),
               " (SE = ", round(se, 3), "), p = ", .format_pval(p_val), " ", sig_mark),
        paste0("\t  ", comparison, " : diff = ", round(estimate, 3),
               " (SE = ", round(se, 3), "), p = ", .format_pval(p_val), " ", sig_mark),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }

    # Afficher CLD si disponible
    if (!is.null(cld_result)) {
      cld_df <- tryCatch({
        # Convertir CLD en data.frame de manière robuste
        df <- as.data.frame(cld_result)
        # S'assurer que les noms de lignes sont dans une colonne
        if (is.null(df[[1]])) {
          df <- cbind(level = rownames(df), df)
        }
        df
      }, error = function(e) {
        # Si la conversion échoue, essayer une approche alternative
        tryCatch({
          data.frame(
            level = names(cld_result),
            emmean = as.numeric(cld_result),
            .group = attr(cld_result, ".group")
          )
        }, error = function(e2) {
          k <<- .vbse(
            paste0("\tWarning: Could not display compact letter display: ", e2$message),
            paste0("\tAttention : Impossible d'afficher les lettres groupes : ", e2$message),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
          return(NULL)
        })
      })

      if (!is.null(cld_df) && nrow(cld_df) > 0) {
        k <- .vbse(
          paste0("\tCompact letter display (groups not sharing a letter differ at alpha = ", alpha, "):"),
          paste0("\tAffichage compact par lettres (groupes ne partageant pas une lettre diffèrent à alpha = ", alpha, ") :"),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        for (i in 1:nrow(cld_df)) {
          level_name <- as.character(cld_df[i, 1])
          emmean <- round(as.numeric(cld_df[i, "emmean"]), 3)
          group_col <- which(colnames(cld_df) == ".group")
          if (length(group_col) > 0) {
            group <- trimws(as.character(cld_df[i, group_col]))
          } else {
            group <- "?"
          }

          k <- .vbse(
            paste0("\t  ", level_name, ": EMM = ", emmean, ", group = ", group),
            paste0("\t  ", level_name, " : EMM = ", emmean, ", groupe = ", group),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }

        compact_letters[[effect]] <- cld_df
      }
    }
  }

  #=============================================================================
  # BLOC 5 : CONSTRUCTION FORMAT STANDARD (.posthoc compatible)
  #=============================================================================

  # Si un seul effet significatif : format simple compatible avec .posthoc()
  # Si plusieurs effets : format étendu avec sous-listes par effet

  if (length(significant_effects) == 1) {
    # CAS SIMPLE : Un seul effet - Format compatible .posthoc()
    effect <- significant_effects[1]

    # Extraire les lettres CLD si disponibles
    groups_df <- NULL
    if (!is.null(compact_letters[[effect]])) {
      cld_data <- compact_letters[[effect]]
      # Construire data.frame au format standard
      groups_df <- data.frame(
        categories = as.character(cld_data[[1]]),
        EMMs_Tukey = trimws(as.character(cld_data[[".group"]])),
        stringsAsFactors = FALSE
      )
    } else if (!is.null(emmeans_summary[[effect]])) {
      # Fallback: juste les catégories sans lettres
      pairs_df <- emmeans_summary[[effect]]
      # Extraire niveaux uniques du facteur
      contrasts <- as.character(pairs_df$contrast)
      levels_list <- unique(unlist(strsplit(contrasts, " - ")))
      groups_df <- data.frame(
        categories = levels_list,
        EMMs_Tukey = "",
        stringsAsFactors = FALSE
      )
    }

    synth <- list(
      groups = groups_df,
      method = "EMMs (emmeans package)",
      effect = effect,
      pairs = emmeans_summary[[effect]],
      k = k
    )

    # Ajouter classe standard
    class(synth) <- "posthoc"

  } else {
    # CAS MULTI-EFFETS : Format étendu avec sous-structure par effet
    synth <- list()

    for (effect in significant_effects) {
      # Extraire les lettres CLD si disponibles
      groups_df <- NULL
      if (!is.null(compact_letters[[effect]])) {
        cld_data <- compact_letters[[effect]]
        groups_df <- data.frame(
          categories = as.character(cld_data[[1]]),
          EMMs_Tukey = trimws(as.character(cld_data[[".group"]])),
          stringsAsFactors = FALSE
        )
      } else if (!is.null(emmeans_summary[[effect]])) {
        # Fallback
        pairs_df <- emmeans_summary[[effect]]
        contrasts <- as.character(pairs_df$contrast)
        levels_list <- unique(unlist(strsplit(contrasts, " - ")))
        groups_df <- data.frame(
          categories = levels_list,
          EMMs_Tukey = "",
          stringsAsFactors = FALSE
        )
      }

      synth[[effect]] <- list(
        groups = groups_df,
        method = "EMMs (emmeans package)",
        pairs = emmeans_summary[[effect]]
      )
      class(synth[[effect]]) <- "posthoc"
    }

    synth$note <- paste0("Multiple effects: ", paste(significant_effects, collapse = ", "))
  }

  return(synth)
}

################################################################################
# FIN .posthoc_mixed_model()
################################################################################

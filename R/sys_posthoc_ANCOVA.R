#' Post-hoc comparisons for ANCOVA (Analysis of Covariance)
#'
#' @description
#' Performs appropriate post-hoc pairwise comparisons following a significant ANCOVA,
#' using estimated marginal means (EMMs) adjusted for covariates. This function follows
#' academic best practices by:
#' \itemize{
#'   \item Computing adjusted means at covariate mean values
#'   \item Performing pairwise comparisons of adjusted means
#'   \item Applying appropriate corrections for multiple comparisons
#'   \item Providing compact letter display for group differences
#' }
#'
#' @details
#' \strong{Academic Rationale:}
#'
#' After a significant ANCOVA, post-hoc comparisons should NOT be performed on raw group means,
#' but on \strong{adjusted means} (also called least-squares means or estimated marginal means).
#' These adjusted means represent the expected response for each group when all covariates
#' are held at their mean values.
#'
#' \strong{Why adjusted means?}
#' \enumerate{
#'   \item Raw group means are confounded by covariate differences between groups
#'   \item ANCOVA removes covariate effects; post-hocs should reflect this adjustment
#'   \item Adjusted means provide fair comparisons controlling for covariate influence
#' }
#'
#' \strong{Recommended Methods:}
#' \itemize{
#'   \item \strong{Tukey HSD} (default): Controls familywise error rate, assumes equal variances
#'   \item \strong{Bonferroni}: More conservative, suitable for few comparisons
#'   \item \strong{Holm}: Sequential Bonferroni, more powerful
#'   \item \strong{Sidak}: Similar to Bonferroni but less conservative
#' }
#'
#' \strong{When homoscedasticity is violated:}
#' Games-Howell test on residuals (experimental feature, use with caution).
#'
#' @param ancova_model Object returned by \code{aov()} or \code{lm()} with ANCOVA formula.
#'   Must include at least one factor and one numeric covariate.
#' @param factor_name Character. Name of the grouping factor for which to perform comparisons.
#'   If NULL (default), uses the first factor in the model.
#' @param alpha Numeric. Significance level (default: 0.05).
#' @param method Character. Adjustment method for multiple comparisons:
#'   \code{"tukey"} (default), \code{"bonferroni"}, \code{"holm"}, \code{"sidak"}, \code{"none"}.
#' @param conf.level Numeric. Confidence level for confidence intervals (default: 0.95).
#' @param verbose Logical. If TRUE, display pedagogical messages (default: TRUE).
#' @param debug Logical. If TRUE, display debugging messages (default: FALSE).
#' @param code Logical. If TRUE, display reproducible R code (default: FALSE).
#' @param k Integer or NULL. Message counter for .vbse() (internal use).
#'
#' @return
#' A list of class \code{"posthoc"} containing:
#' \itemize{
#'   \item \code{groups}: data.frame with group names and compact letter display
#'   \item \code{adjusted_means}: data.frame with adjusted means and standard errors
#'   \item \code{pairwise_comparisons}: data.frame with all pairwise comparisons
#'   \item \code{p.value}: matrix of pairwise p-values
#'   \item \code{method}: character string describing the adjustment method used
#'   \item \code{note}: additional notes or warnings
#' }
#'
#' @references
#' \itemize{
#'   \item Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018). \emph{Designing Experiments
#'     and Analyzing Data: A Model Comparison Perspective} (3rd ed.). Routledge. Chapter 9 (ANCOVA).
#'     ISBN: 978-1138892286
#'   \item Huitema, B. E. (2011). \emph{The Analysis of Covariance and Alternatives:
#'     Statistical Methods for Experiments, Quasi-Experiments, and Single-Case Studies}
#'     (2nd ed.). Wiley. DOI: 10.1002/9781118067475
#'   \item Rutherford, A. (2011). \emph{ANOVA and ANCOVA: A GLM Approach} (2nd ed.). Wiley.
#'     ISBN: 978-0470385555
#'   \item Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population marginal means
#'     in the linear model: An alternative to least squares means. \emph{The American
#'     Statistician}, 34(4), 216-221. DOI: 10.1080/00031305.1980.10483031
#' }
#'
#' @importFrom stats coef model.matrix pairwise.t.test p.adjust terms
#' @importFrom emmeans emmeans contrast
#'
#' @seealso
#' \code{\link[emmeans]{emmeans}}, \code{\link[emmeans]{contrast}},
#' \code{\link{.posthoc}}, \code{\link{.posthoc_MANOVA}}
#'
#' @examples
#' \dontrun{
#' # Example: Effect of treatment on weight, controlling for baseline weight
#' set.seed(123)
#' n <- 60
#' baseline_weight <- rnorm(n, 70, 10)
#' treatment <- factor(rep(c("Control", "Drug_A", "Drug_B"), each = 20))
#' # Treatment effect + covariate effect
#' final_weight <- 50 + 0.6 * baseline_weight +
#'                 5 * (treatment == "Drug_A") +
#'                 10 * (treatment == "Drug_B") +
#'                 rnorm(n, 0, 3)
#'
#' # ANCOVA
#' ancova_model <- aov(final_weight ~ treatment + baseline_weight)
#'
#' # Post-hoc comparisons on adjusted means
#' posthoc_results <- .posthoc_ANCOVA(ancova_model,
#'                                    factor_name = "treatment",
#'                                    alpha = 0.05,
#'                                    method = "tukey",
#'                                    verbose = TRUE)
#' print(posthoc_results)
#' }
#'
#' @keywords internal
#' @export
.posthoc_ANCOVA <- function(ancova_model,
                            factor_name = NULL,
                            alpha = 0.05,
                            method = "tukey",
                            conf.level = 0.95,
                            verbose = TRUE,
                            debug = FALSE,
                            code = FALSE,
                            k = NULL) {

  # NOTE PERSO: Cette fonction implémente les post-hocs appropriés pour ANCOVA
  # selon les recommandations académiques (Maxwell & Delaney, 2018; Huitema, 2011)
  # Format de sortie compatible avec .posthoc() pour intégration dans m.test()

  # ============================================================================
  # BLOC 1: VALIDATION ET INITIALISATION
  # ============================================================================

  .dbg("=== Start .posthoc_ANCOVA() ===",
       "=== Début .posthoc_ANCOVA() ===", debug = debug)

  if (is.null(k)) k <- 0

  # Vérifier que le modèle est bien un objet aov ou lm
  if (!inherits(ancova_model, c("aov", "lm"))) {
    .exit("ancova_model must be an object of class 'aov' or 'lm'.",
          "ancova_model doit être un objet de classe 'aov' ou 'lm'.",
          verbose = verbose, code = code, return = TRUE)
  }

  # Vérifier disponibilité du package emmeans
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    k <- .vbse(
      paste0("ERROR: Package 'emmeans' is required for ANCOVA post-hocs.\n",
             "Please install it with: install.packages('emmeans')\n\n",
             "Academic justification:\n",
             "  ANCOVA post-hocs require adjusted means (EMMs), which are computed\n",
             "  by the 'emmeans' package following best statistical practices.\n",
             "  Raw group means are inappropriate after ANCOVA because they are\n",
             "  confounded by covariate differences between groups.\n\n",
             "References:\n",
             "  Lenth, R. V. (2024). emmeans: Estimated Marginal Means. R package.\n",
             "  https://CRAN.R-project.org/package=emmeans"),
      paste0("ERREUR : Le package 'emmeans' est requis pour les post-hocs ANCOVA.\n",
             "Veuillez l'installer avec : install.packages('emmeans')\n\n",
             "Justification académique :\n",
             "  Les post-hocs ANCOVA nécessitent des moyennes ajustées (EMMs), calculées\n",
             "  par le package 'emmeans' selon les meilleures pratiques statistiques.\n",
             "  Les moyennes brutes sont inappropriées après ANCOVA car elles sont\n",
             "  confondues par les différences de covariables entre groupes.\n\n",
             "Références :\n",
             "  Lenth, R. V. (2024). emmeans: Estimated Marginal Means. R package.\n",
             "  https://CRAN.R-project.org/package=emmeans"),
      verbose = TRUE, k = k, cpt = "on"
    )

    # Retourner objet posthoc minimal
    synth <- list(
      groups = data.frame(
        categories = character(0),
        ANCOVA = character(0)
      ),
      note = "Package 'emmeans' required but not available.",
      method = "none"
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  # Extraire termes du modèle
  model_terms <- terms(ancova_model)
  term_labels <- attr(model_terms, "term.labels")

  # Identifier facteurs et covariables
  model_data <- model.frame(ancova_model)
  factor_vars <- character(0)
  covariate_vars <- character(0)

  for (term in term_labels) {
    # Vérifier si c'est un facteur ou une covariable
    if (term %in% names(model_data)) {
      # IMPORTANT: aov() traite automatiquement character comme factor
      # Donc on doit détecter BOTH factor ET character comme facteurs catégoriques
      if (is.factor(model_data[[term]]) || is.character(model_data[[term]])) {
        factor_vars <- c(factor_vars, term)
      } else if (is.numeric(model_data[[term]])) {
        covariate_vars <- c(covariate_vars, term)
      }
    }
  }

  .dbg(paste0("Factors found: ", paste(factor_vars, collapse = ", ")),
       paste0("Facteurs trouvés : ", paste(factor_vars, collapse = ", ")),
       debug = debug)

  .dbg(paste0("Covariates found: ", paste(covariate_vars, collapse = ", ")),
       paste0("Covariables trouvées : ", paste(covariate_vars, collapse = ", ")),
       debug = debug)

  # Vérifier qu'il y a au moins 1 facteur et 1 covariable
  if (length(factor_vars) == 0) {
    k <- .vbse(
      "ERROR: No factor found in ANCOVA model. At least one factor is required.",
      "ERREUR : Aucun facteur trouvé dans le modèle ANCOVA. Au moins un facteur est requis.",
      verbose = TRUE, k = k, cpt = "on"
    )

    synth <- list(
      groups = data.frame(categories = character(0), ANCOVA = character(0)),
      note = "No factor found in model.",
      method = "none"
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  if (length(covariate_vars) == 0) {
    k <- .vbse(
      paste0("WARNING: No numeric covariate found in model.\n",
             "Are you sure this is an ANCOVA? Consider using .posthoc() instead."),
      paste0("ATTENTION : Aucune covariable numérique trouvée dans le modèle.\n",
             "Êtes-vous sûr qu'il s'agit d'une ANCOVA ? Considérez .posthoc() à la place."),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
  }

  # Déterminer le facteur à analyser
  if (is.null(factor_name)) {
    factor_name <- factor_vars[1]
    .dbg(paste0("factor_name not specified. Using first factor: ", factor_name),
         paste0("factor_name non spécifié. Utilisation du premier facteur : ", factor_name),
         debug = debug)
  }

  if (!factor_name %in% factor_vars) {
    k <- .vbse(
      paste0("ERROR: factor_name '", factor_name, "' not found in model factors: ",
             paste(factor_vars, collapse = ", ")),
      paste0("ERREUR : factor_name '", factor_name, "' non trouvé dans les facteurs du modèle : ",
             paste(factor_vars, collapse = ", ")),
      verbose = TRUE, k = k, cpt = "on"
    )

    synth <- list(
      groups = data.frame(categories = character(0), ANCOVA = character(0)),
      note = paste0("Factor '", factor_name, "' not found."),
      method = "none"
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  # Extraire les niveaux du facteur
  factor_data <- model_data[[factor_name]]
  # IMPORTANT: Si character, le convertir en factor pour obtenir les levels
  if (is.character(factor_data)) {
    factor_data <- as.factor(factor_data)
  }
  factor_levels <- levels(factor_data)
  n_groups <- length(factor_levels)

  .dbg(paste0("Factor '", factor_name, "' has ", n_groups, " levels: ",
              paste(factor_levels, collapse = ", ")),
       paste0("Le facteur '", factor_name, "' a ", n_groups, " niveaux : ",
              paste(factor_levels, collapse = ", ")),
       debug = debug)

  if (n_groups < 2) {
    synth <- list(
      groups = data.frame(categories = factor_levels, ANCOVA = "a"),
      note = "Only one group, no post-hoc needed.",
      method = "none"
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  # ============================================================================
  # BLOC 2: MESSAGES PÉDAGOGIQUES
  # ============================================================================

  k <- .vbse(
    paste0("POST-HOC COMPARISONS FOR ANCOVA [{emmeans}]\n",
           "\tANCOVA post-hocs compare ADJUSTED MEANS, not raw means.\n",
           "\tAdjusted means = expected response when covariates are at their mean.\n",
           "\tThis controls for covariate effects and provides fair group comparisons.\n\n",
           "\tModel: ", deparse(formula(ancova_model)), "\n",
           "\tFactor analyzed: ", factor_name, " (", n_groups, " groups)\n",
           "\tCovariates: ", paste(covariate_vars, collapse = ", "), "\n",
           "\tAdjustment method: ", method, "\n",
           "\tSignificance level: α = ", alpha),
    paste0("COMPARAISONS POST-HOC POUR ANCOVA [{emmeans}]\n",
           "\tLes post-hocs ANCOVA comparent les MOYENNES AJUSTÉES, pas les moyennes brutes.\n",
           "\tMoyennes ajustées = réponse attendue quand les covariables sont à leur moyenne.\n",
           "\tCela contrôle les effets des covariables et fournit des comparaisons équitables.\n\n",
           "\tModèle : ", deparse(formula(ancova_model)), "\n",
           "\tFacteur analysé : ", factor_name, " (", n_groups, " groupes)\n",
           "\tCovariables : ", paste(covariate_vars, collapse = ", "), "\n",
           "\tMéthode d'ajustement : ", method, "\n",
           "\tSeuil de significativité : α = ", alpha),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  # ============================================================================
  # BLOC 3: CALCUL DES MOYENNES AJUSTÉES (EMMs) - Fusionné avec affichage
  # ============================================================================

  if (isTRUE(code)) {
    cat("# ==================================================\n")
    cat("# Post-hocs ANCOVA : Moyennes marginales estimées\n")
    cat("# ==================================================\n\n")
    cat("# Calculer les moyennes ajustées (EMMs) aux moyennes des covariables\n")
    cat("library(emmeans)\n")
    cat("emm_result <- emmeans(ancova_model, specs = \"", factor_name, "\")\n", sep = "")
    cat("summary(emm_result)\n\n")
  }

  .dbg("Computing estimated marginal means...",
       "Calcul des moyennes marginales estimées...", debug = debug)

  # Supprimer le message anglais "NOTE: Results may be misleading..."
  emm_result <- tryCatch({
    suppressMessages(emmeans::emmeans(ancova_model, specs = factor_name))
  }, error = function(e) {
    k <<- .vbse(
      paste0("ERROR computing EMMs: ", e$message),
      paste0("ERREUR lors du calcul des EMMs : ", e$message),
      verbose = TRUE, k = k, cpt = "on"
    )
    return(NULL)
  })

  if (is.null(emm_result)) {
    synth <- list(
      groups = data.frame(categories = factor_levels, ANCOVA = rep("?", n_groups)),
      note = "Error computing adjusted means.",
      method = method
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  # Extraire moyennes ajustées - AFFICHAGE FUSIONNÉ avec étape précédente
  emm_summary <- summary(emm_result)
  adjusted_means <- data.frame(
    group = as.character(emm_summary[[factor_name]]),
    adjusted_mean = emm_summary$emmean,
    SE = emm_summary$SE,
    lower.CL = emm_summary$lower.CL,
    upper.CL = emm_summary$upper.CL
  )

  # Affichage intégré avec ajout "Moyennes ajustées:"
  if (verbose) {
    cat("\n\tMoyennes ajustées (aux moyennes des covariables) :\n")
    print(adjusted_means, row.names = FALSE)
    cat("\n")
  }

  # ============================================================================
  # BLOC 4: COMPARAISONS PAR PAIRES
  # ============================================================================

  if (isTRUE(code)) {
    cat("# Comparaisons par paires des moyennes ajustées\n")
    cat("pairs_result <- pairs(emm_result, adjust = \"", method, "\")\n", sep = "")
    cat("summary(pairs_result)\n\n")
  }

  .dbg(paste0("Performing pairwise comparisons with method: ", method),
       paste0("Réalisation des comparaisons par paires avec méthode : ", method),
       debug = debug)

  # Convertir method en format emmeans
  emmeans_method <- switch(method,
                            "tukey" = "tukey",
                            "bonferroni" = "bonferroni",
                            "holm" = "holm",
                            "sidak" = "sidak",
                            "none" = "none",
                            "tukey")  # default

  pairs_result <- tryCatch({
    emmeans::contrast(emm_result, method = "pairwise", adjust = emmeans_method)
  }, error = function(e) {
    k <<- .vbse(
      paste0("ERROR computing pairwise comparisons: ", e$message),
      paste0("ERREUR lors des comparaisons par paires : ", e$message),
      verbose = TRUE, k = k, cpt = "on"
    )
    return(NULL)
  })

  if (is.null(pairs_result)) {
    synth <- list(
      groups = data.frame(categories = factor_levels, ANCOVA = rep("?", n_groups)),
      adjusted_means = adjusted_means,
      note = "Error computing pairwise comparisons.",
      method = method
    )
    class(synth) <- "posthoc"
    return(synth)
  }

  # Extraire résultats des comparaisons
  pairs_summary <- summary(pairs_result)
  pairwise_comparisons <- data.frame(
    contrast = as.character(pairs_summary$contrast),
    estimate = pairs_summary$estimate,
    SE = pairs_summary$SE,
    df = pairs_summary$df,
    t_ratio = pairs_summary$t.ratio,
    p_value = pairs_summary$p.value
  )

  # Affichage avec analyse des différences significatives
  if (verbose) {
    cat("\n\tComparaisons par paires (ajustement ", emmeans_method, ") :\n", sep = "")
    print(pairwise_comparisons, row.names = FALSE)

    # Conclusion : identifier les différences significatives
    sig_comparisons <- pairwise_comparisons[pairwise_comparisons$p_value < alpha, ]
    if (nrow(sig_comparisons) > 0) {
      cat("\n\t=> Différences significatives détectées :\n")
      for (i in 1:nrow(sig_comparisons)) {
        cat("\t   ", as.character(sig_comparisons$contrast[i]),
            " (p = ", .format_pval(sig_comparisons$p_value[i]), ")\n", sep = "")
      }
    } else {
      cat("\n\t=> Aucune différence significative entre les groupes.\n")
    }
    cat("\n")
  }

  # ============================================================================
  # BLOC 5: COMPACT LETTER DISPLAY (CLD)
  # ============================================================================

  .dbg("Computing compact letter display...",
       "Calcul de l'affichage compact par lettres...", debug = debug)

  # Utiliser emmeans::cld() pour générer compact letter display
  cld_result <- tryCatch({
    emmeans::cld(emm_result, alpha = alpha, Letters = letters, adjust = emmeans_method)
  }, error = function(e) {
    .dbg(paste0("Warning: Could not compute CLD: ", e$message),
         paste0("Attention : Impossible de calculer le CLD : ", e$message),
         debug = debug)
    return(NULL)
  })

  if (!is.null(cld_result)) {
    # Extraire lettres
    groups_df <- data.frame(
      categories = as.character(cld_result[[factor_name]]),
      ANCOVA = trimws(as.character(cld_result$.group))
    )
  } else {
    # Fallback: générer lettres basées sur p-values
    .dbg("Using fallback letter generation based on p-values",
         "Utilisation de la génération de lettres de secours basée sur les p-values",
         debug = debug)

    # Créer matrice de p-values
    p_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
    rownames(p_matrix) <- factor_levels
    colnames(p_matrix) <- factor_levels

    for (i in 1:nrow(pairwise_comparisons)) {
      contrast_str <- as.character(pairwise_comparisons$contrast[i])
      # Parser le contraste (format: "level1 - level2")
      parts <- strsplit(contrast_str, " - ")[[1]]
      if (length(parts) == 2) {
        g1 <- trimws(parts[1])
        g2 <- trimws(parts[2])
        if (g1 %in% factor_levels && g2 %in% factor_levels) {
          p_val <- pairwise_comparisons$p_value[i]
          p_matrix[g1, g2] <- p_val
          p_matrix[g2, g1] <- p_val
        }
      }
    }

    # Générer lettres simples
    letters_assigned <- rep("a", n_groups)
    names(letters_assigned) <- factor_levels

    # Algorithme simple de lettres basé sur alpha
    current_letter <- 1
    for (i in 1:n_groups) {
      if (letters_assigned[i] == "a") {
        same_group <- c(i)
        for (j in (i+1):n_groups) {
          if (j <= n_groups && p_matrix[i, j] > alpha) {
            same_group <- c(same_group, j)
          }
        }
        for (idx in same_group) {
          letters_assigned[idx] <- letters[current_letter]
        }
        current_letter <- current_letter + 1
      }
    }

    groups_df <- data.frame(
      categories = factor_levels,
      ANCOVA = letters_assigned
    )
  }

  # ÉTAPE 13 SUPPRIMÉE : L'affichage compact par lettres est maintenant
  # seulement disponible via la sortie de m.test() (synth$groups)
  # pour éviter redondance dans la sortie verbose

  # ============================================================================
  # BLOC 6: CRÉATION MATRICE P-VALUES
  # ============================================================================

  # Créer matrice de p-values pour compatibilité avec .posthoc()
  p_value_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
  rownames(p_value_matrix) <- factor_levels
  colnames(p_value_matrix) <- factor_levels

  for (i in 1:nrow(pairwise_comparisons)) {
    contrast_str <- as.character(pairwise_comparisons$contrast[i])
    parts <- strsplit(contrast_str, " - ")[[1]]
    if (length(parts) == 2) {
      g1 <- trimws(parts[1])
      g2 <- trimws(parts[2])
      if (g1 %in% factor_levels && g2 %in% factor_levels) {
        p_val <- pairwise_comparisons$p_value[i]
        p_value_matrix[g1, g2] <- p_val
        p_value_matrix[g2, g1] <- p_val
      }
    }
  }

  # ============================================================================
  # BLOC 7: CONSTRUCTION OBJET DE RETOUR
  # ============================================================================

  # Structure STANDARD harmonisée (même ordre que ANOVA)
  # Éléments obligatoires : groups, p.value, global_pvalue, descriptive_stats, method
  # Éléments spécifiques ANCOVA stockés dans $details (pour utilisateurs avancés)

  synth <- list(
    # --- Éléments STANDARD (harmonisés avec ANOVA) ---
    groups = groups_df,
    p.value = p_value_matrix,
    # global_pvalue sera ajouté par m.test()
    # descriptive_stats sera ajouté par m.test() depuis adjusted_means
    method = paste0("ANCOVA post-hoc (", method, " adjustment)"),

    # --- Éléments spécifiques ANCOVA (pour compatibilité/détails) ---
    adjusted_means = adjusted_means,  # Utilisé pour construire descriptive_stats
    details = list(
      pairwise_comparisons = pairwise_comparisons,
      emm_result = emm_result,
      pairs_result = pairs_result,
      note = paste0("Post-hoc comparisons performed on adjusted means. ",
                    "Covariates held at their mean values: ",
                    paste(paste0(covariate_vars, " = ",
                                 round(sapply(covariate_vars, function(cv)
                                   mean(model_data[[cv]], na.rm = TRUE)), 2)),
                          collapse = ", "))
    )
  )

  class(synth) <- "posthoc"

  .dbg("=== End .posthoc_ANCOVA() ===",
       "=== Fin .posthoc_ANCOVA() ===", debug = debug)

  return(synth)
}

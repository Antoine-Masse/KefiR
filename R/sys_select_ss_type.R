#' Sélectionne le type de Sommes des Carrés optimal (Type I/II/III)
#'
#' @description
#' Fonction helper qui détermine le type de sommes des carrés le plus approprié
#' selon le contexte (ANOVA vs ANCOVA, équilibré vs déséquilibré, ordre formule).
#'
#' @param formula Formule du modèle
#' @param data Data frame
#' @param model Modèle lm() ajusté (optionnel, pour tester interactions)
#' @param analysis_type Character. "ANOVA" ou "ANCOVA"
#' @param alpha Seuil de significativité pour tester interactions (défaut: 0.05)
#' @param verbose Logical. Afficher messages détaillés (défaut: FALSE)
#' @param k Compteur d'étapes (défaut: 0)
#' @param debug Logical. Mode debug (défaut: FALSE)
#'
#' @return Liste contenant:
#' \itemize{
#'   \item \code{ss_type}: "I", "II", ou "III"
#'   \item \code{reason}: Raison du choix (character)
#'   \item \code{balanced}: Logical. Plan équilibré?
#'   \item \code{formula_optimal}: Logical. Formule optimale? (pour ANCOVA)
#'   \item \code{suggested_formula}: Formule suggérée si non-optimale (ou NULL)
#'   \item \code{warnings}: Character vector de warnings
#'   \item \code{k}: Compteur mis à jour
#' }
#'
#' @details
#' **RÈGLES DE SÉLECTION**:
#'
#' **ANOVA à deux facteurs (sans covariables)**:
#' 1. Plan équilibré → Type I par défaut
#' 2. Déséquilibré + interaction non-significative → Type II
#' 3. Déséquilibré + interaction significative → Type III
#'
#' **ANCOVA (n facteurs + m covariables)**:
#' 1. Plan équilibré + covariables en premier + bon ordre → Type I
#' 2. Déséquilibre → Type III obligatoire
#' 3. Covariables mal placées → Warning + Type III
#' 4. Covariables mal ordonnées → Warning + Type III
#' 5. Interactions facteur×covariable → Pipeline pentes hétérogènes
#'
#' @references
#' - Maxwell, S. E., & Delaney, H. D. (2004). Designing Experiments and
#'   Analyzing Data (2nd ed.). Psychology Press.
#' - Huitema, B. (2011). The Analysis of Covariance and Alternatives (2nd ed.).
#'   Wiley.
#'
#' @keywords internal
.select_ss_type <- function(formula, data, model = NULL,
                            analysis_type = c("ANOVA", "ANCOVA"),
                            alpha = 0.05, verbose = FALSE, k = 0, debug = FALSE, code = NULL) {

  analysis_type <- match.arg(analysis_type)

  # ==========================================================================
  # EXTRACTION VARIABLES DE LA FORMULE
  # ==========================================================================

  all_vars <- all.vars(formula)
  response_var <- all_vars[1]
  predictor_vars <- all_vars[-1]

  # Identifier facteurs vs covariables
  factor_vars <- character(0)
  numeric_vars <- character(0)

  for (var in predictor_vars) {
    if (is.factor(data[[var]]) || is.character(data[[var]])) {
      factor_vars <- c(factor_vars, var)
    } else if (is.numeric(data[[var]])) {
      numeric_vars <- c(numeric_vars, var)
    }
  }

  n_factors <- length(factor_vars)
  n_covariates <- length(numeric_vars)

  .dbg(paste0("SS Type Selection: ", n_factors, " factors, ", n_covariates, " covariates"),
       paste0("S\u00e9lection Type SS: ", n_factors, " facteurs, ", n_covariates, " covariables"),
       debug = debug)

  # ==========================================================================
  # D\u00c9TECTION \u00c9QUILIBRE DU PLAN FACTORIEL
  # ==========================================================================

  balanced <- FALSE
  if (n_factors > 0) {
    # Cr\u00e9er interaction compl\u00e8te des facteurs
    if (n_factors == 1) {
      g_cat <- data[[factor_vars[1]]]
    } else {
      g_cat <- interaction(data[, factor_vars, drop = FALSE], drop = TRUE)
    }

    # Compter par groupe
    table_counts <- table(g_cat)

    # \u00c9quilibr\u00e9 si toutes les cellules ont m\u00eame taille
    balanced <- (length(unique(table_counts)) == 1)

    .dbg(paste0("Plan balanced: ", balanced, " (counts: ",
                paste(table_counts, collapse = ", "), ")"),
         paste0("Plan \u00e9quilibr\u00e9: ", balanced, " (effectifs: ",
                paste(table_counts, collapse = ", "), ")"),
         debug = debug)
  }

  # Initialiser r\u00e9sultats
  result <- list(
    ss_type = "III",  # D\u00e9faut conservateur
    reason = "",
    balanced = balanced,
    formula_optimal = TRUE,
    suggested_formula = NULL,
    warnings = character(0),
    k = k
  )

  # ==========================================================================
  # CAS 1: ANOVA \u00c0 DEUX FACTEURS (sans covariables)
  # ==========================================================================

  if (analysis_type == "ANOVA" && n_factors == 2 && n_covariates == 0) {

    .dbg("Analyzing ANOVA 2-way for SS type selection",
         "Analyse ANOVA 2-way pour s\u00e9lection type SS",
         debug = debug)

    if (balanced) {
      # Plan \u00e9quilibr\u00e9 \u2192 Type I par d\u00e9faut
      result$ss_type <- "I"
      result$reason <- .msg(
        "Balanced 2-way ANOVA design: Type I Sum of Squares (order-independent for balanced designs)",
        "Plan ANOVA 2-way \u00e9quilibr\u00e9: Sommes des Carr\u00e9s Type I (ind\u00e9pendant de l'ordre pour plans \u00e9quilibr\u00e9s)"
      )

    } else {
      # Plan d\u00e9s\u00e9quilibr\u00e9 \u2192 tester interaction si mod\u00e8le disponible
      if (!is.null(model)) {
        # V\u00e9rifier si formule contient interaction
        formula_str <- deparse1(formula)
        has_interaction <- grepl("\\*|:", formula_str)

        if (has_interaction) {
          # Tester significativit\u00e9 interaction avec anova()
          anova_result <- anova(model)

          # Trouver ligne interaction (derni\u00e8re ligne avant Residuals)
          interaction_row <- nrow(anova_result) - 1

          if (interaction_row > 0) {
            interaction_pval <- anova_result[interaction_row, "Pr(>F)"]

            if (!is.na(interaction_pval) && interaction_pval >= alpha) {
              # Interaction NON significative \u2192 Type II
              result$ss_type <- "II"
              result$reason <- .msg(
                paste0("Unbalanced design with non-significant interaction (p=",
                       round(interaction_pval, 3), "): Type II Sum of Squares"),
                paste0("Plan d\u00e9s\u00e9quilibr\u00e9 avec interaction non-significative (p=",
                       round(interaction_pval, 3), "): Sommes des Carr\u00e9s Type II")
              )
            } else {
              # Interaction significative \u2192 Type III
              result$ss_type <- "III"
              result$reason <- .msg(
                paste0("Unbalanced design with significant interaction (p=",
                       round(interaction_pval, 3), "): Type III Sum of Squares"),
                paste0("Plan d\u00e9s\u00e9quilibr\u00e9 avec interaction significative (p=",
                       round(interaction_pval, 3), "): Sommes des Carr\u00e9s Type III")
              )
            }
          }
        } else {
          # Pas d'interaction dans formule \u2192 Type II
          result$ss_type <- "II"
          result$reason <- .msg(
            "Unbalanced design without interaction: Type II Sum of Squares",
            "Plan d\u00e9s\u00e9quilibr\u00e9 sans interaction: Sommes des Carr\u00e9s Type II"
          )
        }
      } else {
        # Pas de mod\u00e8le pour tester \u2192 Type III conservateur
        result$ss_type <- "III"
        result$reason <- .msg(
          "Unbalanced design (interaction not tested): Type III Sum of Squares (conservative)",
          "Plan d\u00e9s\u00e9quilibr\u00e9 (interaction non test\u00e9e): Sommes des Carr\u00e9s Type III (conservateur)"
        )
      }
    }
  }

  # ==========================================================================
  # CAS 2: ANCOVA (n facteurs + m covariables)
  # ==========================================================================

  if (analysis_type == "ANCOVA" && n_covariates > 0) {

    .dbg("Analyzing ANCOVA for SS type selection",
         "Analyse ANCOVA pour s\u00e9lection type SS",
         debug = debug)

    # TOUJOURS v\u00e9rifier l'ordre des covariables dans la formule
    # (ind\u00e9pendamment de l'\u00e9quilibre du plan)
    # Extraire ordre des termes dans la formule
    formula_terms <- attr(terms(formula), "term.labels")

    # Interaction facteur\u00d7covariable CONSERV\u00c9E (pentes h\u00e9t\u00e9rog\u00e8nes) => Type III
    # L'ordre des termes est alors sans objet, et les tests restent align\u00e9s sur les
    # post-hocs (moyennes ajust\u00e9es \u00e0 la moyenne de la covariable centr\u00e9e).
    has_factor_cov_interaction <- any(vapply(
      strsplit(formula_terms[grepl(":", formula_terms, fixed = TRUE)], ":", fixed = TRUE),
      function(parts) any(parts %in% factor_vars) && any(parts %in% numeric_vars),
      logical(1)
    ))

    if (has_factor_cov_interaction) {
      result$ss_type <- "III"
      result$reason <- .msg(
        "Factor\u00d7covariate interaction kept in model: Type III Sum of Squares",
        "Interaction facteur\u00d7covariable conserv\u00e9e dans le mod\u00e8le : Sommes des Carr\u00e9s Type III"
      )
      return(result)
    }

    # Identifier position des covariables et facteurs
    cov_positions <- which(formula_terms %in% numeric_vars)
    factor_positions <- which(formula_terms %in% factor_vars)

    # V\u00e9rifier que covariables viennent AVANT facteurs
    covariates_first <- length(factor_positions) == 0 ||
                        length(cov_positions) == 0 ||
                        all(cov_positions < min(factor_positions))

    if (!covariates_first) {
      result$formula_optimal <- FALSE
      result$warnings <- c(result$warnings,
        .msg("Type III does not depend on the order of terms: this result stands as is.",
             "Le Type III ne d\u00e9pend pas de l'ordre des termes : ce r\u00e9sultat reste valable tel quel.")
      )

      # Construire formule sugg\u00e9r\u00e9e (covariables puis facteurs)
      suggested_terms <- c(numeric_vars, factor_vars)
      result$suggested_formula <- as.formula(paste0(response_var, " ~ ",
                                                     paste(suggested_terms, collapse = " + ")))

      # Basculer vers Type III
      result$ss_type <- "III"
      result$reason <- .msg(
        "Covariates not placed first in formula: Type III Sum of Squares (order-dependent Type I not optimal)",
        "Covariables non plac\u00e9es en premier dans la formule: Sommes des Carr\u00e9s Type III (Type I d\u00e9pendant de l'ordre non optimal)"
      )
      # NE PAS return ici - continuer pour v\u00e9rifier \u00e9quilibre et ajouter raisons suppl\u00e9mentaires
    }

    # V\u00e9rifier ordre des covariables (plus corr\u00e9l\u00e9e \u00e0 Y en premier)
    if (n_covariates > 1 && result$formula_optimal) {
      # Calculer corr\u00e9lations absolues avec Y
      cors <- sapply(numeric_vars, function(cov) {
        abs(cor(data[[response_var]], data[[cov]], use = "complete.obs"))
      })

      # Ordre optimal (d\u00e9croissant)
      optimal_order <- names(sort(cors, decreasing = TRUE))

      # Ordre actuel dans formule
      current_order <- numeric_vars[order(cov_positions)]

      if (!identical(current_order, optimal_order)) {
        result$formula_optimal <- FALSE
        result$warnings <- c(result$warnings,
          .msg(
            paste0("Covariates not in optimal order (should be sorted by decreasing correlation with DV):\n",
                   "  Current: ", paste(current_order, collapse = ", "), "\n",
                   "  Optimal: ", paste(optimal_order, " (|r|=", round(cors[optimal_order], 3), ")",
                                       collapse = ", ", sep = "")),
            paste0("Covariables non dans l'ordre optimal (devraient \u00eatre tri\u00e9es par corr\u00e9lation d\u00e9croissante avec VD):\n",
                   "  Actuel: ", paste(current_order, collapse = ", "), "\n",
                   "  Optimal: ", paste(optimal_order, " (|r|=", round(cors[optimal_order], 3), ")",
                                       collapse = ", ", sep = ""))
          )
        )

        # Construire formule sugg\u00e9r\u00e9e
        suggested_terms <- c(optimal_order, factor_vars)
        result$suggested_formula <- as.formula(paste0(response_var, " ~ ",
                                                       paste(suggested_terms, collapse = " + ")))

        # Basculer vers Type III
        result$ss_type <- "III"
        result$reason <- .msg(
          "Covariates not in optimal order: Type III Sum of Squares (Type I depends on order)",
          "Covariables non dans l'ordre optimal: Sommes des Carr\u00e9s Type III (Type I d\u00e9pend de l'ordre)"
        )
        # NE PAS return ici - continuer pour v\u00e9rifier \u00e9quilibre
      }
    }

    # V\u00e9rifier si d\u00e9s\u00e9quilibre \u2192 Type III obligatoire
    if (!balanced) {
      result$ss_type <- "III"
      # Type III impos\u00e9 quel que soit l'ordre : ne pas sugg\u00e9rer une formule
      # qui n'aboutirait de toute fa\u00e7on pas \u00e0 un Type I
      result$warnings <- character(0)
      result$suggested_formula <- NULL
      # Ajouter raison d\u00e9s\u00e9quilibre sans \u00e9craser les warnings pr\u00e9c\u00e9dents
      if (nchar(result$reason) == 0) {
        result$reason <- .msg(
          "Unbalanced factorial design: Type III Sum of Squares (mandatory for ANCOVA with imbalance)",
          "Plan factoriel d\u00e9s\u00e9quilibr\u00e9: Sommes des Carr\u00e9s Type III (obligatoire pour ANCOVA avec d\u00e9s\u00e9quilibre)"
        )
      } else {
        # Ajouter info d\u00e9s\u00e9quilibre aux raisons existantes
        result$reason <- paste0(result$reason, " + ", .msg(
          "Unbalanced design",
          "Plan d\u00e9s\u00e9quilibr\u00e9"
        ))
      }
      return(result)
    }

    # Si formule non-optimale mais plan \u00e9quilibr\u00e9 \u2192 garder Type III avec warnings
    if (!result$formula_optimal && balanced) {
      # Les warnings sont d\u00e9j\u00e0 ajout\u00e9s, result$ss_type est d\u00e9j\u00e0 "III"
      return(result)
    }

    # Si tout est optimal \u2192 Type I
    if (result$formula_optimal && balanced) {
      result$ss_type <- "I"
      result$reason <- .msg(
        "Balanced design with optimal formula (covariates first, sorted by correlation): Type I Sum of Squares",
        "Plan \u00e9quilibr\u00e9 avec formule optimale (covariables en premier, tri\u00e9es par corr\u00e9lation): Sommes des Carr\u00e9s Type I"
      )
    }
  }

  return(result)
}

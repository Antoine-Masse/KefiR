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
#' @export
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
       paste0("Sélection Type SS: ", n_factors, " facteurs, ", n_covariates, " covariables"),
       debug = debug)

  # ==========================================================================
  # DÉTECTION ÉQUILIBRE DU PLAN FACTORIEL
  # ==========================================================================

  balanced <- FALSE
  if (n_factors > 0) {
    # Créer interaction complète des facteurs
    if (n_factors == 1) {
      g_cat <- data[[factor_vars[1]]]
    } else {
      g_cat <- interaction(data[, factor_vars, drop = FALSE], drop = TRUE)
    }

    # Compter par groupe
    table_counts <- table(g_cat)

    # Équilibré si toutes les cellules ont même taille
    balanced <- (length(unique(table_counts)) == 1)

    .dbg(paste0("Plan balanced: ", balanced, " (counts: ",
                paste(table_counts, collapse = ", "), ")"),
         paste0("Plan équilibré: ", balanced, " (effectifs: ",
                paste(table_counts, collapse = ", "), ")"),
         debug = debug)
  }

  # Initialiser résultats
  result <- list(
    ss_type = "III",  # Défaut conservateur
    reason = "",
    balanced = balanced,
    formula_optimal = TRUE,
    suggested_formula = NULL,
    warnings = character(0),
    k = k
  )

  # ==========================================================================
  # CAS 1: ANOVA À DEUX FACTEURS (sans covariables)
  # ==========================================================================

  if (analysis_type == "ANOVA" && n_factors == 2 && n_covariates == 0) {

    .dbg("Analyzing ANOVA 2-way for SS type selection",
         "Analyse ANOVA 2-way pour sélection type SS",
         debug = debug)

    if (balanced) {
      # Plan équilibré → Type I par défaut
      result$ss_type <- "I"
      result$reason <- .msg(
        "Balanced 2-way ANOVA design: Type I Sum of Squares (order-independent for balanced designs)",
        "Plan ANOVA 2-way équilibré: Sommes des Carrés Type I (indépendant de l'ordre pour plans équilibrés)"
      )

    } else {
      # Plan déséquilibré → tester interaction si modèle disponible
      if (!is.null(model)) {
        # Vérifier si formule contient interaction
        formula_str <- deparse(formula)
        has_interaction <- grepl("\\*|:", formula_str)

        if (has_interaction) {
          # Tester significativité interaction avec anova()
          anova_result <- anova(model)

          # Trouver ligne interaction (dernière ligne avant Residuals)
          interaction_row <- nrow(anova_result) - 1

          if (interaction_row > 0) {
            interaction_pval <- anova_result[interaction_row, "Pr(>F)"]

            if (!is.na(interaction_pval) && interaction_pval >= alpha) {
              # Interaction NON significative → Type II
              result$ss_type <- "II"
              result$reason <- .msg(
                paste0("Unbalanced design with non-significant interaction (p=",
                       round(interaction_pval, 3), "): Type II Sum of Squares"),
                paste0("Plan déséquilibré avec interaction non-significative (p=",
                       round(interaction_pval, 3), "): Sommes des Carrés Type II")
              )
            } else {
              # Interaction significative → Type III
              result$ss_type <- "III"
              result$reason <- .msg(
                paste0("Unbalanced design with significant interaction (p=",
                       round(interaction_pval, 3), "): Type III Sum of Squares"),
                paste0("Plan déséquilibré avec interaction significative (p=",
                       round(interaction_pval, 3), "): Sommes des Carrés Type III")
              )
            }
          }
        } else {
          # Pas d'interaction dans formule → Type II
          result$ss_type <- "II"
          result$reason <- .msg(
            "Unbalanced design without interaction: Type II Sum of Squares",
            "Plan déséquilibré sans interaction: Sommes des Carrés Type II"
          )
        }
      } else {
        # Pas de modèle pour tester → Type III conservateur
        result$ss_type <- "III"
        result$reason <- .msg(
          "Unbalanced design (interaction not tested): Type III Sum of Squares (conservative)",
          "Plan déséquilibré (interaction non testée): Sommes des Carrés Type III (conservateur)"
        )
      }
    }
  }

  # ==========================================================================
  # CAS 2: ANCOVA (n facteurs + m covariables)
  # ==========================================================================

  if (analysis_type == "ANCOVA" && n_covariates > 0) {

    .dbg("Analyzing ANCOVA for SS type selection",
         "Analyse ANCOVA pour sélection type SS",
         debug = debug)

    # TOUJOURS vérifier l'ordre des covariables dans la formule
    # (indépendamment de l'équilibre du plan)
    # Extraire ordre des termes dans la formule
    formula_terms <- attr(terms(formula), "term.labels")

    # Identifier position des covariables et facteurs
    cov_positions <- which(formula_terms %in% numeric_vars)
    factor_positions <- which(formula_terms %in% factor_vars)

    # Vérifier que covariables viennent AVANT facteurs
    covariates_first <- length(factor_positions) == 0 ||
                        length(cov_positions) == 0 ||
                        all(cov_positions < min(factor_positions))

    if (!covariates_first) {
      result$formula_optimal <- FALSE
      result$warnings <- c(result$warnings,
        .msg("Covariates should appear BEFORE factors in ANCOVA formula",
             "Les covariables doivent apparaître AVANT les facteurs dans la formule ANCOVA")
      )

      # Construire formule suggérée (covariables puis facteurs)
      suggested_terms <- c(numeric_vars, factor_vars)
      result$suggested_formula <- as.formula(paste0(response_var, " ~ ",
                                                     paste(suggested_terms, collapse = " + ")))

      # Basculer vers Type III
      result$ss_type <- "III"
      result$reason <- .msg(
        "Covariates not placed first in formula: Type III Sum of Squares (order-dependent Type I not optimal)",
        "Covariables non placées en premier dans la formule: Sommes des Carrés Type III (Type I dépendant de l'ordre non optimal)"
      )
      # NE PAS return ici - continuer pour vérifier équilibre et ajouter raisons supplémentaires
    }

    # Vérifier ordre des covariables (plus corrélée à Y en premier)
    if (n_covariates > 1 && result$formula_optimal) {
      # Calculer corrélations absolues avec Y
      cors <- sapply(numeric_vars, function(cov) {
        abs(cor(data[[response_var]], data[[cov]], use = "complete.obs"))
      })

      # Ordre optimal (décroissant)
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
            paste0("Covariables non dans l'ordre optimal (devraient être triées par corrélation décroissante avec VD):\n",
                   "  Actuel: ", paste(current_order, collapse = ", "), "\n",
                   "  Optimal: ", paste(optimal_order, " (|r|=", round(cors[optimal_order], 3), ")",
                                       collapse = ", ", sep = ""))
          )
        )

        # Construire formule suggérée
        suggested_terms <- c(optimal_order, factor_vars)
        result$suggested_formula <- as.formula(paste0(response_var, " ~ ",
                                                       paste(suggested_terms, collapse = " + ")))

        # Basculer vers Type III
        result$ss_type <- "III"
        result$reason <- .msg(
          "Covariates not in optimal order: Type III Sum of Squares (Type I depends on order)",
          "Covariables non dans l'ordre optimal: Sommes des Carrés Type III (Type I dépend de l'ordre)"
        )
        # NE PAS return ici - continuer pour vérifier équilibre
      }
    }

    # Vérifier si déséquilibre → Type III obligatoire
    if (!balanced) {
      result$ss_type <- "III"
      # Ajouter raison déséquilibre sans écraser les warnings précédents
      if (nchar(result$reason) == 0) {
        result$reason <- .msg(
          "Unbalanced factorial design: Type III Sum of Squares (mandatory for ANCOVA with imbalance)",
          "Plan factoriel déséquilibré: Sommes des Carrés Type III (obligatoire pour ANCOVA avec déséquilibre)"
        )
      } else {
        # Ajouter info déséquilibre aux raisons existantes
        result$reason <- paste0(result$reason, " + ", .msg(
          "Unbalanced design",
          "Plan déséquilibré"
        ))
      }
      return(result)
    }

    # Vérifier interactions facteur×covariable
    formula_str <- deparse(formula)
    has_factor_cov_interaction <- FALSE

    for (fac in factor_vars) {
      for (cov in numeric_vars) {
        if (grepl(paste0(fac, ":", cov, "|", cov, ":", fac, "|", fac, "\\*", cov, "|", cov, "\\*", fac),
                  formula_str)) {
          has_factor_cov_interaction <- TRUE
          break
        }
      }
      if (has_factor_cov_interaction) break
    }

    if (has_factor_cov_interaction) {
      result$warnings <- c(result$warnings,
        .msg("Factor×Covariate interaction detected: this tests slope homogeneity (ANCOVA assumption)",
             "Interaction Facteur×Covariable détectée: teste l'homogénéité des pentes (assomption ANCOVA)")
      )
      # Note: le pipeline normal gérera cela via test homogénéité pentes
    }

    # Si formule non-optimale mais plan équilibré → garder Type III avec warnings
    if (!result$formula_optimal && balanced) {
      # Les warnings sont déjà ajoutés, result$ss_type est déjà "III"
      return(result)
    }

    # Si tout est optimal → Type I
    if (result$formula_optimal && balanced) {
      result$ss_type <- "I"
      result$reason <- .msg(
        "Balanced design with optimal formula (covariates first, sorted by correlation): Type I Sum of Squares",
        "Plan équilibré avec formule optimale (covariables en premier, triées par corrélation): Sommes des Carrés Type I"
      )
    }
  }

  return(result)
}

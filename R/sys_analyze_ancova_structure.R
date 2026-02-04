#' Analyse de la structure ANCOVA (pentes communes vs hétérogènes)
#'
#' Détermine si un modèle ANCOVA doit être traité comme :
#' - ANCOVA classique (pentes parallèles, pas d'interaction facteur×covariable)
#' - Modèle à pentes hétérogènes (interaction facteur×covariable présente)
#'
#' Effectue un test hiérarchique (anova reduced vs full) pour décider objectivement.
#'
#' @param formula Formule du modèle (peut contenir interactions)
#' @param data Data.frame contenant les variables
#' @param alpha Seuil de significativité pour test interaction (défaut 0.05)
#' @param verbose Logique, afficher messages pédagogiques
#' @param debug Logique, afficher messages debug
#' @param k Compteur messages .vbse()
#'
#' @return Liste avec :
#'   \itemize{
#'     \item model_type: "ANCOVA" ou "Heterogeneous_Slopes"
#'     \item has_interaction: Booléen, interaction facteur×covariable détectée
#'     \item interaction_terms: Vecteur noms des interactions détectées
#'     \item factor_vars: Noms facteurs catégoriques
#'     \item covariate_vars: Noms covariables continues
#'     \item hierarchical_test: Liste résultats anova(reduced, full) si effectué
#'     \item decision: "parallel_slopes" ou "heterogeneous_slopes"
#'     \item k: Compteur messages mis à jour
#'   }
#'
#' @details
#' **STANDARD ACADÉMIQUE** (Maxwell et al., 2018, pp. 423-428) :
#'
#' Lorsqu'une interaction facteur×covariable est présente dans la formule,
#' un test hiérarchique doit être effectué :
#'
#' 1. Ajuster modèle COMPLET (avec interaction)
#' 2. Ajuster modèle RÉDUIT (sans interaction, pentes communes)
#' 3. Comparer via anova(reduced, full)
#' 4. Si p < alpha : Garder interaction (pentes hétérogènes)
#'    Si p >= alpha : Supprimer interaction (pentes communes, ANCOVA classique)
#'
#' **CONSÉQUENCES DU CHOIX** :
#'
#' - ANCOVA classique :
#'   * Assomption : Homogénéité des pentes (parallélisme)
#'   * Post-hocs : Comparaisons moyennes ajustées
#'   * Graphique : Lignes parallèles
#'
#' - Pentes hétérogènes :
#'   * Pas d'assomption homogénéité pentes
#'   * Post-hocs : Comparaisons slopes, simple effects
#'   * Graphique : Lignes NON parallèles (interaction plot)
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018).
#'   Designing Experiments and Analyzing Data (3rd ed.). Routledge.
#'   Chapter 9, pp. 423-428: Testing homogeneity of regression slopes.
#'
#' @keywords internal
#' @export
.analyze_ancova_structure <- function(formula,
                                      data,
                                      alpha = 0.05,
                                      verbose = FALSE,
                                      debug = FALSE,
                                      k = NULL,
                                      code = NULL) {

  # Initialisation
  if (is.null(k)) k <- 0

  .dbg("=== Start .analyze_ancova_structure() ===",
       "=== D\u00e9but .analyze_ancova_structure() ===", debug = debug)

  # Extraire termes de la formule
  model_terms <- terms(formula, data = data)
  term_labels <- attr(model_terms, "term.labels")

  .dbg(paste0("Term labels: ", paste(term_labels, collapse = ", ")),
       paste0("Termes : ", paste(term_labels, collapse = ", ")),
       debug = debug)

  # S\u00e9parer termes simples et interactions
  simple_terms <- term_labels[!grepl(":", term_labels)]
  interaction_terms <- term_labels[grepl(":", term_labels)]

  .dbg(paste0("Simple terms: ", paste(simple_terms, collapse = ", ")),
       paste0("Termes simples : ", paste(simple_terms, collapse = ", ")),
       debug = debug)

  .dbg(paste0("Interactions: ", paste(interaction_terms, collapse = ", ")),
       paste0("Interactions : ", paste(interaction_terms, collapse = ", ")),
       debug = debug)

  # Identifier facteurs et covariables parmi les termes simples
  factor_vars <- character(0)
  covariate_vars <- character(0)

  for (term in simple_terms) {
    if (term %in% names(data)) {
      # Conversion character \u2192 factor comme .posthoc_ANCOVA
      if (is.factor(data[[term]]) || is.character(data[[term]])) {
        factor_vars <- c(factor_vars, term)
      } else if (is.numeric(data[[term]])) {
        covariate_vars <- c(covariate_vars, term)
      }
    }
  }

  .dbg(paste0("Factors: ", paste(factor_vars, collapse = ", ")),
       paste0("Facteurs : ", paste(factor_vars, collapse = ", ")),
       debug = debug)

  .dbg(paste0("Covariates: ", paste(covariate_vars, collapse = ", ")),
       paste0("Covariables : ", paste(covariate_vars, collapse = ", ")),
       debug = debug)

  # D\u00e9tecter interactions facteur\u00d7covariable
  factor_covariate_interactions <- character(0)

  for (interaction in interaction_terms) {
    parts <- strsplit(interaction, ":")[[1]]

    # V\u00e9rifier si interaction contient AU MOINS 1 facteur ET 1 covariable
    has_factor <- any(parts %in% factor_vars)
    has_covariate <- any(parts %in% covariate_vars)

    if (has_factor && has_covariate) {
      factor_covariate_interactions <- c(factor_covariate_interactions, interaction)
    }
  }

  has_interaction <- length(factor_covariate_interactions) > 0

  .dbg(paste0("Factor\u00d7Covariate interactions: ",
              ifelse(has_interaction, paste(factor_covariate_interactions, collapse = ", "), "NONE")),
       paste0("Interactions facteur\u00d7covariable : ",
              ifelse(has_interaction, paste(factor_covariate_interactions, collapse = ", "), "AUCUNE")),
       debug = debug)

  # ===========================================================================
  # CAS 1 : AUCUNE interaction facteur\u00d7covariable => ANCOVA classique directe
  # ===========================================================================

  if (!has_interaction) {
    k <- .vbse(
      paste0("Model structure: Classical ANCOVA (no factor\u00d7covariate interaction)\n",
             "\tParallel slopes assumed (homogeneity of regression slopes)"),
      paste0("Structure du mod\u00e8le : ANCOVA classique (pas d'interaction facteur\u00d7covariable)\n",
             "\tPentes parall\u00e8les assum\u00e9es (homog\u00e9n\u00e9it\u00e9 des pentes de r\u00e9gression)"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    result <- list(
      model_type = "ANCOVA",
      has_interaction = FALSE,
      interaction_terms = character(0),
      factor_vars = factor_vars,
      covariate_vars = covariate_vars,
      hierarchical_test = NULL,
      decision = "parallel_slopes",
      k = k
    )

    .dbg("=== End .analyze_ancova_structure() - Type: ANCOVA ===",
         "=== Fin .analyze_ancova_structure() - Type: ANCOVA ===",
         debug = debug)

    return(result)
  }

  # ===========================================================================
  # CAS 2 : Interaction(s) facteur\u00d7covariable d\u00e9tect\u00e9e(s)
  #         => TEST HI\u00c9RARCHIQUE pour d\u00e9cider
  # ===========================================================================

  # R\u00c9F\u00c9RENCE ACAD\u00c9MIQUE (d\u00e9veloppeurs/documentation uniquement - bp.log 7.4.6.1):
  # Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018).
  # Designing Experiments and Analyzing Data (3rd ed.). Routledge.
  # Chapter 9, pp. 423-428: Testing homogeneity of regression slopes assumption.

  k <- .vbse(
    paste0("DETECTED: Factor\u00d7covariate interaction in formula\n",
           "\tInteraction(s): ", paste(factor_covariate_interactions, collapse = ", "), "\n",
           "\tPerforming hierarchical test to decide between:\n",
           "\t  \u2022 Parallel slopes (classical ANCOVA)\n",
           "\t  \u2022 Heterogeneous slopes (interaction model)"),
    paste0("D\u00c9TECT\u00c9 : Interaction facteur\u00d7covariable dans la formule\n",
           "\tInteraction(s) : ", paste(factor_covariate_interactions, collapse = ", "), "\n",
           "\tTest hi\u00e9rarchique pour d\u00e9cider entre :\n",
           "\t  \u2022 Pentes parall\u00e8les (ANCOVA classique)\n",
           "\t  \u2022 Pentes h\u00e9t\u00e9rog\u00e8nes (mod\u00e8le avec interaction)"),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  # Construire formule R\u00c9DUITE (sans interactions facteur\u00d7covariable)
  reduced_terms <- setdiff(term_labels, factor_covariate_interactions)
  response_var <- as.character(formula[[2]])
  formula_reduced <- as.formula(paste(response_var, "~", paste(reduced_terms, collapse = " + ")))

  # Construire formule COMPL\u00c8TE (avec interactions) - c'est la formule originale
  formula_full <- formula

  .dbg(paste0("Reduced formula: ", deparse(formula_reduced)),
       paste0("Formule r\u00e9duite : ", deparse(formula_reduced)),
       debug = debug)

  .dbg(paste0("Full formula: ", deparse(formula_full)),
       paste0("Formule compl\u00e8te : ", deparse(formula_full)),
       debug = debug)

  # Ajuster les deux mod\u00e8les
  mod_reduced <- aov(formula_reduced, data = data)
  mod_full <- aov(formula_full, data = data)

  # Test hi\u00e9rarchique
  hierarchical_test <- anova(mod_reduced, mod_full)

  # Extraire p-value de l'interaction
  p_interaction <- hierarchical_test[2, "Pr(>F)"]

  .dbg(paste0("Hierarchical test p-value: ", round(p_interaction, 5)),
       paste0("P-value test hi\u00e9rarchique : ", round(p_interaction, 5)),
       debug = debug)

  # D\u00e9cision
  if (p_interaction < alpha) {
    # Interaction SIGNIFICATIVE => Pentes h\u00e9t\u00e9rog\u00e8nes
    decision <- "heterogeneous_slopes"
    model_type <- "Heterogeneous_Slopes"

    k <- .vbse(
      paste0("HIERARCHICAL TEST RESULT:\n",
             "\tInteraction p-value = ", .format_pval(p_interaction), " (< ", alpha, ")\n",
             "\t==> SIGNIFICANT interaction detected\n\n",
             "MODEL SELECTED: Heterogeneous slopes model\n",
             "\tSlopes differ significantly across groups.\n",
             "\tHomogeneity of slopes assumption VIOLATED.\n",
             "\tInterpretation must focus on how the covariate-response\n",
             "\trelationship differs between groups (not main effects)."),
      paste0("R\u00c9SULTAT TEST HI\u00c9RARCHIQUE :\n",
             "\tP-value interaction = ", .format_pval(p_interaction), " (< ", alpha, ")\n",
             "\t==> Interaction SIGNIFICATIVE d\u00e9tect\u00e9e\n\n",
             "MOD\u00c8LE S\u00c9LECTIONN\u00c9 : Mod\u00e8le \u00e0 pentes h\u00e9t\u00e9rog\u00e8nes\n",
             "\tLes pentes diff\u00e8rent significativement selon les groupes.\n",
             "\tAssomption homog\u00e9n\u00e9it\u00e9 des pentes VIOL\u00c9E.\n",
             "\tL'interpr\u00e9tation doit se concentrer sur comment la relation\n",
             "\tcovariable-r\u00e9ponse diff\u00e8re entre groupes (pas effets principaux)."),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

  } else {
    # Interaction NON significative => Pentes communes (ANCOVA classique)
    decision <- "parallel_slopes"
    model_type <- "ANCOVA"

    k <- .vbse(
      paste0("HIERARCHICAL TEST RESULT:\n",
             "\tInteraction p-value = ", .format_pval(p_interaction), " (>= ", alpha, ")\n",
             "\t==> NON-significant interaction\n\n",
             "MODEL SELECTED: Classical ANCOVA (parallel slopes)\n",
             "\tSlopes do NOT differ significantly across groups.\n",
             "\tHomogeneity of slopes assumption MET.\n",
             "\tProceeding with classical ANCOVA (interaction removed)."),
      paste0("R\u00c9SULTAT TEST HI\u00c9RARCHIQUE :\n",
             "\tP-value interaction = ", .format_pval(p_interaction), " (>= ", alpha, ")\n",
             "\t==> Interaction NON significative\n\n",
             "MOD\u00c8LE S\u00c9LECTIONN\u00c9 : ANCOVA classique (pentes parall\u00e8les)\n",
             "\tLes pentes NE diff\u00e8rent PAS significativement selon groupes.\n",
             "\tAssomption homog\u00e9n\u00e9it\u00e9 des pentes RESPECT\u00c9E.\n",
             "\tPoursuite avec ANCOVA classique (interaction supprim\u00e9e)."),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
  }

  # Retour
  result <- list(
    model_type = model_type,
    has_interaction = TRUE,  # Interaction \u00e9tait dans formule (m\u00eame si non significative)
    interaction_terms = factor_covariate_interactions,
    factor_vars = factor_vars,
    covariate_vars = covariate_vars,
    hierarchical_test = list(
      test_result = hierarchical_test,
      p_value = p_interaction,
      alpha = alpha,
      formula_reduced = formula_reduced,
      formula_full = formula_full,
      mod_reduced = mod_reduced,
      mod_full = mod_full
    ),
    decision = decision,
    k = k
  )

  .dbg(paste0("=== End .analyze_ancova_structure() - Type: ", model_type, " ==="),
       paste0("=== Fin .analyze_ancova_structure() - Type: ", model_type, " ==="),
       debug = debug)

  return(result)
}

#' Analyse modèles mixtes (Linear Mixed-Effects Models)
#'
#' Fonction interne pour analyser données avec effets fixes et aléatoires.
#' Utilisée pour designs complexes: mesures répétées déséquilibrées, effets
#' clusters, structures imbriquées/croisées.
#'
#' @param x Vecteur numérique - Variable dépendante
#' @param g Facteur ou data.frame - Variable(s) groupante(s)
#' @param formula Formule (optionnel) - Peut inclure effets aléatoires
#' @param data Data.frame (optionnel) - Données source
#' @param paired Logique - Mesures répétées/appariées
#' @param id Caractère - Nom colonne identifiant sujet (si paired=TRUE)
#' @param within Vecteur caractères - Facteurs intra-sujets
#' @param between Vecteur caractères - Facteurs inter-sujets
#' @param random_slope Logique - Inclure pente aléatoire (défaut FALSE)
#' @param alpha Numérique - Seuil significativité (défaut 0.05)
#' @param k Entier - Compteur messages verbose
#' @param code Logique - Afficher code R reproductible
#' @param debug Logique - Messages débogage détaillés
#' @param verbose Logique - Messages explicatifs pédagogiques
#'
#' @details
#' Cette fonction implémente les recommandations académiques pour modèles mixtes:
#'
#' 1. DÉTECTION AUTOMATIQUE:
#'    - Effets aléatoires depuis formule ou paramètres (id, within, between)
#'    - Type structure: imbriquée vs croisée
#'    - Nécessité pente aléatoire (random slopes)
#'
#' 2. CONSTRUCTION FORMULE:
#'    - Intercept aléatoire: (1 | id)
#'    - Pente aléatoire: (predictor | id)
#'    - Effets croisés: (1 | factor1) + (1 | factor2)
#'
#' 3. VALIDATION MODÈLE:
#'    - Utilise valreg() pour diagnostics complets
#'    - Normalité résidus
#'    - Homoscédasticité
#'    - VIF (multicolinéarité)
#'    - Test rapport vraisemblance (LRT)
#'
#' 4. LIMITES AUTOMATISATION:
#'    - Messages clairs si structure trop complexe
#'    - Suggestions code manuel pour cas particuliers
#'    - Références académiques
#'
#' @return Liste avec éléments:
#'   - model: Objet lmerMod (modèle ajusté)
#'   - formula: Formule utilisée
#'   - validation: Résultats valreg()
#'   - anova_table: Tableau ANOVA type III
#'   - variance_components: Composantes variance (VarCorr)
#'   - random_effects: Effets aléatoires par groupe (ranef)
#'   - fixed_effects: Effets fixes (fixef)
#'   - global_pvalue: P-value globale pour return=FALSE
#'   - k: Compteur messages mis à jour
#'   - assumptions_checked: Liste diagnostics
#'   - warnings: Vecteur avertissements
#'
#' @note
#' Fonction interne - Non exportée
#' Utilisée par m.test() via .multi_factor_analysis()
#'
#' @references
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015).
#' Fitting linear mixed-effects models using lme4.
#' \emph{Journal of Statistical Software}, 67(1), 1-48.
#'
#' Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013).
#' Random effects structure for confirmatory hypothesis testing: Keep it maximal.
#' \emph{Journal of Memory and Language}, 68(3), 255-278.
#'
#' @keywords internal
#' @importFrom lmerTest lmer
#' @importFrom lme4 VarCorr ranef fixef
#' @importFrom stats formula terms
.mixed_model_analysis <- function(x = NULL,
                                  g = NULL,
                                  formula = NULL,
                                  data = NULL,
                                  paired = FALSE,
                                  id = NULL,
                                  within = NULL,
                                  between = NULL,
                                  random_slope = FALSE,
                                  alpha = 0.05,
                                  k = NULL,
                                  code = FALSE,
                                  debug = FALSE,
                                  verbose = FALSE) {

  # Helper function for code mode output
  .code_mixed <- function(step_num, title, code_lines) {
    cat(paste0("# ", step_num, ") ", title, "\n"))
    for (line in code_lines) {
      cat(paste0(line, "\n"))
    }
    cat("\n")
  }

  # Initialize code mode
  if (isTRUE(code)) {
    verbose_original <- verbose
    verbose <- FALSE
    k_code <- 0  # Compteur séparé pour mode code
  }

  # Initialisation k
  if (is.null(k)) k <- 0

  ################################################################################
  # AJUSTEMENT MODÈLE MIXTE
  ################################################################################
  # Références académiques (commentaires code):
  # - Barr et al. (2013). Random effects structure. DOI:10.1016/j.jml.2012.11.001
  # - Pinheiro & Bates (2000). Mixed-Effects Models. DOI:10.1007/978-1-4419-0318-1

  # Extraire info depuis formule si fournie
  has_formula_random <- FALSE
  if (!is.null(formula)) {
    formula_str <- deparse(formula)
    has_formula_random <- grepl("\\|", formula_str)
  }

  # Construire formule si non fournie
  lmer_formula <- NULL

  if (has_formula_random) {
    # Formule contient déjà effets aléatoires
    lmer_formula <- formula

  } else if (!is.null(id)) {
    # Construire formule depuis paramètres

    # Variables dépendante
    if (!is.null(formula)) {
      # Extraire VD depuis formule
      dv_name <- as.character(formula[[2]])
    } else if (!is.null(data) && is.character(x)) {
      dv_name <- x
    } else {
      dv_name <- "y"  # Fallback
    }

    # Effets fixes
    fixed_terms <- c()

    if (!is.null(within)) {
      fixed_terms <- c(fixed_terms, within)
    }

    if (!is.null(between)) {
      fixed_terms <- c(fixed_terms, between)
    }

    if (length(fixed_terms) == 0 && !is.null(formula)) {
      # Extraire depuis formule
      formula_terms <- attr(terms(formula), "term.labels")
      fixed_terms <- formula_terms[!grepl("\\|", formula_terms)]
    }

    # Construire partie fixe
    if (length(fixed_terms) > 0) {
      fixed_part <- paste(fixed_terms, collapse = " + ")
    } else {
      fixed_part <- "1"  # Modèle intercept-only
    }

    # Construire partie aléatoire
    random_model_type <- ""
    if (random_slope && length(within) > 0) {
      # Pente aléatoire sur premier facteur within
      random_part <- paste0("(", within[1], " | ", id, ")")
      random_model_type <- "random slopes"
    } else {
      # Intercept aléatoire seulement
      random_part <- paste0("(1 | ", id, ")")
      random_model_type <- "random intercept"
    }

    # Formule complète
    lmer_formula <- as.formula(paste0(dv_name, " ~ ", fixed_part, " + ", random_part))

  } else {
    # Pas d'effets aléatoires identifiés
    k <- .vbse(
      paste0("ERROR: Cannot build mixed model without random effects.\n",
             "       Please specify 'id' parameter or use formula with random effects (e.g., y ~ x + (1|id))"),
      paste0("ERREUR : Impossible construire modèle mixte sans effets aléatoires.\n",
             "         Veuillez spécifier paramètre 'id' ou utiliser formule avec effets aléatoires (ex: y ~ x + (1|id))"),
      verbose = verbose, code = code, k = k, cpt = "off"
    )

    return(list(
      model = NULL,
      formula = NULL,
      validation = FALSE,
      k = k,
      warnings = "No random effects specified - cannot fit mixed model"
    ))
  }

  # MESSAGE UNIQUE: Ajustement modèle avec toutes les infos
  k <- .vbse(
    paste0("Mixed-effects model fitting [lmer() {lmerTest}]\n",
           "\tFormula: ", deparse(lmer_formula, width.cutoff = 500), "\n",
           "\tRandom effects: ", if(exists("random_model_type")) random_model_type else "from formula", "\n",
           "\tMethod: REML (Restricted Maximum Likelihood)\n",
           "\tDF approximation: Satterthwaite (for p-values)"),
    paste0("Ajustement modèle à effets mixtes [lmer() {lmerTest}]\n",
           "\tFormule : ", deparse(lmer_formula, width.cutoff = 500), "\n",
           "\tEffets aléatoires : ", if(exists("random_model_type")) ifelse(random_model_type == "random intercept", "intercept aléatoire", "pentes aléatoires") else "depuis formule", "\n",
           "\tMéthode : REML (Maximum vraisemblance restreinte)\n",
           "\tApproximation DL : Satterthwaite (pour p-values)"),
    verbose = verbose, code = code, k = k, cpt = "on"
  )
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_mixed(k_code, "Ajustement modèle mixte (effets aléatoires)", c(
      "library(lmerTest)",
      paste0("lmer_formula <- ", deparse(lmer_formula)),
      "mixed_model <- lmer(lmer_formula, data = data, REML = TRUE)",
      "summary(mixed_model)",
      "anova(mixed_model)"
    ))
  }

  # Ajuster modèle
  mixed_model <- NULL
  fit_error <- NULL

  mixed_model <- tryCatch({
    # Supprimer warnings techniques (boundary singular fit, etc.)
    suppressWarnings(suppressMessages(
      lmerTest::lmer(lmer_formula, data = data, REML = TRUE)
    ))
  }, error = function(e) {
    fit_error <<- e$message
    return(NULL)
  })

  if (is.null(mixed_model)) {
    k <- .vbse(
      paste0("ERROR: Model fitting failed.\n",
             "       ", fit_error, "\n",
             "       Possible causes:\n",
             "       - Formula syntax error\n",
             "       - Insufficient data for random effects\n",
             "       - Convergence issues (try simpler random structure)"),
      paste0("ERREUR : échec ajustement modèle.\n",
             "         ", fit_error, "\n",
             "         Causes possibles :\n",
             "         - Erreur syntaxe formule\n",
             "         - Données insuffisantes pour effets aléatoires\n",
             "         - Problèmes convergence (essayer structure aléatoire plus simple)"),
      verbose = verbose, code = code, k = k, cpt = "off"
    )

    return(list(
      model = NULL,
      formula = lmer_formula,
      validation = FALSE,
      k = k,
      warnings = fit_error
    ))
  }

#return(mixed_model)
  ################################################################################
  # VALIDATION MODÈLE VIA valreg()
  ################################################################################

  validation_result <- FALSE
  validation_warnings <- c()

	.dbg(NULL,"On lance valreg()",debug=debug)
	# NOTE: Ne PAS retourner directement le modèle ici - continuer pour valreg() et construction bilan
	# return(mixed_model)  # LIGNE DÉSACTIVÉE - causait erreur "bilan[[4]] : indice hors limites"
  # Appel valreg() - son bilan s'affiche directement avec numérotation continue
  valreg_output <- tryCatch({
    valreg(
      reg = mixed_model,
      verbose = verbose, code = code,  # CRITIQUE: doit être TRUE pour afficher
      alpha = alpha,
      boot = FALSE,
      plot = FALSE,
      k = k,
      tolerance = "extrem",
      orderDW = NULL,
	  debug=debug
    )
  }, error = function(e) {
    k <<- .vbse(
      paste0("ERROR: valreg() validation failed: ", e$message),
      paste0("ERREUR : Validation valreg() échouée : ", e$message),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    validation_warnings <<- c(validation_warnings, paste0("valreg() error: ", e$message))
    return(list(valid = FALSE, k = k))
  })
	.dbg(NULL,"Fin de valreg()",debug=debug)
  # Extract results from valreg() new return format
  validation_result <- valreg_output$valid
  k <- valreg_output$k

  # Message selon validation
  if (!isTRUE(validation_result)) {
    # Référence: Pinheiro & Bates (2000) Mixed-Effects Models in S and S-PLUS
    # Référence: Bolker et al. (2009) Generalized linear mixed models: a practical guide
    k <- .vbse(
      paste0("WARNING: Some validation checks failed (see valreg() output above).\n",
             "         Consider if violations are severe for your analysis:\n",
             "         - Non-significant coefficients: May indicate weak effects or multicollinearity\n",
             "         - Normality violations: Consider data transformation (log, sqrt, Box-Cox) or GLMM\n",
             "         - Heteroscedasticity: Consider nlme::lme() with variance structure\n",
             "         - Model simplification: Reduce random effects structure if convergence issues"),
      paste0("ATTENTION : Certaines vérifications ont échoué (voir bilan valreg() ci-dessus).\n",
             "            Évaluer si les violations sont sévères pour votre analyse :\n",
             "            - Coefficients non-significatifs : Effets faibles possibles ou multicolinéarité\n",
             "            - Violations normalité : Transformation données (log, sqrt, Box-Cox) ou GLMM\n",
             "            - Hétéroscédasticité : Utiliser nlme::lme() avec structure de variance\n",
             "            - Simplification modèle : Réduire structure effets aléatoires si problèmes convergence"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
  }

  ################################################################################
  # BLOC 5: RÉSULTATS STATISTIQUES DU MODÈLE MIXTE
  ################################################################################

  # ANOVA Type III avec stats::anova() (méthode lmerTest pour approximation Satterthwaite)
  # Note: lmerTest provides anova() method via stats::anova() generic
  anova_table <- stats::anova(mixed_model, type = "III")

  # Extraire p-value globale (premier effet fixe non-intercept)
  global_pvalue <- NA
  if (nrow(anova_table) > 0) {
    pval_col <- which(colnames(anova_table) %in% c("Pr(>F)", "p.value", "P(>|F|)"))
    if (length(pval_col) > 0) {
      row_names <- rownames(anova_table)
      valid_rows <- which(!grepl("Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
      if (length(valid_rows) > 0) {
        global_pvalue <- anova_table[valid_rows[1], pval_col[1]]
      }
    }
  }

  # Identifier effets significatifs pour interprétation
  significant_effects <- c()
  if (length(pval_col) > 0) {
    sig_rows <- which(anova_table[, pval_col[1]] < alpha)
    if (length(sig_rows) > 0) {
      significant_effects <- rownames(anova_table)[sig_rows]
    }
  }

  if (verbose) {
    k <- .vbse(
      "Type III ANOVA for mixed model [anova() {lmerTest}]\n\tDF approximation: Satterthwaite",
      "ANOVA Type III pour modèle mixte [anova() {lmerTest}]\n\tApproximation DL : Satterthwaite",
      verbose = verbose, code = code, k = k, cpt = "off"
    )
    cat("\n")
    print(anova_table)
    cat("\n")

    # INTERPRÉTATION immédiate après le tableau ANOVA
    if (length(significant_effects) > 0) {
      k <- .vbse(
        paste0("\tInterpretation: Significant fixed effect(s): ", paste(significant_effects, collapse = ", "), "\n",
               "\t                These factors show population-level effects on the response.\n",
               "\t                Random effects capture subject-specific variations."),
        paste0("\tInterprétation : Effet(s) fixe(s) significatif(s) : ", paste(significant_effects, collapse = ", "), "\n",
               "\t                 Ces facteurs montrent des effets au niveau populationnel.\n",
               "\t                 Les effets aléatoires capturent les variations inter-sujets."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    } else {
      k <- .vbse(
        paste0("\tInterpretation: No significant fixed effects at alpha = ", alpha),
        paste0("\tInterprétation : Aucun effet fixe significatif à alpha = ", alpha),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
    cat("\n")
  }

  # Composantes variance avec lme4::VarCorr()
  var_comp_S4 <- lme4::VarCorr(mixed_model)
  # Convertir en dataframe pour éviter erreurs S4 as.vector()
  var_comp <- as.data.frame(var_comp_S4)

  if (verbose) {
    k <- .vbse(
      "Variance components [VarCorr() {lme4}]:",
      "Composantes variance [VarCorr() {lme4}] :",
      verbose = verbose, code = code, k = k, cpt = "off"
    )
    cat("\n")
    print(var_comp_S4)  # Afficher format original pour lisibilité
    cat("\n")
  }

  # Effets fixes avec lme4::fixef()
  fixed_eff <- lme4::fixef(mixed_model)

  if (verbose) {
    k <- .vbse(
      "Fixed effects coefficients [fixef() {lme4}]:",
      "Coefficients effets fixes [fixef() {lme4}] :",
      verbose = verbose, code = code, k = k, cpt = "off"
    )
    cat("\n")
    print(fixed_eff)
    cat("\n")
  }

  # Effets aléatoires par groupe avec lme4::ranef()
  random_eff <- lme4::ranef(mixed_model)

  ################################################################################
  # BLOC 7: MODE CODE
  ################################################################################

  # NOTE: Centralized code_str generation removed (bp.log 7.4.6).
  # Mixed model code generation will be implemented via coordinated messages in future.

  ################################################################################
  # BLOC 8: RETURN
  ################################################################################

  # Préparer infos pour post-hoc
  posthoc_applicable <- length(significant_effects) > 0

  result <- list(
    model = mixed_model,
    formula = lmer_formula,
    validation = validation_result,
    anova_table = anova_table,
    variance_components = var_comp,
    random_effects = random_eff,
    fixed_effects = fixed_eff,
    global_pvalue = global_pvalue,
    significant_effects = significant_effects,
    posthoc_applicable = posthoc_applicable,
    k = k,
    assumptions_checked = list(
      validation_passed = validation_result
    ),
    warnings = validation_warnings
  )

  return(result)
}

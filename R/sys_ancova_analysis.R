#' Analysis of Covariance (ANCOVA) - Complete Pipeline
#'
#' Performs ANCOVA with exhaustive assumption checks and automatic robust
#' methods if violations detected. Follows academic standards for covariate
#' analysis in factorial designs.
#'
#' @param x Numeric vector. Dependent variable.
#' @param g Data frame. Contains categorical factors AND continuous covariates.
#' @param formula Formula specifying the model (e.g., y ~ factor1 * factor2 + covariate).
#' @param data Data frame containing all variables.
#' @param alpha Numeric. Significance level (default = 0.05).
#' @param k Integer. Message counter for verbose output.
#' @param code Logical. If TRUE, displays R code for reproducibility.
#' @param debug Logical. If TRUE, displays debug messages.
#' @param verbose Logical. If TRUE, provides detailed output.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{model}: Fitted ANCOVA model (or NULL if robust)
#'     \item \code{assumptions_checked}: List of all assumption test results
#'     \item \code{robust}: Logical indicating if robust methods were used
#'     \item \code{robust_results}: Results from robust analysis if applicable
#'     \item \code{k}: Updated message counter
#'   }
#'
#' @details
#' **ANCOVA ASSUMPTIONS (Academic Standard)**
#'
#' Following Maxwell, Delaney & Kelley (2018, Chapter 9), ANCOVA requires:
#'
#' 1. **Independence of observations** (design-level, cannot be tested statistically)
#' 2. **Normality of residuals** (Shapiro-Wilk, Jarque-Bera)
#' 3. **Homoscedasticity** (Bartlett on factor interaction, NOT including covariates)
#' 4. **Linearity** (DV ~ covariate relationship must be linear)
#' 5. **Homogeneity of regression slopes** (factor:covariate interaction non-significant)
#'
#' **References:**
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018). Designing Experiments
#' and Analyzing Data: A Model Comparison Perspective (3rd ed.). Routledge.
#' https://doi.org/10.4324/9781315642956
#'
#' Huitema, B. E. (2011). The Analysis of Covariance and Alternatives:
#' Statistical Methods for Experiments, Quasi-Experiments, and Single-Case Studies
#' (2nd ed.). Wiley. https://doi.org/10.1002/9781118067475
#'
#' @export
# Fonction auxiliaire pour afficher code R avec numérotation
.code_ancova <- function(step_num, title, code_lines) {
  cat(paste0("# ", step_num, ") ", title, "\n"))
  for (line in code_lines) {
    cat(paste0(line, "\n"))
  }
  cat("\n")
}

.ancova_analysis <- function(
    x = NULL,
    g = NULL,
    formula = NULL,
    data = NULL,
    paired = FALSE,
    id = NULL,
    alpha = 0.05,
    k = NULL,
    code = FALSE,
    debug = FALSE,
    verbose = FALSE
) {

  if (is.null(k)) k <- 0

  # Si code==TRUE, désactiver verbose pour éviter mélange texte/code
  # et initialiser compteur séparé pour numérotation des étapes code
  if (isTRUE(code)) {
    verbose_original <- verbose
    verbose <- FALSE
    k_code <- 0  # Compteur séparé pour mode code
  }

  .dbg("=== Start .ancova_analysis() ===",
       "=== Début de .ancova_analysis() ===", debug = debug)

  # ==========================================================================
  # DÉTECTION CAS LIMITES : PAIRED ANCOVA / RANDOM EFFECTS
  # ==========================================================================
  # Référence: Maxwell et al. (2018), Chapter 15, pp. 739-789
  # Paired ANCOVA requires mixed models (random intercepts for subjects)

  if (paired || !is.null(id)) {
    k <- .vbse(
      paste0("DETECTED: Paired/Repeated measures ANCOVA\n",
             "\tm.test() does not yet support paired ANCOVA or ANCOVA with random effects.\n\n",
             "RECOMMENDED APPROACH:\n",
             "\t1. Use mixed models with lmerTest package:\n",
             "\t   install.packages('lmerTest')\n",
             "\t   library(lmerTest)\n",
             "\t   model <- lmer(DV ~ factor + covariate + (1|subject_id), data=your_data)\n",
             "\t   anova(model)  # Type III tests with Satterthwaite DF\n",
             "\t2. Alternative: Repeated measures ANCOVA in SPSS (GLM Repeated Measures)\n",
             "\t3. Alternative: afex::aov_ez() for repeated measures designs\n\n",
             "WANT THIS FEATURE?\n",
             "\tContact: antoine.masse@u-bordeaux.fr\n",
             "\tPlease include:\n",
             "\t  • Your code line that triggered this message\n",
             "\t  • Your dataset (anonymized if confidential)\n",
             "\t  • Brief description of your research context"),
      paste0("DÉTECTÉ : ANCOVA appariée / Mesures répétées\n",
             "\tm.test() ne supporte pas encore l'ANCOVA appariée ou avec effets aléatoires.\n\n",
             "APPROCHE RECOMMANDÉE :\n",
             "\t1. Utiliser modèles mixtes avec package lmerTest :\n",
             "\t   install.packages('lmerTest')\n",
             "\t   library(lmerTest)\n",
             "\t   modele <- lmer(VD ~ facteur + covariable + (1|sujet_id), data=vos_donnees)\n",
             "\t   anova(modele)  # Tests Type III avec degrés liberté Satterthwaite\n",
             "\t2. Alternative : ANCOVA mesures répétées dans SPSS (MLG Mesures répétées)\n",
             "\t3. Alternative : afex::aov_ez() pour plans mesures répétées\n\n",
             "VOUS VOULEZ CETTE FONCTIONNALITÉ ?\n",
             "\tContact : antoine.masse@u-bordeaux.fr\n",
             "\tVeuillez inclure :\n",
             "\t  • Votre ligne de code ayant déclenché ce message\n",
             "\t  • Votre jeu de données (anonymisé si confidentiel)\n",
             "\t  • Brève description de votre contexte de recherche"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    .exit("Paired ANCOVA not supported. Use mixed models (lmerTest).",
          "ANCOVA appariée non supportée. Utiliser modèles mixtes (lmerTest).")
  }

  # ==========================================================================
  # DÉTECTION CAS LIMITES : MANCOVA (Multiple DVs)
  # ==========================================================================
  # Référence: Tabachnick & Fidell (2013), Chapter 7
  # MANCOVA = Multivariate ANCOVA (multiple dependent variables + covariates)

  # Check if formula left side has multiple DVs (cbind())
  # NOTE: as.character(formula[[2]]) returns vector for I(), log(), etc.
  # Example: I(log(A)) -> c("I", "log(A)") -> length = 2
  # Only cbind() indicates multiple DVs, not I() or other transformations
  formula_lhs <- deparse(formula[[2]])  # deparse gives single string

  if (grepl("^cbind\\s*\\(", formula_lhs, ignore.case = TRUE)) {
    k <- .vbse(
      paste0("DETECTED: MANCOVA (Multivariate ANCOVA)\n",
             "\tm.test() does not yet support MANCOVA (multiple dependent variables + covariates).\n\n",
             "RECOMMENDED APPROACH:\n",
             "\t1. Use stats::manova() with covariates:\n",
             "\t   # Example:\n",
             "\t   model <- lm(cbind(DV1, DV2, DV3) ~ factor + covariate, data=your_data)\n",
             "\t   manova_result <- manova(model)\n",
             "\t   summary(manova_result)  # Multivariate tests (Wilks, Pillai, etc.)\n",
             "\t2. Check assumptions: Box's M test (homogeneity of covariance matrices)\n",
             "\t3. For post-hocs: Use discriminant analysis or separate ANCOVAs with correction\n\n",
             "WANT THIS FEATURE?\n",
             "\tContact: antoine.masse@u-bordeaux.fr\n",
             "\tPlease include:\n",
             "\t  • Your code line that triggered this message\n",
             "\t  • Your dataset (anonymized if confidential)\n",
             "\t  • Brief description of your research context"),
      paste0("DÉTECTÉ : MANCOVA (ANCOVA multivariée)\n",
             "\tm.test() ne supporte pas encore MANCOVA (plusieurs variables dépendantes + covariables).\n\n",
             "APPROCHE RECOMMANDÉE :\n",
             "\t1. Utiliser stats::manova() avec covariables :\n",
             "\t   # Exemple :\n",
             "\t   modele <- lm(cbind(VD1, VD2, VD3) ~ facteur + covariable, data=vos_donnees)\n",
             "\t   resultat_manova <- manova(modele)\n",
             "\t   summary(resultat_manova)  # Tests multivariés (Wilks, Pillai, etc.)\n",
             "\t2. Vérifier assomptions : Test M de Box (homogénéité matrices covariances)\n",
             "\t3. Pour post-hocs : Analyse discriminante ou ANCOVAs séparées avec correction\n\n",
             "VOUS VOULEZ CETTE FONCTIONNALITÉ ?\n",
             "\tContact : antoine.masse@u-bordeaux.fr\n",
             "\tVeuillez inclure :\n",
             "\t  • Votre ligne de code ayant déclenché ce message\n",
             "\t  • Votre jeu de données (anonymisé si confidentiel)\n",
             "\t  • Brève description de votre contexte de recherche"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    .exit("MANCOVA not supported. Use stats::manova().",
          "MANCOVA non supportée. Utiliser stats::manova().")
  }

  # ==========================================================================
  # PHASE 1: ANALYSE STRUCTURE MODÈLE (ANCOVA vs PENTES HÉTÉROGÈNES)
  # ==========================================================================
  # Appel nouvelle fonction qui :
  # 1. Détecte interactions facteur×covariable
  # 2. Effectue test hiérarchique si interaction présente
  # 3. Décide ANCOVA classique vs pentes hétérogènes

  ancova_structure <- .analyze_ancova_structure(
    formula = formula,
    data = data,
    alpha = alpha,
    verbose = verbose, code = code,
    debug = debug,
    k = k
  )

  k <- ancova_structure$k  # Mettre à jour compteur messages

  # Extraire infos pour utilisation ultérieure
  model_type <- ancova_structure$model_type  # "ANCOVA" ou "Heterogeneous_Slopes"
  factor_vars <- ancova_structure$factor_vars
  numeric_vars <- ancova_structure$covariate_vars
  has_factor_cov_interaction <- ancova_structure$has_interaction
  hierarchical_test_result <- ancova_structure$hierarchical_test

  # ==========================================================================
  # VÉRIFICATION PRÉCOCE : TAILLE DES GROUPES
  # ==========================================================================
  # Vérifier qu'il y a au moins 2 observations par groupe pour éviter erreurs
  # dans les tests statistiques ultérieurs (Bartlett, Levene, etc.)

  if (length(factor_vars) > 0) {
    # Créer interaction de tous les facteurs pour vérification
    if (length(factor_vars) > 1) {
      g_check <- interaction(data[, factor_vars, drop = FALSE], drop = TRUE)
    } else {
      g_check <- factor(data[[factor_vars[1]]])
    }

    group_counts <- table(g_check)
    groups_with_lt2 <- names(group_counts)[group_counts < 2]

    if (length(groups_with_lt2) > 0) {
      k <- .vbse(
        paste0("WARNING: Insufficient observations in some groups.\n",
               "\tGroups with < 2 observations: ", paste(groups_with_lt2, collapse = ", "), "\n",
               "\tGroup sizes: ", paste(paste0(names(group_counts), "=", group_counts), collapse = ", "), "\n",
               "\tSome statistical tests (Bartlett, Levene) require at least 2 observations per group.\n",
               "\tAnalysis will continue but some tests may be skipped or may fail."),
        paste0("ATTENTION : Observations insuffisantes dans certains groupes.\n",
               "\tGroupes avec < 2 observations : ", paste(groups_with_lt2, collapse = ", "), "\n",
               "\tEffectifs par groupe : ", paste(paste0(names(group_counts), "=", group_counts), collapse = ", "), "\n",
               "\tCertains tests statistiques (Bartlett, Levene) nécessitent au moins 2 observations par groupe.\n",
               "\tL'analyse continuera mais certains tests pourraient être ignorés ou échouer."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }
  }

  # Si pentes hétérogènes, utiliser formule complète (avec interaction)
  # Si ANCOVA classique + interaction était présente, utiliser formule réduite
  if (model_type == "Heterogeneous_Slopes") {
    # Garder formule originale (avec interaction)
    formula_to_use <- formula
  } else if (model_type == "ANCOVA" && !is.null(hierarchical_test_result)) {
    # Interaction était présente mais non significative -> utiliser formule réduite
    formula_to_use <- hierarchical_test_result$formula_reduced
  } else {
    # Pas d'interaction dès le départ
    formula_to_use <- formula
  }

  .dbg(paste0("Model type determined: ", model_type),
       paste0("Type de modèle déterminé : ", model_type),
       debug = debug)

  .dbg(paste0("Formula to use: ", deparse(formula_to_use)),
       paste0("Formule à utiliser : ", deparse(formula_to_use)),
       debug = debug)

  # Message résumé structure (déjà affiché par .analyze_ancova_structure)
  # Déterminer le nombre d'assomptions selon le type de modèle
  n_assumptions <- if(has_factor_cov_interaction) 4 else 5

  k <- .vbse(
    paste0("\tFactors: ", paste(factor_vars, collapse = ", "), "\n",
           "\t\tCovariates: ", paste(numeric_vars, collapse = ", "), "\n",
           if(has_factor_cov_interaction) {
             paste0("\tInteraction in model: Heterogeneous slopes allowed\n",
                    "\t==> Complete assumption checking (", n_assumptions, " assumptions: NO slopes homogeneity test).")
           } else {
             paste0("\t==> Complete assumption checking (", n_assumptions, " assumptions).")
           }),
    paste0("\tFacteurs : ", paste(factor_vars, collapse = ", "), "\n",
           "\t\tCovariables : ", paste(numeric_vars, collapse = ", "), "\n",
           if(has_factor_cov_interaction) {
             paste0("\tInteraction dans modèle : Pentes hétérogènes autorisées\n",
                    "\t==> Vérifications complètes des assomptions (", n_assumptions, " assomptions : PAS de test homogénéité pentes).")
           } else {
             paste0("\t==> Vérifications complètes des assomptions (", n_assumptions, " assomptions).")
           }),
    verbose = verbose, code = code, k = k, cpt = "off"
  )

  # CODE: Étape 1 - Structure du modèle
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_ancova(k_code, "Structure du modèle ANCOVA", c())
  }

  # ==========================================================================
  # INITIALISATION DES VARIABLES D'ASSOMPTIONS (GLOBAL SCOPE)
  # ==========================================================================
  # Initialiser TOUTES les variables d'assomptions au début pour éviter
  # erreurs "objet introuvable" si early stop avant leur définition
  check_normality <- TRUE
  check_variance_equal <- TRUE
  check_linearity <- TRUE
  check_slopes_homogeneity <- TRUE
  linearity_results <- list()
  slopes_results <- list()
  goto_robust_ancova <- FALSE

  # Créer interaction des facteurs pour tests
  if (length(factor_vars) > 1) {
    g_cat <- interaction(data[, factor_vars, drop = FALSE], drop = TRUE)
  } else {
    g_cat <- data[[factor_vars[1]]]
  }

  # ==========================================================================
  # ASSOMPTION 1: INDÉPENDANCE DES OBSERVATIONS
  # ==========================================================================
  # Référence: Maxwell et al. (2018), Chapter 9, pp. 395-396

  k <- .vbse(
    paste0("ASSUMPTION 1/", n_assumptions, ": Independence of observations\n",
           "\tThis is a DESIGN-LEVEL assumption that cannot be tested statistically.\n",
           "\tVerify that:\n",
           "\t  • No repeated measures (each observation from different subject)\n",
           "\t  • No cluster effects (observations not grouped/nested)\n",
           "\t  • No carryover effects (order of measurements doesn't influence results)"),
    paste0("ASSOMPTION 1/", n_assumptions, " : Indépendance des observations\n",
           "\tC'est une assomption de PLAN qui ne peut pas être testée statistiquement.\n",
           "\tVérifiez que :\n",
           "\t  • Pas de mesures répétées (chaque observation d'un sujet différent)\n",
           "\t  • Pas d'effets cluster (observations non groupées/emboîtées)\n",
           "\t  • Pas d'effets report (ordre des mesures n'influence pas les résultats)"),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  # CODE: Étape 2 - Indépendance des observations
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_ancova(k_code, "ASSOMPTION 1 : Indépendance des observations", c())
  }

  assumptions_checked <- list(
    independence = list(
      tested = FALSE,
      note = "Design-level assumption - user verification required"
    )
  )

  # ==========================================================================
  # AJUSTEMENT DU MODÈLE (nécessaire pour assomptions suivantes)
  # ==========================================================================

  model <- tryCatch({
    # Construire liste de contrastes nommée pour Type III SS
    contrast_list <- list()
    for (fvar in factor_vars) {
      contrast_list[[fvar]] <- "contr.sum"
    }

    lm(formula, data = data, contrasts = contrast_list)
  }, error = function(e) {
    k <<- .vbse(
      paste0("Error fitting model: ", e$message),
      paste0("Erreur lors de l'ajustement du modèle : ", e$message),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    return(NULL)
  })

  if (is.null(model)) {
    return(list(
      model = NULL,
      assumptions_checked = assumptions_checked,
      robust = TRUE,
      g_cat = g_cat,
      k = k,
      error = "Model fitting failed"
    ))
  }

  # ==========================================================================
  # SÉLECTION INTELLIGENTE DU TYPE DE SOMMES DES CARRÉS (AVANT message ajustement)
  # ==========================================================================
  # Déterminer le type de SS optimal selon contexte (équilibre, ordre formule, etc.)
  ss_selection <- .select_ss_type(
    formula = formula,
    data = data,
    model = model,
    analysis_type = "ANCOVA",
    alpha = alpha,
    verbose = FALSE,  # Ne pas afficher messages ici, on les affiche après
    code = code,
    k = k,
    debug = debug
  )

  k <- ss_selection$k
  selected_type <- ss_selection$ss_type
  residuals_model <- residuals(model)

  # Références: Maxwell et al. (2018) Ch.9, DOI:10.4324/9781315642956
  #             Huitema (2011), DOI:10.1002/9781118067475

  # Construire texte méthode
  method_text_en <- if (selected_type == "I") {
    "anova() {stats}"
  } else {
    paste0("car::Anova(type='", selected_type, "')")
  }
  method_text_fr <- if (selected_type == "I") {
    "anova() {stats}"
  } else {
    paste0("car::Anova(type='", selected_type, "')")
  }

  # Déterminer si covariables mal placées (pour message restructuré)
  has_warnings <- length(ss_selection$warnings) > 0

  # MESSAGE RESTRUCTURÉ ÉTAPE 3 : Ajustement modèle ANCOVA
  # Structure : Formule → Type SS (si warning) → Contrastes → Modèle ajusté → Avertissements détaillés
  if (has_warnings && selected_type == "III") {
    # Cas covariables mal placées : Type III forcé
    k <- .vbse(
      paste0("Model fitting: ANCOVA with ", length(factor_vars), " factor(s) [lm()]\n",
             "\tFormula: ", deparse(formula), "\n",
             "\t==> Covariate(s) not placed first in formula:\n",
             "\t    ==> Current analysis will use Type III (robust to order)\n",
             "\t    --> Type III Sum of Squares (order-dependent Type I not optimal)\n",
             "\tContrasts: contr.sum (effects coding for unbiased estimates)\n",
             "\t    ==> Model fitted successfully (n=", nrow(data), ", residuals computed)\n",
             "\t        Type of Sums of Squares selected: Type ", selected_type, " [", method_text_en, "]"),
      paste0("Ajustement du modèle ANCOVA à ", length(factor_vars), " facteur(s) [lm()]\n",
             "\tFormule : ", deparse(formula), "\n",
             "\t==> Covariable(s) non placée(s) en premier dans la formule :\n",
             "\t    ==> Analyse actuelle utilisera Type III (robuste à l'ordre)\n",
             "\t    --> Sommes des Carrés Type III (Type I dépendant de l'ordre non optimal)\n",
             "\tContrastes : contr.sum (codage effets pour estimations non biaisées)\n",
             "\t    ==> Modèle ajusté avec succès (n=", nrow(data), ", résidus calculés)\n",
             "\t        Type de Sommes de Carrés sélectionné : Type ", selected_type, " [", method_text_fr, "]"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
  } else {
    # Cas normal (pas de warning ou Type I)
    k <- .vbse(
      paste0("Model fitting: ANCOVA with ", length(factor_vars), " factor(s) [lm()]\n",
             "\tFormula: ", deparse(formula), "\n",
             "\tContrasts: contr.sum (effects coding for unbiased estimates)\n",
             "\t==> Model fitted successfully (n=", nrow(data), ", residuals computed)\n",
             "\t    Type of Sums of Squares selected: Type ", selected_type, " [", method_text_en, "]"),
      paste0("Ajustement du modèle ANCOVA à ", length(factor_vars), " facteur(s) [lm()]\n",
             "\tFormule : ", deparse(formula), "\n",
             "\tContrastes : contr.sum (codage effets pour estimations non biaisées)\n",
             "\t==> Modèle ajusté avec succès (n=", nrow(data), ", résidus calculés)\n",
             "\t    Type de Sommes de Carrés sélectionné : Type ", selected_type, " [", method_text_fr, "]"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
  }

  # Afficher warnings détaillés si formule non-optimale
  if (has_warnings) {
    k <- .vbse(
      paste0("FORMULA OPTIMIZATION WARNINGS:\n",
             paste0("\t", ss_selection$warnings, collapse = "\n")),
      paste0("AVERTISSEMENTS OPTIMISATION FORMULE :\n",
             paste0("\t", ss_selection$warnings, collapse = "\n")),
      verbose = verbose, k = k, cpt = "off"
    )

    # Suggérer formule améliorée
    if (!is.null(ss_selection$suggested_formula)) {
      k <- .vbse(
        paste0("SUGGESTED FORMULA (for optimal Type I analysis):\n",
               "\t", deparse(ss_selection$suggested_formula), "\n",
               "--> Rerun m.test() with this formula for Type I Sum of Squares"),
        paste0("FORMULE SUGGÉRÉE (pour analyse Type I optimale) :\n",
               "\t", deparse(ss_selection$suggested_formula), "\n",
               "--> Relancer m.test() avec cette formule pour Sommes des Carrés Type I"),
        verbose = verbose, k = k, cpt = "off"
      )
    }
  }

  # CODE: Étape 3 - Ajustement modèle ANCOVA
  if (isTRUE(code)) {
    formula_str <- deparse(formula)
    contrast_code <- paste0("contrasts_list <- list(",
      paste(sapply(factor_vars, function(f) paste0(f, " = 'contr.sum'")), collapse = ", "),
      ")")

    k_code <- k_code + 1
    .code_ancova(k_code, "Ajustement du modèle ANCOVA", c(
      contrast_code,
      paste0("ancova_model <- lm(", formula_str, ", data = dt, contrasts = contrasts_list)"),
      "residuals_model <- residuals(ancova_model)",
      "",
      paste0("# Type de SS sélectionné : Type ", selected_type),
      if(selected_type == "I") "# Utiliser anova() pour Type I" else paste0("# Utiliser car::Anova(type='", selected_type, "') pour Type ", selected_type)
    ))
  }

  # ==========================================================================
  # DIAGNOSTIC: VIF (MULTICOLLINÉARITÉ)
  # ==========================================================================
  # Référence: Maxwell et al. (2018), pp. 57-59, 428
  # Checking multicollinearity only if ≥2 covariates

  vif_results <- NULL
  vif_warning <- FALSE

  if (length(numeric_vars) >= 2) {
    # Try loading car package for VIF
    car_available <- requireNamespace("car", quietly = TRUE)

    if (car_available) {
      k <- .vbse(
        paste0("DIAGNOSTIC: Multicollinearity detection [vif() {car}]\n",
               "\tVariance Inflation Factor measures correlation among covariates\n",
               "\tThresholds: VIF ≤ 5 (OK), 5 < VIF ≤ 10 (moderate), VIF > 10 (serious)"),
        paste0("DIAGNOSTIC : Détection multicolinéarité [vif() {car}]\n",
               "\tFacteur d'Inflation de Variance mesure corrélation entre covariables\n",
               "\tSeuils : VIF ≤ 5 (OK), 5 < VIF ≤ 10 (modéré), VIF > 10 (sérieux)"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      vif_values <- tryCatch({
        suppressWarnings(suppressMessages(car::vif(model)))
      }, error = function(e) {
        .dbg("VIF calculation failed", "Calcul VIF échoué", debug = debug)
        NULL
      })

      if (!is.null(vif_values)) {
        # Handle case where VIF returns matrix (for factors) or vector (for continuous)
        if (is.matrix(vif_values)) {
          # Extract GVIF^(1/(2*Df)) for comparability
          vif_values <- vif_values[, "GVIF^(1/(2*Df))"]
        }

        # Filter to only covariates (continuous predictors)
        vif_covariates <- vif_values[names(vif_values) %in% numeric_vars]

        vif_results <- list(
          values = vif_covariates,
          max_vif = max(vif_covariates),
          variables = names(vif_covariates)
        )

        # Check thresholds
        serious_vif <- vif_covariates > 10
        moderate_vif <- vif_covariates > 5 & vif_covariates <= 10

        if (any(serious_vif)) {
          vif_warning <- TRUE
          vars_serious <- names(vif_covariates[serious_vif])

          k <- .vbse(
            paste0("  WARNING: SERIOUS multicollinearity detected (VIF > 10)\n",
                   "\tAffected covariate(s): ", paste(vars_serious, collapse = ", "), "\n",
                   "\tVIF values: ", paste(sprintf("%s=%.2f", vars_serious, vif_covariates[serious_vif]), collapse = ", "), "\n",
                   "\tCoefficients may be UNSTABLE. Consider:\n",
                   "\t  • Removing one of the correlated covariates\n",
                   "\t  • Centering variables\n",
                   "\t  • Using Principal Component Analysis (PCA)"),
            paste0("  ATTENTION : Multicolinéarité SÉRIEUSE détectée (VIF > 10)\n",
                   "\tCovariable(s) affectée(s) : ", paste(vars_serious, collapse = ", "), "\n",
                   "\tValeurs VIF : ", paste(sprintf("%s=%.2f", vars_serious, vif_covariates[serious_vif]), collapse = ", "), "\n",
                   "\tLes coefficients peuvent être INSTABLES. Envisagez :\n",
                   "\t  • Supprimer une des covariables corrélées\n",
                   "\t  • Centrer les variables\n",
                   "\t  • Utiliser une Analyse en Composantes Principales (ACP)"),
            verbose = verbose, code = code, k = k, cpt = "on"
          )

        } else if (any(moderate_vif)) {
          vif_warning <- TRUE
          vars_moderate <- names(vif_covariates[moderate_vif])

          k <- .vbse(
            paste0("\tNOTE: Moderate multicollinearity detected (5 < VIF ≤ 10)\n",
                   "\t  Affected covariate(s): ", paste(vars_moderate, collapse = ", "), "\n",
                   "\t  VIF values: ", paste(sprintf("%s=%.2f", vars_moderate, vif_covariates[moderate_vif]), collapse = ", "), "\n",
                   "\t  Coefficients may be less stable. Monitor interpretation carefully."),
            paste0("\tNOTE : Multicolinéarité modérée détectée (5 < VIF ≤ 10)\n",
                   "\t  Covariable(s) affectée(s) : ", paste(vars_moderate, collapse = ", "), "\n",
                   "\t  Valeurs VIF : ", paste(sprintf("%s=%.2f", vars_moderate, vif_covariates[moderate_vif]), collapse = ", "), "\n",
                   "\t  Les coefficients peuvent être moins stables. Surveillez l'interprétation."),
            verbose = verbose, code = code, k = k, cpt = "off"
          )

        } else {
          k <- .vbse(
            paste0("\t==> All VIF values ≤ 5: No multicollinearity concern\n",
                   "\t    VIF values: ", paste(sprintf("%s=%.2f", names(vif_covariates), vif_covariates), collapse = ", ")),
            paste0("\t==> Toutes valeurs VIF ≤ 5 : Pas de préoccupation multicolinéarité\n",
                   "\t    Valeurs VIF : ", paste(sprintf("%s=%.2f", names(vif_covariates), vif_covariates), collapse = ", ")),
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      }

    } else {
      # car package not available
      k <- .vbse(
        paste0("DIAGNOSTIC: Multicollinearity check (VIF) - SKIPPED\n",
               "\tPackage 'car' not available. VIF cannot be calculated.\n",
               "\tTo enable: install.packages('car')"),
        paste0("DIAGNOSTIC : Contrôle multicolinéarité (VIF) - IGNORÉ\n",
               "\tPackage 'car' non disponible. VIF ne peut être calculé.\n",
               "\tPour activer : install.packages('car')"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }

    assumptions_checked$multicollinearity <- list(
      tested = car_available,
      vif_results = vif_results,
      warning = vif_warning,
      note = if(car_available) "VIF check performed on covariates" else "car package not available"
    )
  }

  # ==========================================================================
  # TEST INDÉPENDANCE COVARIABLES-FACTEURS (APPROCHE ROBUSTE)
  # ==========================================================================
  # Référence académique:
  #   Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018).
  #   Designing experiments and analyzing data: A model comparison perspective
  #   (3rd ed.). Routledge. https://doi.org/10.4324/9781315642956
  #   Chapter 9, pp. 401-405, 419-420: ANCOVA assumptions and covariate independence
  #
  # Rationale: ANCOVA assumes covariates are NOT affected by treatment (factor).
  # If violated, causal interpretation is compromised (covariate may be a mediator).
  #
  # Test strategy (inspired by .one_factor_analysis):
  #   1. Check normality of covariate by factor groups
  #   2. Check homogeneity of variances
  #   3. Select appropriate test:
  #      - ANOVA if assumptions met
  #      - Welch ANOVA if variances unequal
  #      - Kruskal-Wallis if normality violated

  covariate_independence_tests <- list()
  independence_violations <- FALSE
  independence_details <- c()

  for (cov in numeric_vars) {
    for (fac in factor_vars) {
      # Extract covariate data by factor levels
      cov_data <- data[[cov]]
      fac_data <- data[[fac]]

      # Remove NA values
      valid_idx <- !is.na(cov_data) & !is.na(fac_data)
      cov_clean <- cov_data[valid_idx]
      fac_clean <- fac_data[valid_idx]

      n_levels <- length(unique(fac_clean))
      test_method <- "aov()"

      .dbg(paste0("Testing independence: ", cov, " ~ ", fac, " (", n_levels, " levels)"),
           paste0("Test indépendance : ", cov, " ~ ", fac, " (", n_levels, " niveaux)"),
           debug = debug)

      # -----------------------------------------------------------------------
      # STEP 1: Check normality
      # -----------------------------------------------------------------------
      if (n_levels == 2) {
        # For 2 groups: test normality of each group
        normality_result <- .normality(cov_clean, fac_clean, alpha = alpha, k = 0, verbose = FALSE)
        is_normal <- normality_result[[1]]
      } else {
        # For >2 groups: test normality of residuals
        temp_aov <- aov(cov_clean ~ fac_clean)
        normality_result <- .normality(temp_aov$residuals, alpha = alpha, k = 0, verbose = FALSE)
        is_normal <- normality_result[[1]]
      }

      .dbg(paste0("  Normality: ", ifelse(is_normal, "YES", "NO")),
           paste0("  Normalité : ", ifelse(is_normal, "OUI", "NON")),
           debug = debug)

      # -----------------------------------------------------------------------
      # STEP 2: Check variance homogeneity (if normal)
      # -----------------------------------------------------------------------
      var_equal <- TRUE
      if (is_normal) {
        if (n_levels == 2) {
          # Fisher-Snedecor test for 2 groups
          var_test_result <- var.test(cov_clean ~ fac_clean)
          var_pval <- var_test_result$p.value
          var_equal <- var_pval >= alpha
        } else {
          # Bartlett test for >2 groups
          bartlett_result <- bartlett.test(cov_clean ~ fac_clean)
          var_pval <- bartlett_result$p.value
          var_equal <- var_pval >= alpha
        }

        .dbg(paste0("  Variance homogeneity: ", ifelse(var_equal, "YES", "NO"), " (p = ", round(var_pval, 4), ")"),
             paste0("  Homogénéité variances : ", ifelse(var_equal, "OUI", "NON"), " (p = ", round(var_pval, 4), ")"),
             debug = debug)
      }

      # -----------------------------------------------------------------------
      # STEP 3: Select and perform appropriate test
      # -----------------------------------------------------------------------
      test_result <- NULL
      p_value <- NA

      if (!is_normal) {
        # Non-normal → Kruskal-Wallis (non-parametric)
        test_method <- "kruskal.test()"
        test_result <- tryCatch({
          kruskal.test(cov_clean ~ fac_clean)
        }, error = function(e) NULL)

        if (!is.null(test_result)) {
          p_value <- test_result$p.value
        }

      } else if (!var_equal) {
        # Normal but heterogeneous variances → Welch ANOVA
        test_method <- "oneway.test()"
        test_result <- tryCatch({
          oneway.test(cov_clean ~ fac_clean, var.equal = FALSE)
        }, error = function(e) NULL)

        if (!is.null(test_result)) {
          p_value <- test_result$p.value
        }

      } else {
        # Normal and homogeneous → Classical ANOVA
        test_method <- "aov()"
        test_result <- tryCatch({
          aov(cov_clean ~ fac_clean)
        }, error = function(e) NULL)

        if (!is.null(test_result)) {
          test_summary <- summary(test_result)
          p_value <- test_summary[[1]]["Pr(>F)"][1, 1]
        }
      }

      # -----------------------------------------------------------------------
      # STEP 4: Interpret and store results
      # -----------------------------------------------------------------------
      if (!is.na(p_value)) {
        is_independent <- p_value >= alpha

        covariate_independence_tests[[paste0(cov, "_vs_", fac)]] <- list(
          covariate = cov,
          factor = fac,
          p_value = p_value,
          independent = is_independent,
          test_type = test_method,
          assumptions = list(
            normal = is_normal,
            var_equal = var_equal
          )
        )

        if (!is_independent) {
          independence_violations <- TRUE
          independence_details <- c(independence_details,
                                   paste0("'", cov, "' vs '", fac, "': p = ", .format_pval(p_value), " [", test_method, "]"))
        } else {
          .dbg(paste0("  Result: INDEPENDENT (p = ", round(p_value, 4), ") via ", test_method),
               paste0("  Résultat : INDÉPENDANTE (p = ", round(p_value, 4), ") via ", test_method),
               debug = debug)
        }
      } else {
        .dbg(paste0("Independence test failed for ", cov, " ~ ", fac),
             paste0("Test indépendance échoué pour ", cov, " ~ ", fac),
             debug = debug)
      }
    }
  }

  # Message consolidé: Test + Résultat + Décision
  # Référence (commentaire uniquement): Maxwell et al. (2018), pp. 401-405
  k <- .vbse(
    paste0("DIAGNOSTIC: Covariate-factor independence test [robust approach]\n",
           "\tTests if each covariate differs across factor levels\n",
           "\tMethod: For each covariate × factor combination:\n",
           "\t  1. Test normality of covariate by factor groups\n",
           "\t  2. Test variance homogeneity [var.test() or bartlett.test()]\n",
           "\t  3. Select appropriate test:\n",
           "\t     • Normal + homogeneous variances: ANOVA [aov() on 1 factor]\n",
           "\t     • Normal + unequal variances: Welch ANOVA [oneway.test()]\n",
           "\t     • Non-normal: Kruskal-Wallis [kruskal.test()]\n",
           if (independence_violations) {
             paste0("\t==> WARNING: Covariate depends on factor (violation detected)\n",
                    "\t    ", paste(independence_details, collapse = "\n\t    "), "\n",
                    "\t    Causal interpretation compromised: the factor may CAUSE changes in\n",
                    "\t    the covariate, making it a mediator rather than a true control variable.\n",
                    "\t    ANCOVA assumes covariates are pre-existing and independent of treatment.\n",
                    "\t==> Analysis continues with caution.")
           } else {
             paste0("\t==> All covariates independent from factor levels (all p >= ", alpha, ")\n",
                    "\t    Covariate independence assumption MET")
           }),
    paste0("DIAGNOSTIC : Test indépendance covariables-facteurs [approche robuste]\n",
           "\tTeste si chaque covariable diffère selon niveaux facteurs\n",
           "\tMéthode : Pour chaque combinaison covariable × facteur :\n",
           "\t  1. Tester normalité covariable par groupes de chaque facteur\n",
           "\t  2. Tester homogénéité variances [var.test() ou bartlett.test()]\n",
           "\t  3. Sélectionner test approprié :\n",
           "\t     • Normal + variances homogènes : ANOVA [aov() sur 1 facteur]\n",
           "\t     • Normal + variances inégales : Welch ANOVA [oneway.test()]\n",
           "\t     • Non-normal : Kruskal-Wallis [kruskal.test()]\n",
           if (independence_violations) {
             paste0("\t==> ATTENTION : Covariable dépend du facteur (violation détectée)\n",
                    "\t    ", paste(independence_details, collapse = "\n\t    "), "\n",
                    "\t    Interprétation causale compromise : le facteur peut CAUSER des changements\n",
                    "\t    dans la covariable, la rendant une variable médiatrice plutôt qu'un vrai\n",
                    "\t    contrôle.\n",
                    "\t    Rappel : l'ANCOVA suppose que les covariables sont préexistantes et\n",
                    "\t        indépendantes du traitement.\n",
                    "\t--> L'analyse continue avec précaution.")
           } else {
             paste0("\t==> Toutes covariables indépendantes des niveaux facteurs (tous p >= ", alpha, ")\n",
                    "\t    Assomption d'indépendance des covariables RESPECTÉE")
           }),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  assumptions_checked$covariate_independence <- list(
    tested = TRUE,
    tests = covariate_independence_tests,
    violations = independence_violations,
    note = "Robust testing: ANOVA/Welch/Kruskal-Wallis based on assumptions (Maxwell et al., 2018)"
  )

  # CODE: Étape 4 - Test indépendance covariables-facteurs
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_ancova(k_code, "Test indépendance covariables-facteurs", c(
      "# Appliquer KefiR::m.test() à chaque couple de covariable-facteur possible"
    ))
  }

  # ==========================================================================
  # DIAGNOSTIC: INFLUENCE DIAGNOSTICS (MODEL-BASED)
  # ==========================================================================
  # Références: Cook (1977), Fox (2016), Belsley et al. (1980)
  # NOTE: Outliers marginaux (IQR-based) NOT relevant in ANCOVA since we adjust
  #       for covariates. Only model influence diagnostics matter.

  influence_results <- NULL

  # -----------------------------------------------------------------------
  # DIAGNOSTIC INFLUENCE (Cook's D, Leverage, DFBETAS)
  # -----------------------------------------------------------------------

  influence_results <- .diagnostic_influence(
    model = model,
    data = data,
    alpha = alpha,
    verbose = FALSE,
    k = k,
    debug = debug
  )

  # -----------------------------------------------------------------------
  # MESSAGE: Influence Diagnostics
  # -----------------------------------------------------------------------

  # Préparer tableau synthétique des combinaisons de seuils franchis
  if (influence_results$n_influential > 0) {
    # Identifier quels seuils franchis pour chaque obs
    all_influential <- unique(influence_results$influential_any)
    has_cook <- all_influential %in% influence_results$influential_cook
    has_leverage <- all_influential %in% influence_results$influential_leverage
    has_dfbetas <- all_influential %in% influence_results$influential_dfbetas

    # Compter combinaisons
    n_cook_only <- sum(has_cook & !has_leverage & !has_dfbetas)
    n_leverage_only <- sum(!has_cook & has_leverage & !has_dfbetas)
    n_dfbetas_only <- sum(!has_cook & !has_leverage & has_dfbetas)
    n_cook_leverage <- sum(has_cook & has_leverage & !has_dfbetas)
    n_cook_dfbetas <- sum(has_cook & !has_leverage & has_dfbetas)
    n_leverage_dfbetas <- sum(!has_cook & has_leverage & has_dfbetas)
    n_all_three <- sum(has_cook & has_leverage & has_dfbetas)

    # Construire tableau synthétique
    summary_table <- c()
    if (n_cook_only > 0) summary_table <- c(summary_table, paste0("Cook only: ", n_cook_only))
    if (n_leverage_only > 0) summary_table <- c(summary_table, paste0("Leverage only: ", n_leverage_only))
    if (n_dfbetas_only > 0) summary_table <- c(summary_table, paste0("DFBETAS only: ", n_dfbetas_only))
    if (n_cook_leverage > 0) summary_table <- c(summary_table, paste0("Cook + Leverage: ", n_cook_leverage))
    if (n_cook_dfbetas > 0) summary_table <- c(summary_table, paste0("Cook + DFBETAS: ", n_cook_dfbetas))
    if (n_leverage_dfbetas > 0) summary_table <- c(summary_table, paste0("Leverage + DFBETAS: ", n_leverage_dfbetas))
    if (n_all_three > 0) summary_table <- c(summary_table, paste0("All three criteria: ", n_all_three))

    influence_summary_en <- paste0("\t==> ", influence_results$n_influential, " observation(s) exceed threshold(s):\n",
                                   "\t    ", paste(summary_table, collapse = ", "), "\n",
                                   if (influence_results$n_critical > 0) {
                                     paste0("\t    ⚠ WARNING: ", influence_results$n_critical, " observation(s) with Cook > 1 (critical)\n",
                                            "\t    Consider sensitivity analysis\n")
                                   } else {
                                     "\t    All Cook < 1 (no critical concerns)\n"
                                   },
                                   "\t==> Model ", if (influence_results$n_critical == 0) "robust" else "SENSITIVE", " to individual observations")

    summary_table_fr <- gsub("Cook only", "Cook seul", summary_table)
    summary_table_fr <- gsub("Leverage only", "Leverage seul", summary_table_fr)
    summary_table_fr <- gsub("DFBETAS only", "DFBETAS seul", summary_table_fr)
    summary_table_fr <- gsub("All three criteria", "Les trois critères", summary_table_fr)

    influence_summary_fr <- paste0("\t==> ", influence_results$n_influential, " observation(s) dépassent seuil(s) :\n",
                                   "\t    ", paste(summary_table_fr, collapse = ", "), "\n",
                                   if (influence_results$n_critical > 0) {
                                     paste0("\t    ⚠ ATTENTION : ", influence_results$n_critical, " observation(s) avec Cook > 1 (critique)\n",
                                            "\t    Considérer analyse de sensibilité\n")
                                   } else {
                                     "\t    Toutes Cook < 1 (aucune préoccupation critique)\n"
                                   },
                                   "\t==> Modèle ", if (influence_results$n_critical == 0) "robuste" else "SENSIBLE", " aux observations individuelles")
  } else {
    influence_summary_en <- "\t==> No influential observations\n\t==> Model robust to individual observations"
    influence_summary_fr <- "\t==> Aucune observation influente\n\t==> Modèle robuste aux observations individuelles"
  }

  k <- .vbse(
    paste0("DIAGNOSTIC: Model Influence [Cook's distance, Leverage, DFBETAS]\n",
           "\tDetects observations with strong impact on regression coefficients (β)\n",
           "\tThresholds: Cook = ", round(influence_results$thresholds$cook, 3), " (critical if > 1), Leverage = ", round(influence_results$thresholds$leverage, 3), "\n",
           influence_summary_en),
    paste0("DIAGNOSTIC : Influence sur le Modèle [Distance de Cook, Leverage, DFBETAS]\n",
           "\tDétecte observations ayant impact fort sur coefficients de régression (β)\n",
           "\tSeuils : Cook = ", round(influence_results$thresholds$cook, 3), " (critique si > 1), Leverage = ", round(influence_results$thresholds$leverage, 3), "\n",
           influence_summary_fr),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  assumptions_checked$influence <- influence_results

  # CODE: Étape 5 - Diagnostic d'influence
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_ancova(k_code, "Diagnostic influence sur le modèle", c(
      "# Distance de Cook (influence combinée)",
      "cook_d <- cooks.distance(ancova_model)",
      paste0("cook_threshold <- 4 / nrow(dt)  # = ", round(influence_results$thresholds$cook, 3)),
      "",
      "# Leverage (distance dans l'espace prédicteurs)",
      "leverage <- hatvalues(ancova_model)",
      paste0("leverage_threshold <- 2 * length(coef(ancova_model)) / nrow(dt)  # = ", round(influence_results$thresholds$leverage, 3)),
      "",
      "# DFBETAS (changement coefficients si observation retirée)",
      "dfbetas_vals <- dfbetas(ancova_model)",
      "dfbetas_threshold <- 2 / sqrt(nrow(dt))",
      "",
      "# Identifier observations influentes",
      "influential <- which(cook_d > cook_threshold | leverage > leverage_threshold |",
      "                     apply(abs(dfbetas_vals) > dfbetas_threshold, 1, any))"
    ))
  }

  # ==========================================================================
  # ASSOMPTION 2: NORMALITÉ DES RÉSIDUS
  # ==========================================================================

  # Déterminer méthode selon taille échantillon
  n_residuals <- length(residuals_model)
  normality_method <- if(n_residuals <= 50) "shapiro.test() {stats}" else if(n_residuals <= 500) "jb.norm.test() {KefiR}" else "skewness/kurtosis {KefiR}"

  k <- .vbse(
    paste0("ASSUMPTION 2/", n_assumptions, ": Normality of residuals (ACADEMIC control)\n",
           "\tTest selection based on sample size:\n",
           "\t  • n ≤ 50: Shapiro-Wilk [shapiro.test() {stats}]\n",
           "\t  • 50 < n ≤ 500: Jarque-Bera [jb.norm.test() {KefiR}]\n",
           "\t  • n > 500: Skewness/Kurtosis [{KefiR}]\n",
           "\tTest used: [", normality_method, "]"),
    paste0("ASSOMPTION 2/", n_assumptions, " : Normalité des résidus (contrôle ACADÉMIQUE)\n",
           "\tSélection du test selon taille échantillon :\n",
           "\t  • n ≤ 50 : Shapiro-Wilk [shapiro.test() {stats}]\n",
           "\t  • 50 < n ≤ 500 : Jarque-Bera [jb.norm.test() {KefiR}]\n",
           "\t  • n > 500 : Skewness/Kurtosis [{KefiR}]\n",
           "\tTest utilisé : [", normality_method, "]"),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  normality_result <- .normality(
    residuals_model,
    alpha = alpha,
    verbose = verbose,
    k = k,
    debug = debug,
    cpt = "off"
  )

  k <- normality_result[[2]]
  check_normality <- normality_result[[1]]

  assumptions_checked$normality <- list(
    passed = check_normality,
    method = if(n_residuals <= 50) "Shapiro-Wilk" else if(n_residuals <= 500) "Jarque-Bera" else "Skewness/Kurtosis"
  )

  # NOTE: Normalité peut bénéficier d'un contrôle tolérant (CLT si n>30)
  # Pas d'early stop ici - décision finale après vérification autres assumptions
  goto_robust_ancova <- FALSE

  # CODE: Étape 6 - Normalité des résidus
  if (isTRUE(code)) {
    test_used <- if(n_residuals <= 50) "shapiro.test" else if(n_residuals <= 500) "jb.norm.test" else "skewness/kurtosis"
    k_code <- k_code + 1
    .code_ancova(k_code, "Test de normalité des résidus", c(
      paste0("# Test utilisé: ", test_used),
      "residuals_model <- residuals(ancova_model)",
      if(n_residuals <= 50) "normality_test <- shapiro.test(residuals_model)" else "# library(KefiR); normality_test <- jb.norm.test(residuals_model)"
    ))
  }

  # ==========================================================================
  # ASSOMPTION 3: HOMOGÉNÉITÉ DES VARIANCES (HOMOSCÉDASTICITÉ)
  # ==========================================================================
  # Référence: Maxwell et al. (2018), pp. 401-402
  # Skip if already violated
  if (!goto_robust_ancova) {

  # Déterminer méthode selon normalité
  variance_method <- if(check_normality) "bartlett.test() {stats}" else "leveneTest() {car}"

  k <- .vbse(
    paste0("ASSUMPTION 3/", n_assumptions, ": Homoscedasticity (ACADEMIC control)\n",
           "\tIMPORTANT: Test on FACTOR INTERACTION ONLY (covariates excluded)\n",
           "\tRationale: Covariates are continuous => homogeneity tested on factor groups only\n",
           "\tTest used: [", variance_method, "] (selected based on normality result)\n",
           "\t  • Residuals normal: Bartlett [bartlett.test() {stats}]\n",
           "\t  • Residuals non-normal: Levene [leveneTest() {car}]"),
    paste0("ASSOMPTION 3/", n_assumptions, " : Homoscédasticité (contrôle ACADÉMIQUE)\n",
           "\tIMPORTANT : Test sur INTERACTION DES FACTEURS UNIQUEMENT (pas les covariables)\n",
           "\tJustification : Les covariables sont continues => homogénéité testée sur groupes facteurs uniquement\n",
           "\tTest utilisé : [", variance_method, "] (sélectionné selon résultat normalité)\n",
           "\t  • Résidus normaux : Bartlett [bartlett.test() {stats}]\n",
           "\t  • Résidus non-normaux : Levene [leveneTest() {car}]"),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  variance_result <- .variance(
    x = residuals_model,
    g = g_cat,
    check_normality = check_normality,
    alpha = alpha,
    paired = FALSE,
    verbose = FALSE,  # Ne pas afficher header redondant
    k = k,
    debug = debug,
    cpt = "off"
  )

  k <- variance_result[[2]]
  check_variance_equal <- variance_result[[1]]

  # Afficher juste le résultat (sans header redondant)
  # Si normalité violée MAIS homoscédasticité OK => message d'espoir pour re-test normalité
  if (check_variance_equal && !check_normality) {
    k <- .vbse(
      paste0("\t==> Homogeneous variances (p >= ", alpha, ")\n",
             "\t--> Homoscedasticity invites a more tolerant normality check."),
      paste0("\t==> Variances homogènes (p >= ", alpha, ")\n",
             "\t--> L'homoscédasticité invite à un contrôle plus tolérant de la normalité."),
      verbose = verbose, code = code, k = k, cpt = "off"
    )
  } else {
    k <- .vbse(
      if (check_variance_equal) {
        paste0("\t==> Homogeneous variances (p >= ", alpha, ").")
      } else {
        paste0("\t==> Heterogeneous variances (p < ", alpha, ").")
      },
      if (check_variance_equal) {
        paste0("\t==> Variances homogènes (p >= ", alpha, ").")
      } else {
        paste0("\t==> Variances hétérogènes (p < ", alpha, ").")
      },
      verbose = verbose, code = code, k = k, cpt = "off"
    )
  }

  assumptions_checked$homoscedasticity <- list(
    passed = check_variance_equal,
    method = if(check_normality) "Bartlett" else "Levene",
    note = "Test on factor interaction only (covariates excluded)"
  )

  # CODE: Étape 7 - Homoscédasticité
  if (isTRUE(code)) {
    k_code <- k_code + 1
    .code_ancova(k_code, "Test d'homoscédasticité", c(
      "# Test utilisé si résidus normaux: bartlett.test",
      if(check_normality) {
        paste0("variance_test <- bartlett.test(residuals_model ~ interaction(", paste(factor_vars, collapse = ", "), "))")
      } else {
        "# library(car); variance_test <- leveneTest(residuals_model ~ interaction(...))"
      }
    ))
  }

  # EARLY STOP si homoscédasticité violée
  # Référence : Maxwell et al. (2018), pp. 400-401
  if (!check_variance_equal) {
    k <- .vbse(
      paste0("  ==> CRITICAL VIOLATIONS detected: ",
             if(!check_normality && !check_variance_equal) "Normality + Homoscedasticity" else "Homoscedasticity",
             "\n\tSTOPPING assumption checks (Linearity/Slopes homogeneity NOT tested)\n",
             "\tSWITCHING to: Robust regression"),
      paste0("  ==> VIOLATIONS CRITIQUES détectées : ",
             if(!check_normality && !check_variance_equal) "Normalité + Homoscédasticité" else "Homoscédasticité",
             "\n\tARRÊT contrôles assomptions (Linéarité/Homogénéité pentes NON testées)\n",
             "\tPASSAGE vers : Régression robuste"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    goto_robust_ancova <- TRUE
  }

  } # Fin if (!goto_robust_ancova) pour assumption 3

  # ==========================================================================
  # RETEST NORMALITÉ EN MODE "EXTREM" (si normalité échouée MAIS homoscédasticité ok)
  # ==========================================================================
  # Référence: Blanca et al. (2017) DOI:10.7334/psicothema2016.383
  # Standard académique: Donner une seconde chance avec critères stricts (skewness/kurtosis)

  if (!check_normality && check_variance_equal && !goto_robust_ancova) {
    k <- .vbse(
      paste0("NORMALITY RETEST: Applying STRICT criteria\n",
             "\tContext: Non-normality detected BUT homoscedasticity confirmed\n",
             "\tSecond chance: Testing with strict skewness/kurtosis thresholds..."),
      paste0("RE-TEST NORMALITÉ : Application critères STRICTS\n",
             "\tContexte : Non-normalité détectée MAIS homoscédasticité confirmée\n",
             "\tSeconde chance : Test avec seuils stricts skewness/kurtosis..."),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # Appel .normality() en mode tolerance="extrem"
    normality_retest <- .normality(
      residuals_model,
      g = NULL,  # Résidus globaux, pas par groupe
      alpha = alpha,
      tolerance = "extrem",  # Critères stricts (skewness/kurtosis stricts)
      verbose = verbose, code = code,
      k = k,
      debug = debug,
      cpt = "off"
    )

    k <- normality_retest[[2]]
    retest_passed <- normality_retest[[1]]

    if (retest_passed) {
      # Normalité passée avec critères stricts → upgrade à paramétrique
      check_normality <- TRUE

      k <- .vbse(
        paste0("  ==> STRICT normality test PASSED\n",
               "\tResiduals acceptable under strict criteria.\n",
               "\tProceeding with parametric ANCOVA."),
        paste0("  ==> Test normalité strict RÉUSSI\n",
               "\tRésidus acceptables selon critères stricts.\n",
               "\tPoursuite avec ANCOVA paramétrique."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      assumptions_checked$normality$retest_extrem <- TRUE
      assumptions_checked$normality$retest_passed <- TRUE
      assumptions_checked$normality$passed <- TRUE

    } else {
      # Normalité échouée même avec critères stricts → confirme passage robuste
      k <- .vbse(
        paste0("  ==> STRICT normality test FAILED\n",
               "\tResiduals severely non-normal even under lenient criteria.\n",
               "\tWill switch to robust ANCOVA after completing all checks."),
        paste0("  ==> Test normalité strict ÉCHOUÉ\n",
               "\tRésidus sévèrement non-normaux même selon critères tolérants.\n",
               "\tPassage vers ANCOVA robuste après complétion des contrôles."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      assumptions_checked$normality$retest_extrem <- TRUE
      assumptions_checked$normality$retest_passed <- FALSE
    }
  }

  # ==========================================================================
  # ASSOMPTION 4: LINÉARITÉ (DV ~ COVARIABLE)
  # ==========================================================================
  # Référence: Huitema (2011), pp. 89-92
  # Initialiser AVANT le bloc conditionnel pour éviter erreurs si skip
  check_linearity <- TRUE  # Par défaut, on assume linéarité
  linearity_results <- list()
  check_slopes_homogeneity <- TRUE  # Par défaut, on assume pentes homogènes
  slopes_results <- list()

  # Skip if already violated
  if (!goto_robust_ancova) {

  # Réinitialiser pour le test
  check_linearity <- TRUE
  linearity_results <- list()
  linearity_details <- c()

  for (cov in numeric_vars) {
    # Construire formule avec terme quadratique
    formula_quad_str <- paste(
      as.character(formula)[2],
      "~",
      paste(c(factor_vars, cov, paste0("I(", cov, "^2)")), collapse = " + ")
    )
    formula_quad <- as.formula(formula_quad_str)

    model_quad <- tryCatch({
      lm(formula_quad, data = data)
    }, error = function(e) NULL)

    if (!is.null(model_quad)) {
      quad_term <- paste0("I(", cov, "^2)")
      if (quad_term %in% rownames(summary(model_quad)$coefficients)) {
        p_quad <- summary(model_quad)$coefficients[quad_term, "Pr(>|t|)"]

        linearity_results[[cov]] <- list(
          p_value = p_quad,
          linear = p_quad >= alpha
        )

        if (p_quad < alpha) {
          linearity_details <- c(linearity_details,
                                 paste0("'", cov, "': quadratic p = ", .format_pval(p_quad), " (non-linear)"))
          check_linearity <- FALSE
        } else {
          linearity_details <- c(linearity_details,
                                 paste0("'", cov, "': quadratic p = ", .format_pval(p_quad), " (linear)"))
        }
      }
    }
  }

  # Message consolidé: Test + Résultats
  k <- .vbse(
    paste0("ASSUMPTION 4/", n_assumptions, ": Linearity of DV ~ covariate relationship\n",
           "\tTest: Include quadratic term I(covariate²) and check significance [lm() with I(covariate^2)].\n",
           "\tIf quadratic significant => non-linear relationship.\n",
           paste0("\t  ", paste(linearity_details, collapse = "\n\t  ")), "\n",
           if(!check_linearity) {
             "\t==> WARNING: Non-linear relationship detected\n\t    Recommendation: Transform covariate or use polynomial ANCOVA"
           } else {
             "\t==> Linear relationships confirmed"
           }),
    paste0("ASSOMPTION 4/", n_assumptions, " : Linéarité de la relation Variable Dépendante ~ covariable\n",
           "\tTest : Inclure terme quadratique I(covariable²) et vérifier significativité [lm() avec I(covariable^2)].\n",
           "\tSi terme quadratique significatif => relation non-linéaire.\n",
           paste0("\t  ", paste(linearity_details, collapse = "\n\t  ")), "\n",
           if(!check_linearity) {
             "\t==> ATTENTION : Relation non-linéaire détectée\n\t    Considérer : transformation covariable ou ANCOVA polynomiale"
           } else {
             "\t==> Relations linéaires confirmées"
           }),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  assumptions_checked$linearity <- list(
    passed = check_linearity,
    details = linearity_results
  )

  # CODE: Étape 8 - Linéarité DV ~ covariable
  if (isTRUE(code)) {
    lin_lines <- c("# Tester linéarité en incluant terme quadratique")
    for (cov in numeric_vars) {
      lin_lines <- c(lin_lines,
        paste0("# Test pour ", cov, ":"),
        paste0("model_quad_", cov, " <- lm(", deparse(formula[[2]]), " ~ ", cov, " + I(", cov, "^2) + ", paste(factor_vars, collapse = " + "), ", data = dt)"),
        paste0("lin_test_", cov, " <- anova(ancova_model, model_quad_", cov, ")"),
        paste0("# Si p < 0.05 => terme quadratique significatif (relation non-linéaire)")
      )
    }
    k_code <- k_code + 1
    .code_ancova(k_code, "Test de linéarité", lin_lines)
  }

  # ==========================================================================
  # ASSOMPTION 5: HOMOGÉNÉITÉ DES PENTES DE RÉGRESSION (UNIQUEMENT ANCOVA CLASSIQUE)
  # ==========================================================================
  # Référence: Maxwell et al. (2018), pp. 405-409
  # NOTE IMPORTANTE : Ce test N'EST PAS nécessaire si l'interaction facteur:covariable
  # est déjà INCLUSE dans le modèle (ex: y ~ F * cov). Dans ce cas, le modèle
  # PERMET explicitement des pentes différentes (modèle linéaire général).
  # (Initialisation déjà faite avant le bloc conditionnel)

  if (has_factor_cov_interaction) {
    # L'interaction est DANS le modèle => pas de test d'homogénéité requis
    # Ce modèle est à 4 ASSOMPTIONS (pas d'assomption homogénéité pentes)
    k <- .vbse(
      paste0("Homogeneity of regression slopes - NOT REQUIRED\n",
             "\t==> Factor:covariate interaction in model (heterogeneous slopes model)\n",
             "\t    This is a 4-assumption model (no slopes homogeneity assumption)\n",
             "\t    Slopes can differ across groups - this is the model's purpose"),
      paste0("Homogénéité pentes de régression - NON REQUISE\n",
             "\t==> Interaction facteur:covariable dans modèle (modèle pentes hétérogènes)\n",
             "\t    C'est un modèle à 4 assomptions (pas d'assomption homogénéité pentes)\n",
             "\t    Les pentes peuvent différer selon groupes - c'est le but du modèle"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    assumptions_checked$slopes_homogeneity <- list(
      tested = FALSE,
      note = "Not applicable - heterogeneous slopes model (4 assumptions only)"
    )

  } else {
    # ANCOVA CLASSIQUE (sans interaction) => tester homogénéité (5ème assomption)
    slopes_details <- c()

  for (f in factor_vars) {
    for (cov in numeric_vars) {
      # Construire formule avec interaction factor:covariate
      formula_int_str <- paste(
        as.character(formula)[2],
        "~",
        paste(c(factor_vars, numeric_vars, paste0(f, ":", cov)), collapse = " + ")
      )
      formula_int <- as.formula(formula_int_str)

      model_int <- tryCatch({
        lm(formula_int, data = data)
      }, error = function(e) NULL)

      if (!is.null(model_int)) {
        anova_int <- anova(model_int)
        int_term <- paste0(f, ":", cov)

        if (int_term %in% rownames(anova_int)) {
          p_int <- anova_int[int_term, "Pr(>F)"]

          slopes_results[[paste0(f, "_x_", cov)]] <- list(
            p_value = p_int,
            homogeneous = p_int >= alpha
          )

          if (p_int < alpha) {
            slopes_details <- c(slopes_details,
                               paste0(f, " × ", cov, ": p = ", .format_pval(p_int), " (heterogeneous)"))
            check_slopes_homogeneity <- FALSE
          } else {
            slopes_details <- c(slopes_details,
                               paste0(f, " × ", cov, ": p = ", .format_pval(p_int), " (homogeneous)"))
          }
        }
      }
    }
  }

  # Message consolidé
  k <- .vbse(
    paste0("ASSUMPTION 5/", n_assumptions, ": Homogeneity of regression slopes\n",
           "\tTest: factor:covariate interaction must be NON-significant [lm() + anova()]\n",
           "\tIf significant => different slopes across groups (violation).\n",
           paste0("\t  ", paste(slopes_details, collapse = "\n\t  ")), "\n",
           if(!check_slopes_homogeneity) {
             "\t==> WARNING: Heterogeneous slopes detected\n\t    Recommendation: Use Johnson-Neyman technique or heterogeneous slopes model"
           } else {
             "\t==> Homogeneous slopes confirmed (parallel regression lines)"
           }),
    paste0("ASSOMPTION 5/", n_assumptions, " : Homogénéité des pentes de régression\n",
           "\tTest : interaction facteur:covariable doit être NON-significative [lm() + anova()]\n",
           "\tSi significative => pentes de régression différentes entre groupes (violation).\n",
           paste0("\t  ", paste(slopes_details, collapse = "\n\t  ")), "\n",
           if(!check_slopes_homogeneity) {
             "\t==> ATTENTION : Pentes hétérogènes détectées\n\t    Considérer : technique Johnson-Neyman ou modèle pentes hétérogènes"
           } else {
             "\t==> Pentes homogènes confirmées (droites régression parallèles)"
           }),
    verbose = verbose, code = code, k = k, cpt = "on"
  )

  assumptions_checked$slopes_homogeneity <- list(
    passed = check_slopes_homogeneity,
    details = slopes_results
  )

  # CODE: Étape 9 - Homogénéité pentes de régression
  if (isTRUE(code) && !has_factor_cov_interaction) {
    slopes_lines <- c("# Tester si pentes de régression diffèrent selon groupes facteurs")
    for (f in factor_vars) {
      for (cov in numeric_vars) {
        slopes_lines <- c(slopes_lines,
          paste0("# Test interaction ", f, " × ", cov, ":"),
          paste0("model_interaction <- lm(", deparse(formula[[2]]), " ~ ", paste(c(factor_vars, numeric_vars), collapse = " + "),
                 " + ", f, ":", cov, ", data = dt)"),
          paste0("slopes_test <- anova(ancova_model, model_interaction)"),
          paste0("# Si p < 0.05 => pentes hétérogènes (violation assomption ANCOVA)")
        )
      }
    }
    k_code <- k_code + 1
    .code_ancova(k_code, "Test homogénéité pentes de régression", slopes_lines)
  }

  } # Fin else (ANCOVA classique, test homogénéité requis)

  } # Fin if (!goto_robust_ancova) pour assumptions 4 et 5

  # ==========================================================================
  # DÉCISION: PARAMÉTRIQUE VS ROBUSTE
  # ==========================================================================

  # Si early stop activé, forcer robust
  if (goto_robust_ancova) {
    all_assumptions_met <- FALSE
  } else {
    all_assumptions_met <- (
    check_normality &&
    check_variance_equal &&
    check_linearity &&
    check_slopes_homogeneity
    )
  }

  # ==========================================================================
  # CONTRÔLE TOLÉRANT DE LA NORMALITÉ (si seule violation)
  # ==========================================================================
  # Si SEULE la normalité est violée ET variances homogènes ET n≥30 par groupe
  # → Le Théorème Central Limite (TCL) garantit la robustesse des tests F
  #
  # Références académiques:
  # - Lumley et al. (2002) DOI:10.1146/annurev.publhealth.23.100901.140546
  # - Blanca et al. (2017) DOI:10.7334/psicothema2016.383
  # - Maxwell et al. (2018), pp. 133-135

  check_normality_tolerant <- FALSE

  if (!all_assumptions_met &&
      !check_normality &&
      check_variance_equal &&
      check_linearity &&
      check_slopes_homogeneity) {

    # Vérifier tailles de groupes
    group_sizes <- table(g_cat)

    if (all(group_sizes >= 30)) {
      k <- .vbse(
        paste0("TOLERANT NORMALITY CHECK: Central Limit Theorem application\n",
               "\tContext: Only normality violated, but homoscedasticity confirmed + all group sizes ≥ 30\n",
               "\tCLT guarantees F-test robustness despite non-normal residuals\n",
               "\t==> Proceeding with parametric ANCOVA"),
        paste0("CONTRÔLE TOLÉRANT DE LA NORMALITÉ : Application du Théorème Central Limite\n",
               "\tContexte : Seule la normalité violée, mais homoscédasticité confirmée + tailles de groupes ≥ 30\n",
               "\tLe TCL garantit la robustesse du test F malgré résidus non normaux\n",
               "\t==> Poursuite avec ANCOVA paramétrique"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      check_normality_tolerant <- TRUE
      all_assumptions_met <- TRUE  # Permettre ANCOVA paramétrique
      assumptions_checked$normality$tolerant_applied <- TRUE
      assumptions_checked$normality$group_sizes <- as.vector(group_sizes)
    }
  }

  if (all_assumptions_met) {
    # ANCOVA PARAMÉTRIQUE
    # NOTE: Le type de Sommes des Carrés a déjà été sélectionné à l'étape 3
    # via .select_ss_type() et warnings affichés si nécessaire

    k <- .vbse(
      paste0("Proceeding with parametric ANCOVA (assumptions met)\n",
             "\t==> Displaying results:"),
      paste0("Poursuite avec ANCOVA paramétrique (assomptions respectées)\n",
             "\t==> Affichage du bilan :"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # CODE: Étape 10 - ANCOVA paramétrique
    if (isTRUE(code)) {
      if (selected_type == "I") {
        k_code <- k_code + 1
        .code_ancova(k_code, "ANCOVA paramétrique - Type I", c(
          "# Type I: Sommes des carrés séquentielles",
          "anova_result <- anova(ancova_model)",
          "print(anova_result)"
        ))
      } else {
        k_code <- k_code + 1
        .code_ancova(k_code, paste0("ANCOVA paramétrique - Type ", selected_type), c(
          paste0("# Type ", selected_type, ": Sommes des carrés marginales"),
          "library(car)",
          paste0("anova_result <- car::Anova(ancova_model, type = '", selected_type, "')"),
          "print(anova_result)"
        ))
      }
    }

    if (verbose) {
      # Type I utilise anova() de base, Type II/III utilisent car::Anova()
      if (selected_type == "I") {
        # Type I: Sommes des carrés séquentielles (anova de base)
        anova_result <- anova(model)
        cat("\n")
        print(anova_result)
        cat("\n")

      } else if (selected_type %in% c("II", "III") && requireNamespace("car", quietly = TRUE)) {
        # Type II/III: Tentative avec car::Anova (peut échouer si coefficients aliasés)
        anova_result <- tryCatch({
          car::Anova(model, type = selected_type)
        }, error = function(e) {
          if (grepl("aliased", conditionMessage(e), ignore.case = TRUE)) {
            # Coefficients aliasés (typique avec interactions facteur×facteur et contr.sum)
            # Basculer vers Type II si on était en Type III
            if (selected_type == "III") {
              k <<- .vbse(
                paste0("\tNote: Aliased coefficients detected (factor×factor interactions with sum contrasts)\n",
                       "\tSwitching to Type II Sum of Squares [car::Anova(type='II')]\n",
                       "\tType II: Each effect adjusted for all others at same or lower order"),
                paste0("\tNote : Coefficients aliasés détectés (interactions facteur×facteur avec contrastes somme)\n",
                       "\tPassage vers Sommes des Carrés Type II [car::Anova(type='II')]\n",
                       "\tType II : Chaque effet ajusté pour tous les autres de même ordre ou inférieur"),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
              car::Anova(model, type = "II")
            } else {
              stop(e)  # Re-throw si déjà Type II
            }
          } else {
            stop(e)  # Re-throw si autre erreur
          }
        })
        cat("\n")
        print(anova_result)
        cat("\n")
      }

      # CONCLUSION: Identifier effets significatifs (distinguer facteurs vs covariables)
      pval_col <- which(colnames(anova_result) %in% c("Pr(>F)", "p.value", "P(>|F|)", "Pr(>Chisq)"))
      if (length(pval_col) > 0) {
        row_names <- rownames(anova_result)
        valid_rows <- which(!grepl("Residuals|^\\(Intercept\\)", row_names, ignore.case = TRUE))
        if (length(valid_rows) > 0) {
          sig_effects <- row_names[valid_rows][anova_result[valid_rows, pval_col[1]] < alpha]
          if (length(sig_effects) > 0) {
            # Séparer facteurs et covariables
            sig_factors <- sig_effects[sig_effects %in% factor_vars]
            sig_covariates <- sig_effects[sig_effects %in% numeric_vars]

            # Message détaillé
            if (length(sig_factors) > 0 && length(sig_covariates) > 0) {
              k <- .vbse(
                paste0("\tSignificant effects detected:\n",
                       "\t  Factors: ", paste(sig_factors, collapse = ", "), "\n",
                       "\t  Covariates: ", paste(sig_covariates, collapse = ", ")),
                paste0("\tEffets significatifs détectés :\n",
                       "\t  Facteurs : ", paste(sig_factors, collapse = ", "), "\n",
                       "\t  Covariables : ", paste(sig_covariates, collapse = ", ")),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            } else if (length(sig_factors) > 0) {
              k <- .vbse(
                paste0("\tSignificant effects detected:\n",
                       "\t  Factors: ", paste(sig_factors, collapse = ", ")),
                paste0("\tEffets significatifs détectés :\n",
                       "\t  Facteurs : ", paste(sig_factors, collapse = ", ")),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            } else if (length(sig_covariates) > 0) {
              k <- .vbse(
                paste0("\tSignificant effects detected:\n",
                       "\t  Covariates: ", paste(sig_covariates, collapse = ", ")),
                paste0("\tEffets significatifs détectés :\n",
                       "\t  Covariables : ", paste(sig_covariates, collapse = ", ")),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            } else {
              # Interactions ou autres effets
              k <- .vbse(
                paste0("\tSignificant effects detected: ", paste(sig_effects, collapse = ", ")),
                paste0("\tEffets significatifs détectés : ", paste(sig_effects, collapse = ", ")),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
          } else {
            k <- .vbse(
              paste0("\tNo significant effects at alpha = ", alpha),
              paste0("\tAucun effet significatif à alpha = ", alpha),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          }
        }
      }
    }

    return(list(
      model = model,
      assumptions_checked = assumptions_checked,
      robust = FALSE,
      g_cat = g_cat,
      k = k
    ))

  } else {
    # ANCOVA ROBUSTE
    # Ne pas afficher de message ici si early stop déjà activé
    # (le message a déjà été affiché lors de l'early stop)
    if (!goto_robust_ancova) {
      # Violations détectées APRÈS early stop (linéarité ou pentes)
      k <- .vbse(
        paste0("==> ASSUMPTION VIOLATIONS DETECTED: Switching to ROBUST ANCOVA\n",
               "\tViolations: ",
               if(!check_normality) "Normality " else "",
               if(!check_variance_equal) "Homoscedasticity " else "",
               if(!check_linearity) "Linearity " else "",
               if(!check_slopes_homogeneity) "Slopes_homogeneity " else ""),
        paste0("==> VIOLATIONS D'ASSOMPTIONS DÉTECTÉES : Passage vers ANCOVA ROBUSTE\n",
               "\tViolations : ",
               if(!check_normality) "Normalité " else "",
               if(!check_variance_equal) "Homoscédasticité " else "",
               if(!check_linearity) "Linéarité " else "",
               if(!check_slopes_homogeneity) "Homogénéité_pentes " else ""),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }

    # Appeler logique robuste ANCOVA (déjà implémentée en session 15)
    n_categorical_factors <- length(factor_vars)
    n_continuous_covariates <- length(numeric_vars)

    robust_results <- list(
      method = NULL,
      test_result = NULL,
      assumptions_checked = assumptions_checked,
      warnings = character(),
      posthoc_applicable = FALSE
    )

    if (n_categorical_factors == 1 && n_continuous_covariates >= 1) {
      # Régression robuste
      # RÉFÉRENCE ACADÉMIQUE (développeurs/documentation uniquement - bp.log 7.4.6.1):
      # Wilcox, R. R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.).
      # Academic Press. ISBN: 978-0128047330. Chapitre 7 (pp. 423-456): Robust ANCOVA methods.

      k <- .vbse(
        paste0("Robust method selected: Robust regression [rlm() {MASS}]\n",
               "\tJustification: 1 factor + covariate(s) + assumption violations\n",
               "\tMethod: MM-estimation (resistant to outliers and violations)\n",
               "\t(Alternatives: Quade test if paired data, GLM if specific distribution)\n",
               "\t==> Fitting robust ANCOVA model..."),
        paste0("Sélection méthode robuste : Régression robuste [rlm() {MASS}]\n",
               "\tJustification : 1 facteur + covariable(s) + violations assomptions\n",
               "\tMéthode : MM-estimation (résistant outliers et violations)\n",
               "\t(Alternatives : Test Quade si paired, GLM si distribution spécifique)\n",
               "\t==> Ajustement du modèle ANCOVA robuste..."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      if (requireNamespace("MASS", quietly = TRUE)) {
        robust_model <- tryCatch({
          MASS::rlm(formula, data = data, method = "MM")
        }, error = function(e) {
          k <<- .vbse(
            paste0("Error: ", e$message),
            paste0("Erreur : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          return(NULL)
        })

        if (!is.null(robust_model)) {
          if (verbose) {
            # IMPORTANT: summary.rlm() peut échouer si résidus contiennent NA
            # Wrapper dans tryCatch() pour éviter crash
            tryCatch({
              print(summary(robust_model))
              cat("\n")
            }, error = function(e) {
              cat("Robust regression model fitted (summary unavailable due to NA residuals)\n\n")
            })
          }
          robust_results$method <- "Robust_Regression_RLM"
          robust_results$model <- robust_model
          robust_results$posthoc_applicable <- TRUE
        }
      }

    } else if (n_categorical_factors >= 2 && n_continuous_covariates >= 1) {
      #==========================================================================
      # CAS 2: PLUSIEURS FACTEURS + COVARIABLES → White's Heteroscedasticity-Consistent SEs
      #==========================================================================
      # NOTE ACADÉMIQUE: aovp() non adapté pour ANCOVA car permute observations entières,
      # ce qui brise la relation covariate-DV. Solution: erreurs standard robustes (HC3).
      #
      # Référence: Long & Ervin (2000). Using heteroscedasticity consistent standard errors.
      #            DOI:10.1080/00031305.2000.10474549

      k <- .vbse(
        paste0("Method selected: ANCOVA with heteroscedasticity-robust standard errors\n",
               "\tApproach: White's correction (HC3 estimator)\n",
               "\tReason: Permutation tests invalid for ANCOVA (break covariate-DV relationship)\n",
               "\tRationale: HC3 provides valid inference without normality/homoscedasticity assumptions"),
        paste0("Méthode sélectionnée : ANCOVA avec erreurs standard robustes à l'hétéroscédasticité\n",
               "\tApproche : Correction de White (estimateur HC3)\n",
               "\tRaison : Tests de permutation invalides pour ANCOVA (brisent relation covariable-Variable Dépendante)\n",
               "\tRationale : HC3 fournit inférence valide sans hypothèses de normalité/homoscédasticité"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      if (requireNamespace("sandwich", quietly = TRUE) &&
          requireNamespace("lmtest", quietly = TRUE)) {

        # Test avec erreurs standard robustes HC3
        robust_vcov <- tryCatch({
          sandwich::vcovHC(model, type = "HC3")
        }, error = function(e) {
          k <<- .vbse(
            paste0("Error computing robust covariance matrix: ", e$message),
            paste0("Erreur calcul matrice covariance robuste : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          return(NULL)
        })

        if (!is.null(robust_vcov)) {
          # Tests de Wald avec SEs robustes (supprimer messages anglais du package car)
          robust_tests <- tryCatch({
            suppressMessages(car::Anova(model, vcov. = robust_vcov, type = "III"))
          }, error = function(e) {
            if (grepl("aliased", conditionMessage(e), ignore.case = TRUE)) {
              # Coefficients aliasés avec robust vcov → car::Anova() échoue même avec Type II
              # Tentative Type II d'abord
              type_ii_result <- tryCatch({
                suppressMessages(car::Anova(model, vcov. = robust_vcov, type = "II"))
              }, error = function(e2) {
                if (grepl("aliased", conditionMessage(e2), ignore.case = TRUE)) {
                  # Type II échoue aussi → basculer vers rlm() comme fallback
                  k <<- .vbse(
                    paste0("\tNote: Aliased coefficients incompatible with robust standard errors\n",
                           "\tFalling back to robust regression [rlm() {MASS}]"),
                    paste0("\tNote : Coefficients aliasés incompatibles avec erreurs standard robustes\n",
                           "\tRetour vers régression robuste [rlm() {MASS}]"),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )

                  # Appeler rlm() comme fallback
                  if (requireNamespace("MASS", quietly = TRUE)) {
                    robust_model <- MASS::rlm(formula, data = data, method = "MM")
                    robust_results$method <<- "Robust_Regression_RLM_Fallback"
                    robust_results$model <<- robust_model
                    robust_results$posthoc_applicable <<- TRUE

                    if (verbose) {
                      tryCatch({
                        print(summary(robust_model))
                        cat("\n")
                      }, error = function(e3) {
                        cat("Robust regression model fitted (summary unavailable due to NA residuals)\n\n")
                      })
                    }
                  }
                  return(NULL)  # Pas de robust_tests, utilise rlm() à la place
                } else {
                  stop(e2)  # Autre erreur, propager
                }
              })

              # Si Type II réussit, l'utiliser
              if (!is.null(type_ii_result)) {
                k <<- .vbse(
                  paste0("\tNote: Aliased coefficients detected, using Type II"),
                  paste0("\tNote : Coefficients aliasés détectés, utilisation Type II"),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
                return(type_ii_result)
              } else {
                return(NULL)
              }
            } else {
              k <<- .vbse(
                paste0("Error computing robust Wald tests: ", conditionMessage(e)),
                paste0("Erreur calcul tests de Wald robustes : ", conditionMessage(e)),
                verbose = verbose, code = code, k = k, cpt = "on"
              )
              return(NULL)
            }
          })

          if (!is.null(robust_tests)) {
            if (verbose) {
              cat("\n=== ANCOVA Type III with Robust Standard Errors (HC3) ===\n\n")
              print(robust_tests)
              cat("\n")

              # CONCLUSION: Identifier effets significatifs (distinguer facteurs vs covariables)
              pval_col <- which(colnames(robust_tests) %in% c("Pr(>F)", "p.value", "P(>|F|)", "Pr(>Chisq)"))
              if (length(pval_col) > 0) {
                row_names <- rownames(robust_tests)
                valid_rows <- which(!grepl("Residuals|^\\(Intercept\\)", row_names, ignore.case = TRUE))
                if (length(valid_rows) > 0) {
                  sig_effects <- row_names[valid_rows][robust_tests[valid_rows, pval_col[1]] < alpha]
                  if (length(sig_effects) > 0) {
                    # Séparer facteurs et covariables
                    sig_factors <- sig_effects[sig_effects %in% factor_vars]
                    sig_covariates <- sig_effects[sig_effects %in% numeric_vars]

                    # Message détaillé
                    if (length(sig_factors) > 0 && length(sig_covariates) > 0) {
                      k <- .vbse(
                        paste0("\tSignificant effects detected:\n",
                               "\t  Factors: ", paste(sig_factors, collapse = ", "), "\n",
                               "\t  Covariates: ", paste(sig_covariates, collapse = ", ")),
                        paste0("\tEffets significatifs détectés :\n",
                               "\t  Facteurs : ", paste(sig_factors, collapse = ", "), "\n",
                               "\t  Covariables : ", paste(sig_covariates, collapse = ", ")),
                        verbose = verbose, code = code, k = k, cpt = "off"
                      )
                    } else if (length(sig_factors) > 0) {
                      k <- .vbse(
                        paste0("\tSignificant effects detected:\n",
                               "\t  Factors: ", paste(sig_factors, collapse = ", ")),
                        paste0("\tEffets significatifs détectés :\n",
                               "\t  Facteurs : ", paste(sig_factors, collapse = ", ")),
                        verbose = verbose, code = code, k = k, cpt = "off"
                      )
                    } else if (length(sig_covariates) > 0) {
                      k <- .vbse(
                        paste0("\tSignificant effects detected:\n",
                               "\t  Covariates: ", paste(sig_covariates, collapse = ", ")),
                        paste0("\tEffets significatifs détectés :\n",
                               "\t  Covariables : ", paste(sig_covariates, collapse = ", ")),
                        verbose = verbose, code = code, k = k, cpt = "off"
                      )
                    } else {
                      # Interactions ou autres effets
                      k <- .vbse(
                        paste0("\tSignificant effects detected: ", paste(sig_effects, collapse = ", ")),
                        paste0("\tEffets significatifs détectés : ", paste(sig_effects, collapse = ", ")),
                        verbose = verbose, code = code, k = k, cpt = "off"
                      )
                    }
                  } else {
                    k <- .vbse(
                      paste0("\tNo significant effects at alpha = ", alpha),
                      paste0("\tAucun effet significatif à alpha = ", alpha),
                      verbose = verbose, code = code, k = k, cpt = "off"
                    )
                  }
                }
              }
            }

            robust_results$method <- "ANCOVA_White_HC3"
            robust_results$test_result <- robust_tests
            robust_results$vcov_robust <- robust_vcov
            robust_results$posthoc_applicable <- TRUE
          }
        }

      } else {
        k <- .vbse(
          paste0("Warning: 'sandwich' or 'lmtest' package not available.\n",
                 "\tRobust ANCOVA cannot be performed.\n",
                 "\tInstall with: install.packages(c('sandwich', 'lmtest'))"),
          paste0("Attention : Package 'sandwich' ou 'lmtest' non disponible.\n",
                 "\tANCOVA robuste impossible.\n",
                 "\tInstaller avec : install.packages(c('sandwich', 'lmtest'))"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        robust_results$method <- "ANCOVA_Robust_Unavailable"
      }
    }

    return(list(
      model = NULL,
      assumptions_checked = assumptions_checked,
      robust = TRUE,
      robust_results = robust_results,
      g_cat = g_cat,
      k = k
    ))
  }
}

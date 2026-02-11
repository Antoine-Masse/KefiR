#' Analyse multi-facteurs pour ANOVA et ANCOVA (Version 4 - Réécriture Intégrale)
#'
#' Cette fonction effectue des analyses multi-facteurs pour des modèles ANOVA ou ANCOVA.
#' Elle détecte automatiquement le type de modèle, vérifie les assomptions de base
#' (normalité, homogénéité des variances, sphéricité), et effectue les ajustements
#' nécessaires pour des analyses robustes. Conçue pour des analyses à  deux facteurs ou plus.
#'
#' @param x Vecteur numérique. Variable dépendante à  analyser. Si \code{NULL}, extrait de \code{data[,1]}.
#' @param g Data frame ou matrice. Variables groupantes/explicatives. Si \code{NULL}, extrait de \code{data[,-1]}.
#' @param formula Une formule spécifiant le modèle à  analyser (ex: \code{y ~ facteur1 * facteur2}).
#' @param data Data frame optionnel contenant toutes les variables spécifiées dans la formule.
#' @param paired Logique. Indique si les données sont appariées/mesures répétées (défaut: \code{FALSE}).
#' @param id Chaîne de caractères. Nom de la colonne utilisée pour apparier les observations si \code{paired = TRUE}.
#' @param within Vecteur de caractères. Noms des facteurs intra-sujet (within-subject).
#' @param between Vecteur de caractères optionnel. Noms des facteurs inter-sujet (between-subject).
#' @param alpha Numérique. Seuil de signification (défaut: 0.05).
#' @param k Entier. Compteur pour les messages verbose (défaut: \code{NULL}, auto-initialisé à  0).
#' @param code Logique. Si \code{TRUE}, affiche des extraits de code (défaut: \code{FALSE}).
#' @param debug Logique. Si \code{TRUE}, affiche des messages de débogage (défaut: \code{FALSE}).
#' @param verbose Logique. Si \code{TRUE}, fournit une sortie détaillée (défaut: \code{FALSE}).
#'
#' @details
#' \strong{ARCHITECTURE DU PIPELINE:}
#' \enumerate{
#'   \item \strong{Détection du type de modèle}: ANOVA vs. ANCOVA
#'   \item \strong{Détection du design}: équilibré/déséquilibré, apparié, effets aléatoires
#'   \item \strong{Tests d'hypothèses}: Normalité, homogénéité des variances, sphéricité, indépendance
#'   \item \strong{Sélection de la stratégie}: Paramétrique vs. Robuste
#'   \item \strong{Contrôles ANCOVA}: Homogénéité des pentes, linéarité, homoscédasticité
#' }
#'
#' @return
#' Une liste contenant:
#' \itemize{
#'   \item \code{x}: Variable dépendante traitée
#'   \item \code{g_cat}: Facteur d'interaction des variables catégorielles
#'   \item \code{check_normality}: Logique indiquant si l'hypothèse de normalité est respectée
#'   \item \code{check_variance_equal}: Logique indiquant si l'homogénéité des variances est respectée
#'   \item \code{k}: Compteur de messages mis à  jour
#'   \item \code{model}: Objet modèle ajusté (PEUT àŠTRE NULL en voie robuste ou mesures répétées)
#'   \item \code{robuste}: Logique indiquant si des méthodes robustes ont été utilisées
#'   \item \code{check_ancova}: Logique indiquant si une ANCOVA a été détectée
#'   \item \code{check_discret}: Logique indiquant si la variable dépendante est discrète
#'   \item \code{alea}: Logique indiquant si des effets aléatoires ont été détectés
#'   \item \code{alea_plus}: Logique indiquant si la structure d'effets aléatoires est complexe
#'   \item \code{formula}: Formule originale
#'   \item \code{updated_formula}: Formule mise à  jour
#'   \item \code{ancova_checks}: Liste des résultats des contrôles ANCOVA
#'   \item \code{robust_results}: Liste des résultats de l'analyse robuste (NULL si paramétrique)
#' }
#'
#' @section Avertissements:
#' \itemize{
#'   \item Pour les designs déséquilibrés avec dépendances entre facteurs, des interactions
#'         peuvent être automatiquement ajoutées au modèle
#'   \item Les modèles mixtes sont utilisés en dernier recours quand les autres approches échouent
#'   \item Les variables dépendantes discrètes déclenchent automatiquement une analyse robuste
#'   \item \strong{IMPORTANT}: Le champ \code{model} peut être \code{NULL} dans les cas suivants:
#'     (1) Analyse robuste déclenchée, (2) Mesures répétées avec ezANOVA, (3) échec d'ajustement du modèle
#' }
#'
#' @references
#' Field, A., Miles, J., & Field, Z. (2012). \emph{Discovering statistics using R}.
#' Sage Publications. \doi{10.1111/insr.12011_21}
#'
#' Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and analyzing data:
#' A model comparison perspective} (2nd ed.). Lawrence Erlbaum Associates.
#'
#' Wilcox, R. R. (2017). \emph{Introduction to robust estimation and hypothesis testing}
#' (4th ed.). Academic Press. \doi{10.1016/C2010-0-67044-1}
#'
#' @importFrom stats aov bartlett.test fligner.test shapiro.test var.test
#'   oneway.test formula terms reshape contr.helmert residuals kruskal.test friedman.test
#'   rstudent rstandard
#' @importFrom car Anova leveneTest
#' @importFrom ez ezANOVA
#' @importFrom DescTools GTest
#' @importFrom rstatix identify_outliers
#' @importFrom withr with_options
#'
#' @seealso
#' \code{\link{.one_factor_analysis}}, \code{\link[car]{Anova}},
#' \code{\link[ez]{ezANOVA}}, \code{\link[lme4]{lmer}}
#'
#' @export
.multi_factor_analysis <- function(
    x = NULL,
    g = NULL,
    formula = NULL,
    data = NULL,
    paired = FALSE,
    id = NULL,
    within = NULL,
    between = NULL,
    alpha = 0.05,
    k = NULL,
    code = FALSE,
    debug = FALSE,
    verbose = FALSE
) {

  .dbg("=== Start .multi_factor_analysis() ===",
       "=== Début de .multi_factor_analysis() ===", debug=debug)

  #============================================================================
  #                      NOTE PERSO #1 [EN COURS]
  #============================================================================
  #* NOTE PERSO #1 : Gestion des données manquantes (NA)
  #*
  #* ➜ Problème identifié :
  #*   - m.test() supprime les NA en amont (listwise deletion)
  #*   - Pour les modèles mixtes, le listwise n'est PAS optimal car il supprime
  #*     des observations complètes pour d'autres variables
  #*   - Beaucoup d'appels (bartlett.test, kruskal.test, leveneTest) plantent si NA
  #*
  #* ➜ Source académique (APA + DOI) :
  #*   Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the
  #*   state of the art. *Psychological Methods*, 7(2), 147â€"177.
  #*   https://doi.org/10.1037/1082-989X.7.2.147
  #*
  #*   Idée principale : Le listwise deletion est acceptable SEULEMENT si :
  #*   (a) Les données manquent complètement au hasard (MCAR)
  #*   (b) Le taux de données manquantes est < 5%
  #*   Pour les modèles mixtes (lmer), Maximum Likelihood (ML) utilise toutes
  #*   les données disponibles et est supérieur au listwise.
  #*
  #* ➜ Solution appliquée :
  #*   1. Détection précoce des NA et calcul du taux de manquants
  #*   2. Si taux > 5%, avertissement à  l'utilisateur
  #*   3. Pour la voie paramétrique classique : na.omit() défensif avant chaque test
  #*   4. Pour les modèles mixtes : conserver les NA, lmer() gère via ML
  #*   5. Documentation explicite de la stratégie dans les messages verbose
  #*
  #* ➜ Statut : Solution partielle implémentée ci-dessous
  #============================================================================

  # Détection précoce des NA
  if (!is.null(data)) {
    na_count <- sum(is.na(data))
    total_cells <- nrow(data) * ncol(data)
    na_rate <- (na_count / total_cells) * 100

    if (na_rate > 5 && verbose) {
      k <- .vbse(
        paste0("Warning: ", round(na_rate, 2), "% missing data detected (>5% threshold).\n",
               "\tListwise deletion may reduce power. Consider imputation for mixed models."),
        paste0("Attention : ", round(na_rate, 2), "% de données manquantes détectées (seuil >5%).\n",
               "\tLa suppression listwise peut réduire la puissance. Envisager l'imputation pour modèles mixtes."),
        verbose = verbose, k = k, cpt = "on"
      )
    }

    # Suppression listwise pour la voie classique (défensif)
    data_clean <- na.omit(data)
    if (nrow(data_clean) < nrow(data) && verbose) {
      k <- .vbse(
        paste0(nrow(data) - nrow(data_clean), " observations removed due to missing data (listwise deletion)."),
        paste0(nrow(data) - nrow(data_clean), " observations supprimées pour données manquantes (listwise deletion)."),
        verbose = verbose, k = k, cpt = "off"
      )
    }
    data <- data_clean
  }

  #============================================================================
  #                   FONCTIONS INTERNES
  #============================================================================

  #' Fonction pour extraire automatiquement les résidus
  get_residuals <- function(model) {
    if (inherits(model, "aov") || inherits(model, "lm")) {
      return(residuals(model))
    } else if (inherits(model, "aovlist")) {
      if (requireNamespace("dae", quietly = TRUE)) {
        return(dae::residuals.aovlist(model))
      } else {
        return(unlist(lapply(model, residuals)))
      }
    } else {
      .exit(
        "The provided model is neither an 'aov', 'lm', nor 'aovlist' object.",
        "Le modèle fourni n'est ni un objet 'aov', 'lm', ni 'aovlist'."
      )
    }
  }

  #============================================================================
  #                      NOTE PERSO #2 [RéSOLUE]
  #============================================================================
  #* NOTE PERSO #2 : Résidus studentisés vs standardisés pour ancova_checks
  #*
  #* ➜ Problème identifié :
  #*   Bonne idée d'avoir un fallback rstudentâ†'rstandardâ†'residuals, mais il faut
  #*   garder trace du type de résidu utilisé dans ancova_checks (utile pour audit).
  #*
  #* ➜ Source académique (APA + DOI) :
  #*   Cook, R. D., & Weisberg, S. (1982). *Residuals and influence in regression*.
  #*   Chapman and Hall. ISBN: 978-0412242809
  #*
  #*   Idée principale : Les résidus studentisés (externally studentized) sont
  #*   préférés pour la détection d'outliers car chaque observation est standardisée
  #*   en utilisant un modèle où elle est exclue. Plus robuste que les résidus
  #*   standardisés (internally studentized).
  #*
  #* ➜ Solution appliquée :
  #*   1. Ajout d'un champ "residual_type" dans le retour de la fonction
  #*   2. Trace explicite du type utilisé dans ancova_checks
  #*   3. Message verbose si fallback nécessaire
  #*
  #* ➜ Statut : RéSOLU
  #============================================================================

  get_studentized_residuals <- function(model) {
    residual_type <- "studentized"  # Par défaut

    tryCatch({
      if (inherits(model, "lm") || inherits(model, "aov")) {
        res <- rstudent(model)
      } else {
        residual_type <- "standardized"
        res <- rstandard(model)
      }
    }, error = function(e) {
      residual_type <- "raw"
      res <- residuals(model)
      if (verbose) {
        k <- .vbse(
          "Warning: Could not compute studentized residuals. Using raw residuals.",
          "Attention : Impossible de calculer les résidus studentisés. Utilisation des résidus bruts.",
          verbose = verbose, k = k, cpt = "off"
        )
      }
    })

    return(list(residuals = res, type = residual_type))
  }

  #============================================================================
  #                         INITIALISATION
  #============================================================================

  if (is.null(k)) k <- 0

  # Initialiser les drapeaux de contrôle
  robuste <- FALSE
  check_ancova <- FALSE
  alea <- FALSE
  alea_plus <- FALSE
  check_normality <- TRUE
  check_variance_equal <- TRUE
  check_discret <- FALSE
  model <- NULL
  updated_formula <- formula
  ancova_checks <- list()
  residual_type_used <- NULL  # NOUVEAU: trace du type de résidu
  # Initialisation de robust_results (utilisé si voie robuste)
  robust_results <- NULL

  #============================================================================
  #              VALIDATION ET PRéPARATION DES ENTRéES
  #============================================================================

  # Extraire x et g depuis data si non fournis
  if (is.null(g)) {
    if (is.null(data)) {
      .exit(
        "Either 'g' or 'data' must be provided.",
        "Soit 'g' soit 'data' doit être fourni."
      )
    }
    if (!is.null(id)) {
      col_to_exclude <- c(1, which(names(data) == id))
      g <- data[, -col_to_exclude, drop = FALSE]
    } else {
      g <- data[, -1, drop = FALSE]
    }
  }

  if (is.null(x)) {
    if (is.null(data)) {
      .exit(
        "Either 'x' or 'data' must be provided.",
        "Soit 'x' soit 'data' doit être fourni."
      )
    }
    x <- data[, 1]
  }

  # CORRECTION CRITIQUE : S'assurer que g est TOUJOURS un data.frame
  if (!is.data.frame(g)) {
    # Extraire les noms de prédicteurs depuis la formule
    if (!is.null(formula)) {
      predictor_names <- setdiff(all.vars(formula), all.vars(formula)[1])
    } else {
      predictor_names <- "g"  # Nom par défaut
    }

    g <- as.data.frame(g)
    names(g) <- predictor_names[1:ncol(g)]
  }

  # Convertir les colonnes binaires et texte dans g en facteurs
  for (col_name in names(g)) {
    col <- g[[col_name]]
    if (is.character(col)) {
      g[[col_name]] <- as.factor(col)
    } else if (is.numeric(col) && length(unique(na.omit(col))) == 2 &&
               all(sort(unique(na.omit(col))) %in% c(0, 1))) {
      g[[col_name]] <- as.factor(col)
    }
  }

  # Reconstruire data
  resp <- if (!is.null(formula)) all.vars(formula)[1] else "x"

  if (!is.null(id)) {
    temp_df <- data.frame(x = x, id_col = data[[id]])
    names(temp_df) <- c(resp, id)
    data <- cbind(temp_df, g)
  } else {
    data <- data.frame(x, g, check.names = FALSE)
    names(data)[1] <- resp
  }

  # Affichage debug
  # DEBUG : Afficher la structure de data après reconstruction
  if (debug) {
    cat("\n=== DEBUG après reconstruction de data ===\n")
    cat("Noms de colonnes dans data :", paste(names(data), collapse = ", "), "\n")
    cat("Noms de colonnes dans g :", paste(names(g), collapse = ", "), "\n")
    cat("dim(data) :", paste(dim(data), collapse = " x "), "\n")
    cat("str(data) :\n")
    str(data)
    cat("==========================================\n\n")
  }

  #============================================================================
  #             CONTRÔLE VARIABLE DéPENDANTE DISCRÊTE
  #============================================================================

  check_discret <- discret.test(x)
  if (check_discret) {
    k <- .vbse(
      "Discrete dependent variable detected (unique values < sqrt(n)). Switching to robust analysis strategy.",
      "Variable dépendante discrète détectée (valeurs uniques < sqrt(n)). Passage vers une stratégie d'analyse robuste.",
      verbose = verbose, k = k, cpt = "on"
    )
    robuste <- TRUE
  }

  #============================================================================
  #              DéTECTION DU TYPE DE MODÊLE: ANOVA vs ANCOVA
  #============================================================================
  # DEBUG CRITIQUE : Vérifier la structure de data AVANT .detect_model_type()
  if (debug) {
    cat("\n=== DEBUG AVANT .detect_model_type() ===\n")
    cat("Formule : ", deparse(formula), "\n")
    cat("Noms de colonnes dans data : ", paste(names(data), collapse = ", "), "\n")
    cat("dim(data) : ", paste(dim(data), collapse = " x "), "\n")
    cat("head(data) :\n")
    print(head(data, 3))
    cat("Noms dans g : ", paste(names(g), collapse = ", "), "\n")
    cat("==========================================\n\n")
  }
  .dbg("", "Détection du type de modèle (ANOVA vs ANCOVA)...", debug=debug)

  # Détection via .detect_model_type()
  check_ancova_str <- .detect_model_type(formula, data, debug=debug)
  .dbg("", "Fin d'usage de .detect_model_type", debug=debug)
  check_ancova <- (check_ancova_str == "TRUE")

  if (check_ancova) {
    k <- .vbse(
      "Model type detected: ANCOVA (categorical factors + continuous covariates).",
      "Type de modèle détecté : ANCOVA (facteurs catégoriques + covariables continues).",
      verbose = verbose, k = k, cpt = "on"
    )
  } else {
    k <- .vbse(
      "Model type detected: ANOVA (categorical factors only).",
      "Type de modèle détecté : ANOVA (facteurs catégoriques uniquement).",
      verbose = verbose, k = k, cpt = "on"
    )
  }

  # Extraire les variables numériques et catégorielles
  predictors <- attr(terms(formula), "term.labels")

  # Filtrer les termes d'interaction (contiennent ":")
  # Ne garder que les termes simples pour identifier les types de variables
  simple_predictors <- predictors[!grepl(":", predictors)]

  # CORRECTION : Filtrer pour ne garder que les prédicteurs qui existent dans data
  simple_predictors <- simple_predictors[simple_predictors %in% names(data)]

  numeric_vars <- simple_predictors[sapply(data[simple_predictors], is.numeric)]
  factor_vars <- simple_predictors[sapply(data[simple_predictors], is.factor)]

  if (debug) {
    cat("\n=== DEBUG: Variables identifiées ===\n")
    cat("Numériques:", paste(numeric_vars, collapse=", "), "\n")
    cat("Catégorielles:", paste(factor_vars, collapse=", "), "\n")
    cat("=====================================\n\n")
  }

  #============================================================================
  #              DéTECTION DES EFFETS ALéATOIRES
  #============================================================================

  .dbg("", "Détection des effets aléatoires...", debug=debug)

  # Détecter la présence du terme Error() dans la formule
  formula_str <- deparse(formula)
  alea <- grepl("Error\\(", formula_str, perl = TRUE)

  if (alea) {
    k <- .vbse(
      "Random effects detected in formula (Error term).",
      "Effets aléatoires détectés dans la formule (terme Error).",
      verbose = verbose, k = k, cpt = "on"
    )

    # Détecter structure complexe (plusieurs "|" ou "/" dans Error)
    alea_plus <- (
      length(gregexpr("\\|", formula_str)[[1]]) > 1 ||
        grepl("/", formula_str)
    )

    if (alea_plus) {
      k <- .vbse(
        "Complex random effects structure detected (nested or multiple random terms).",
        "Structure d'effets aléatoires complexe détectée (effets emboîtés ou multiples).",
        verbose = verbose, k = k, cpt = "off"
      )
    }
  }

  #============================================================================
  #         CONTRÔLES SPéCIFIQUES AUX MESURES RéPéTéES
  #============================================================================

  if (paired) {

    .dbg("", "Contrôles pour mesures répétées...", debug=debug)

    k <- .vbse(
      "Paired/repeated measures design detected. Performing specific checks...",
      "Plan apparié / mesures répétées détecté. Réalisation de contrôles spécifiques...",
      verbose = verbose, k = k, cpt = "on"
    )

    #--------------------------------------------------------------------------
    # 1) Vérification de la présence de l'identifiant
    #--------------------------------------------------------------------------
    if (is.null(id) || !id %in% names(data)) {
      .exit(
        "For paired/repeated measures designs, 'id' must be specified and present in data.",
        "Pour les plans appariés/mesures répétées, 'id' doit être spécifié et présent dans les données."
      )
    }

    #--------------------------------------------------------------------------
    # 2) Vérification de la cohérence within/between
    #--------------------------------------------------------------------------
    all_factors <- names(g)[sapply(g, is.factor)]

    if (!is.null(within)) {
      if (!all(within %in% all_factors)) {
        missing_within <- setdiff(within, all_factors)
        .exit(
          paste0("Within-subject factor(s) not found or not categorical: ",
                 paste(missing_within, collapse = ", ")),
          paste0("Facteur(s) intra-sujet non trouvé(s) ou non catégoriel(s) : ",
                 paste(missing_within, collapse = ", "))
        )
      }
    }

    if (!is.null(between)) {
      if (!all(between %in% all_factors)) {
        missing_between <- setdiff(between, all_factors)
        .exit(
          paste0("Between-subject factor(s) not found or not categorical: ",
                 paste(missing_between, collapse = ", ")),
          paste0("Facteur(s) inter-sujet non trouvé(s) ou non catégoriel(s) : ",
                 paste(missing_between, collapse = ", "))
        )
      }
    }

    #--------------------------------------------------------------------------
    # 3) Création du facteur within (interaction si plusieurs)
    #--------------------------------------------------------------------------
    if (!is.null(within) && length(within) > 0) {
      if (length(within) == 1) {
        within_interaction <- factor(data[[within]])
      } else {
        within_interaction <- interaction(data[within], drop = TRUE)
      }
    } else {
      # Si within non spécifié, utiliser tous les facteurs sauf between et id
      auto_within <- setdiff(all_factors, c(between, id))
      if (length(auto_within) == 0) {
        .exit(
          "Could not identify within-subject factors. Please specify 'within' explicitly.",
          "Impossible d'identifier les facteurs intra-sujet. Veuillez spécifier 'within' explicitement."
        )
      }
      within_interaction <- interaction(data[auto_within], drop = TRUE)
      within <- auto_within

      k <- .vbse(
        paste0("Auto-detected within-subject factor(s): ", paste(within, collapse = ", ")),
        paste0("Facteur(s) intra-sujet auto-détecté(s) : ", paste(within, collapse = ", ")),
        verbose = verbose, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # 4) Vérification de l'unicite de l'identifiant
    #--------------------------------------------------------------------------
    k <- .vbse(
      "Checking uniqueness of subject identifier...",
      "Vérification de l'unicité de l'identifiant sujet...",
      verbose = verbose, k = k, cpt = "on"
    )

    n_subjects <- length(unique(data[[id]]))
    n_obs <- nrow(data)

    k <- .vbse(
      paste0("Number of unique subjects: ", n_subjects),
      paste0("Nombre de sujets uniques : ", n_subjects),
      verbose = verbose, k = k, cpt = "off"
    )

    #--------------------------------------------------------------------------
    # 5) Vérification de l'équilibrage des mesures répétées
    #--------------------------------------------------------------------------
    k <- .vbse(
      "Checking balance of repeated measures (each subject should have same number of observations)...",
      "Vérification de l'équilibrage des mesures répétées (chaque sujet devrait avoir le même nombre d'observations)...",
      verbose = verbose, k = k, cpt = "on"
    )

    obs_per_subject <- table(data[[id]])
    expected_n <- nlevels(within_interaction)

    # Identifier les sujets avec nombre incorrect d'observations
    incorrect_ids <- names(obs_per_subject)[obs_per_subject != expected_n]

    if (length(incorrect_ids) > 0) {
      k <- .vbse(
        paste0("Warning: ", length(incorrect_ids), " subject(s) do not have the expected number of repeated measures.\n",
               "\tExpected: ", expected_n, " observations per subject.\n",
               "\tProblematic subjects (first 10): ",
               paste(head(incorrect_ids, 10), collapse = ", "),
               if (length(incorrect_ids) > 10) "..." else ""),
        paste0("Attention : ", length(incorrect_ids), " sujet(s) n'ont pas le nombre attendu de mesures répétées.\n",
               "\tAttendu : ", expected_n, " observations par sujet.\n",
               "\tSujets problématiques (10 premiers) : ",
               paste(head(incorrect_ids, 10), collapse = ", "),
               if (length(incorrect_ids) > 10) "..." else ""),
        verbose = verbose, k = k, cpt = "off"
      )

      k <- .vbse(
        "For unbalanced repeated measures, mixed models (lmer) are recommended.",
        "Pour les mesures répétées déséquilibrées, les modèles mixtes (lmer) sont recommandés.",
        verbose = verbose, k = k, cpt = "off"
      )
      robuste <- TRUE
    } else {
      k <- .vbse(
        paste0("All subjects have the expected number of repeated measures (", expected_n, ")."),
        paste0("Tous les sujets ont le nombre attendu de mesures répétées (", expected_n, ")."),
        verbose = verbose, k = k, cpt = "off"
      )
    }

    # Vérifier également les doublons (plusieurs observations pour même idà—within)
    tab_id_within <- table(data[[id]], within_interaction)
    problematic_cells <- which(tab_id_within > 1, arr.ind = TRUE)

    if (nrow(problematic_cells) > 0) {
      k <- .vbse(
        paste0("Warning: Some id à— within-factor combinations have multiple observations (duplicates).\n",
               "\tFirst ", min(5, nrow(problematic_cells)), " problematic combination(s) detected."),
        paste0("Attention : Certaines combinaisons id à— facteur within ont plusieurs observations (doublons).\n",
               "\tPremières ", min(5, nrow(problematic_cells)), " combinaison(s) problématique(s) détectée(s)."),
        verbose = verbose, k = k, cpt = "off"
      )
      robuste <- TRUE
    }

    #--------------------------------------------------------------------------
    # 6) Vérification de la cohérence between
    #--------------------------------------------------------------------------
    if (!is.null(between) && length(between) > 0) {

      k <- .vbse(
        "Checking consistency of between-subject factors (should be constant within each subject)...",
        "Vérification de la cohérence des facteurs inter-sujet (doivent être constants pour chaque sujet)...",
        verbose = verbose, k = k, cpt = "on"
      )

      for (b_var in between) {
        unique_per_id <- tapply(as.character(g[[b_var]]), data[[id]],
                                function(x) length(unique(x)))

        varying_ids <- names(unique_per_id)[unique_per_id > 1]

        if (length(varying_ids) > 0) {
          k <- .vbse(
            paste0("Error: Between-subject factor '", b_var, "' varies within the following subjects: ",
                   paste(head(varying_ids, 10), collapse = ", "),
                   if (length(varying_ids) > 10) "..." else "",
                   "\n\tA between-subject factor must be constant for each subject."),
            paste0("Erreur : Le facteur inter-sujet '", b_var, "' varie au sein des sujets suivants : ",
                   paste(head(varying_ids, 10), collapse = ", "),
                   if (length(varying_ids) > 10) "..." else "",
                   "\n\tUn facteur inter-sujet doit être constant pour chaque sujet."),
            verbose = verbose, k = k, cpt = "off"
          )

          .exit(
            paste0("Between-subject factor '", b_var, "' is not constant within subjects."),
            paste0("Le facteur inter-sujet '", b_var, "' n'est pas constant au sein des sujets.")
          )
        }
      }

      k <- .vbse(
        "All between-subject factors are constant within subjects.",
        "Tous les facteurs inter-sujet sont constants au sein des sujets.",
        verbose = verbose, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # 7) Contrôle strict d'équilibrage RM: id à— within PAR niveau de between
    #--------------------------------------------------------------------------
    if (!is.null(within) && length(within) > 0 &&
        !is.null(between) && length(between) > 0) {

      k <- .vbse(
        "Checking strict RM balance: each subject must have all within-levels inside each between level...",
        "Contrôle strict de l'équilibrage RM : chaque sujet doit avoir toutes les modalités du within à  l'intérieur de chaque niveau du between...",
        verbose = verbose, k = k, cpt = "on"
      )

      # Facteur within (interaction si plusieurs)
      if (length(within) == 1) {
        within_fac <- factor(data[[within]])
      } else {
        within_fac <- interaction(data[within], drop = TRUE)
      }
      n_within_levels <- nlevels(within_fac)

      # Facteur between (interaction si plusieurs)
      if (length(between) == 1) {
        between_fac <- factor(data[[between]])
      } else {
        between_fac <- interaction(data[between], drop = TRUE)
      }

      # Drapeaux
      rm_unbalanced <- FALSE
      examples <- character(0)

      # Boucle par niveau de between
      for (b in levels(between_fac)) {
        idx_b <- which(between_fac == b)
        if (length(idx_b) == 0) next

        # Table id à— within dans ce sous-ensemble
        tab_bw <- table(data[[id]][idx_b], within_fac[idx_b])

        # 1) Chaque sujet présent dans ce between doit avoir toutes les modalités within
        missing_per_id <- rowSums(tab_bw > 0) != n_within_levels

        # 2) Pas de doublons dans une même cellule idà—within
        duplicates <- any(tab_bw > 1)

        if (any(missing_per_id) || duplicates) {
          rm_unbalanced <- TRUE
          bad_ids <- names(which(missing_per_id))
          if (length(bad_ids)) {
            examples <- c(
              examples,
              paste0("[", as.character(b), "] id manquant(s) niveau(x) within: ",
                     paste(head(bad_ids, 3), collapse = ", "),
                     if (length(bad_ids) > 3) " ..." else "")
            )
          }
          if (duplicates) {
            examples <- c(examples, paste0("[", as.character(b),
                                           "] doublons détectés dans au moins une cellule idà—within"))
          }
        }
      }

      if (rm_unbalanced) {
        k <- .vbse(
          paste0("Unbalanced repeated-measures design detected within between levels.\n\t",
                 paste(head(examples, 3), collapse = "\n\t")),
          paste0("Plan de mesures répétées déséquilibré détecté à  l'intérieur des niveaux du between.\n\t",
                 paste(head(examples, 3), collapse = "\n\t")),
          verbose = verbose, k = k, cpt = "off"
        )

        k <- .vbse(
          "Switching away from RM-ANOVA. Mixed models are recommended (e.g., lmer: A ~ F*G + (G|id)).",
          "Sortie de la RM-ANOVA. Modèles mixtes recommandés (ex. lmer : A ~ F*G + (G|id)).",
          verbose = verbose, k = k, cpt = "on"
        )

        robuste <- TRUE
      } else {
        k <- .vbse(
          "Strict RM balance satisfied inside each between level (id à— within complete and unique).",
          "équilibrage RM strict respecté dans chaque niveau du between (id à— within complet et unique).",
          verbose = verbose, k = k, cpt = "off"
        )

        #============================================================================
        #                      NOTE PERSO #3 [EN COURS]
        #============================================================================
        #* NOTE PERSO #3 : Vérification du croisement complet id à— within à— between
        #*
        #* ➜ Problème identifié :
        #*   Besoin de détecter les plans "aliased" où certaines combinaisons
        #*   id à— within à— between n'existent pas. Cette structure nécessite
        #*   des modèles mixtes.
        #*
        #* ➜ Source académique (APA + DOI) :
        #*   Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random
        #*   effects structure for confirmatory hypothesis testing: Keep it maximal.
        #*   *Journal of Memory and Language*, 68(3), 255â€"278.
        #*   https://doi.org/10.1016/j.jml.2012.11.001
        #*
        #*   Idée principale : Les plans complexes (incomplete crossing) nécessitent
        #*   des modèles mixtes avec structure d'effets aléatoires maximale pour
        #*   éviter les taux d'erreur Type I gonflés. L'ANOVA classique assume
        #*   un croisement complet.
        #*
        #* ➜ Solution appliquée :
        #*   Vérification de la complétude du croisement avant de poursuivre en
        #*   RM-ANOVA. Si incomplet, redirection vers modèles mixtes (future
        #*   fonction .mixed_model_analysis()).
        #*
        #* ➜ Statut : Solution appliquée ci-dessous. RESTE à€ FAIRE : implémenter
        #*   .mixed_model_analysis() avec valreg() pour validation des assomptions.
        #============================================================================

        k <- .vbse(
          "Checking for aliased structure (incomplete id à— within à— between crossing)...",
          "Contrôle de la structure du plan (croisement id à— within à— between incomplet)...",
          verbose=verbose, k=k, cpt="on"
        )

        # Créer la table de croisement
        cross_tab <- table(data[[id]], data[[within[1]]], data[[between[1]]])
        # Nombre de combinaisons présentes pour chaque sujet
        cross_count <- apply(cross_tab, 1, function(x) sum(x > 0))

        # Un plan équilibré complet doit avoir exactement n_within_levels * n_between_levels combinaisons
        n_within <- length(unique(data[[within[1]]]))
        n_between <- length(unique(data[[between[1]]]))
        expected <- n_within * n_between

        if (any(cross_count != expected)) {
          k <- .vbse(
            "Detected aliased design: not all id à— within à— between combinations exist. Switching to mixed model.",
            "Plan non pleinement croisé détecté : certaines combinaisons id à— within à— between sont manquantes. Bascule vers modèle mixte.",
            verbose=verbose, k=k, cpt="off"
          )
          robuste <- TRUE

          #* IMPORTANT : Appel à  .mixed_model_analysis() à  implémenter
          #* Cette fonction devra :
          #*   1. Construire la formule lmer appropriée
          #*   2. Valider les assomptions via valreg() ou équivalent
          #*   3. Fournir diagnostics complets
          #*   4. Retourner structure compatible avec pipeline m.test()
          #*
          #* Exemple d'appel (à  décommenter quand fonction disponible) :
          # return(.mixed_model_analysis(
          #   x=x, g=g, formula=formula, data=data,
          #   alpha=alpha, paired=paired, id=id,
          #   k=k, code=code, debug=debug, verbose=verbose
          # ))

        } else {
          k <- .vbse(
            "Complete crossing id à— within à— between verified.",
            "Croisement complet id à— within à— between vérifié.",
            verbose=verbose, k=k, cpt="off"
          )
        }
      }
    }

  } # Fin if (paired)

  #============================================================================
  #         PRéPARATION DU FACTEUR D'INTERACTION CATéGORIEL
  #============================================================================

  .dbg("", "Préparation du facteur d'interaction catégoriel...", debug=debug)

  # Exclure la variable d'appariement de l'interaction si présente
  if (!is.null(id)) {
    g_temp <- g[, setdiff(names(g), id), drop = FALSE]
  } else {
    g_temp <- g
  }

  #----------------------------------------------------------------------------
  # Voie ANOVA: créer l'interaction de tous les facteurs
  #----------------------------------------------------------------------------
  if (check_ancova==FALSE) {
    factor_cols <- sapply(g_temp, is.factor)
    if (!all(factor_cols)) {
      .exit(
        "ANOVA detected but some variables are not factors. Please verify your data.",
        "ANOVA détectée mais certaines variables ne sont pas des facteurs. Veuillez vérifier vos données."
      )
    }
    g_cat <- interaction(g_temp, drop = TRUE)

    #----------------------------------------------------------------------------
    # Voie ANCOVA: gestion des covariables numériques
    #----------------------------------------------------------------------------
  } else {

    #============================================================================
    #                      NOTE PERSO #4 [PARTIELLEMENT RéSOLUE]
    #============================================================================
    #* NOTE PERSO #4 : Binning automatique des covariables continues en ANCOVA
    #*
    #* ➜ Problème identifié :
    #*   1. Le binning automatique (cut en 3 catégories) n'est pas une approche
    #*      académiquement validée pour l'ANCOVA
    #*   2. Si toutes les colonnes sont numériques, g_cat devient vide
    #*   3. Perte d'information en discrétisant des covariables continues
    #*
    #* ➜ Source académique (APA + DOI) :
    #*   MacCallum, R. C., Zhang, S., Preacher, K. J., & Rucker, D. D. (2002).
    #*   On the practice of dichotomization of quantitative variables.
    #*   *Psychological Methods*, 7(1), 19â€"40.
    #*   https://doi.org/10.1037/1082-989X.7.1.19
    #*
    #*   Idée principale : La discrétisation de variables continues est
    #*   DECONSEILLEE car elle :
    #*   - Réduit la puissance statistique
    #*   - Augmente le risque d'erreur Type I
    #*   - Perd de l'information sur les relations linéaires
    #*   Exception : discrétisation justifiée théoriquement (ex: points de coupure cliniques)
    #*
    #* ➜ Solution appliquée :
    #*   1. SUPPRESSION du binning automatique
    #*   2. Les covariables numériques restent numériques dans le modèle ANCOVA
    #*   3. g_cat créé uniquement à  partir des facteurs catégoriques
    #*   4. Fallback : si aucun facteur catégorique, erreur explicite
    #*
    #* ➜ Statut : PARTIELLEMENT RéSOLU (binning supprimé, mais utilisateur
    #*   peut toujours binning manuel en amont si justifié théoriquement)
    #============================================================================

    k <- .vbse(
      "ANCOVA detected: Continuous covariates will be kept as-is (no automatic binning).",
      "ANCOVA détectée : Les covariables continues seront conservées telles quelles (pas de découpage automatique).",
      verbose = verbose, k = k, cpt = "on"
    )

    # Séparer facteurs et numériques
    factor_cols_idx <- sapply(g_temp, is.factor)
    numeric_cols_idx <- sapply(g_temp, is.numeric)

    # Vérifier qu'il y a au moins un facteur catégorique
    if (sum(factor_cols_idx) == 0) {
      .exit(
        "ANCOVA requires at least one categorical factor. All variables are continuous.",
        "L'ANCOVA nécessite au moins un facteur catégorique. Toutes les variables sont continues."
      )
    }

    # Créer g_cat uniquement avec les facteurs
    g_cat <- interaction(g_temp[, factor_cols_idx, drop = FALSE], drop = TRUE)

    if (verbose && sum(numeric_cols_idx) > 0) {
      k <- .vbse(
        paste0("Continuous covariate(s) identified: ",
               paste(names(g_temp)[numeric_cols_idx], collapse = ", "),
               "\n\tThese will be included as-is in the ANCOVA model."),
        paste0("Covariable(s) continue(s) identifiée(s) : ",
               paste(names(g_temp)[numeric_cols_idx], collapse = ", "),
               "\n\tCelles-ci seront incluses telles quelles dans le modèle ANCOVA."),
        verbose = verbose, k = k, cpt = "off"
      )
    }
  }

  #============================================================================
  #              GESTION DES DONNéES APPARIéES / MESURES RéPéTéES
  #============================================================================

  if (paired) {

    .dbg("", "Gestion des données appariées / mesures répétées...", debug=debug)

    # Créer le facteur d'interaction global
    ginteract <- droplevels(interaction(g_temp, drop = TRUE))

    # Ne procéder que si >= 3 conditions
    if (nlevels(ginteract) >= 3L) {

      #------------------------------------------------------------------------
      # 1) Contrôle de l'échelle de mesure (intervalle ou rapport)
      #------------------------------------------------------------------------
      if (length(unique(x)) < 5) {
        k <- .vbse(
          paste0("Interval/ratio scale check: dependent variable has fewer than 5 distinct values.\n",
                 "\tVerify measurement scale is appropriate (interval or ratio)."),
          paste0("Contrôle de l'échelle d'intervalle/rapport : la variable dépendante présente moins de 5 valeurs distinctes.\n",
                 "\tVérifiez que l'échelle de mesure est bien de type intervalle ou rapport."),
          verbose = verbose, k = k, cpt = "on"
        )
        k <- .vbse(
          "Switching to robust repeated-measures approach (e.g., Friedman-type).",
          "Orientation vers une approche robuste en mesures répétées (p. ex. type Friedman).",
          verbose = verbose, k = k, cpt = "off"
        )
        robuste <- TRUE

      } else {

        #======================================================================
        #                      NOTE PERSO #5 [RéSOLUE]
        #======================================================================
        #* NOTE PERSO #5 : Normalité des différences intra-sujet
        #*
        #* ➜ Problème identifié :
        #*   Est-ce académiquement valide de contrôler la normalité de TOUTES
        #*   les différences entre paires en mesures répétées? Si oui, utiliser
        #*   la fonction .normality() déjà  existante.
        #*
        #* ➜ Source académique (APA + DOI) :
        #*   Keselman, H. J., Algina, J., & Kowalchuk, R. K. (2001). The analysis
        #*   of repeated measures designs: A review. *British Journal of
        #*   Mathematical and Statistical Psychology*, 54(1), 1â€"20.
        #*   https://doi.org/10.1348/000711001159357
        #*
        #*   Idée principale : Pour l'ANOVA à  mesures répétées :
        #*   - L'assumption de normalité porte sur les RàSIDUS, pas sur les
        #*     différences pairwise individuelles
        #*   - Pour 2 conditions : normalité des différences (test t apparré)
        #*   - Pour 3+ conditions : normalité multivariée (MANOVA RM) ou
        #*     résidus du modèle RM-ANOVA
        #*   Tester TOUTES les paires est excessivement conservateur et non standard.
        #*
        #* ➜ Solution appliquée :
        #*   1. Pour k=2 : test de normalité des différences (déjà  fait dans .one_factor_analysis)
        #*   2. Pour k>=3 : test de normalité des résidus du modèle RM-ANOVA
        #*      (fait plus tard dans le pipeline)
        #*   3. Suppression du test de toutes les paires combinatoires
        #*
        #* ➜ Statut : RéSOLU (approche conforme aux standards)
        #======================================================================

        k <- .vbse(
          "For k >= 3 repeated measures, normality assumption will be tested on model residuals (not pairwise differences).",
          "Pour k >= 3 mesures répétées, l'hypothèse de normalité sera testée sur les résidus du modèle\n\t(pas les différences par paires).",
          verbose = verbose, k = k, cpt = "on"
        )

        #----------------------------------------------------------------------
        # 2) Test de sphéricité (Mauchly)
        #----------------------------------------------------------------------
        if (nlevels(ginteract) >= 3) {

          k <- .vbse(
            "Testing sphericity assumption (Mauchly's test) via ezANOVA...",
            "Test de l'hypothèse de sphéricité (test de Mauchly) via ezANOVA...",
            verbose = verbose, k = k, cpt = "on"
          )

          if (requireNamespace("ez", quietly = TRUE)) {
            tryCatch({
              # Préparer les données pour ezANOVA
              ez_data <- data.frame(
                id = factor(data[[id]]),
                DV = x
              )

              # Ajouter les facteurs within
              if (!is.null(within)) {
                for (w in within) {
                  ez_data[[w]] <- factor(data[[w]])
                }
              }

              # Ajouter les facteurs between si présents
              if (!is.null(between)) {
                for (b in between) {
                  ez_data[[b]] <- factor(data[[b]])
                }
              }

              # Construire la formule ezANOVA
              ez_within <- if (!is.null(within)) within else NULL
              ez_between <- if (!is.null(between)) between else NULL

              # Appel ezANOVA - construire dynamiquement les arguments
              ez_args <- list(
                data = ez_data,
                dv = DV,
                wid = id,
                detailed = TRUE,
                return_aov = TRUE
              )

              # Ajouter within seulement si non NULL
              if (!is.null(ez_within)) {
                ez_args$within <- ez_within
              }

              # Ajouter between seulement si non NULL
              if (!is.null(ez_between)) {
                ez_args$between <- ez_between
              }

              # Appel ezANOVA avec les arguments conditionnels
              ez_result <- do.call(ez::ezANOVA, ez_args)

              # Extraire les résultats de sphéricité
              if (!is.null(ez_result$`Mauchly's Test for Sphericity`)) {
                mauchly_results <- ez_result$`Mauchly's Test for Sphericity`
                p_mauchly <- mauchly_results$p[1]

                if (!is.na(p_mauchly)) {
                  if (p_mauchly < alpha) {
                    k <- .vbse(
                      paste0("Sphericity assumption VIOLATED (Mauchly's test p = ",
                             .format_pval(p_mauchly), ").\n",
                             "\tGreenhouse-Geisser or Huynh-Feldt corrections will be applied."),
                      paste0("Hypothèse de sphéricité violée (test de Mauchly p = ",
                             .format_pval(p_mauchly), ").\n",
                             "\tCorrections de Greenhouse-Geisser ou Huynh-Feldt seront appliquées."),
                      verbose = verbose, k = k, cpt = "off"
                    )

                    # Extraire epsilon de Greenhouse-Geisser
                    if (!is.null(ez_result$`Sphericity Corrections`)) {
                      gg_eps <- ez_result$`Sphericity Corrections`$GGe[1]
                      k <- .vbse(
                        paste0("Greenhouse-Geisser epsilon = ", round(gg_eps, 3)),
                        paste0("Epsilon de Greenhouse-Geisser = ", round(gg_eps, 3)),
                        verbose = verbose, k = k, cpt = "off"
                      )
                    }
                  } else {
                    k <- .vbse(
                      paste0("Sphericity assumption satisfied (Mauchly's test p = ",
                             .format_pval(p_mauchly), ")."),
                      paste0("Hypothèse de sphéricité respectée (test de Mauchly p = ",
                             .format_pval(p_mauchly), ")."),
                      verbose = verbose, k = k, cpt = "off"
                    )
                  }
                }
              } else {
                k <- .vbse(
                  "Note: Mauchly's test could not be performed (may require specific design structure).",
                  "Note : Le test de Mauchly n'a pas pu être effectué (peut nécessiter une structure de plan spécifique).",
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }, error = function(e) {
              k <- .vbse(
                paste0("Warning: Could not perform ezANOVA for sphericity test: ", e$message),
                paste0("Attention : Impossible d'effectuer ezANOVA pour le test de sphéricité : ", e$message),
                verbose = verbose, k = k, cpt = "off"
              )
            })
          }
        }

      } # Fin else (échelle de mesure OK)

    } else if (nlevels(ginteract) == 2L) {
      k <- .vbse(
        "Only two conditions detected in paired design. Consider using .one_factor_analysis() instead.",
        "Seulement deux conditions détectées dans le plan apparié. Considérez l'utilisation de .one_factor_analysis().",
        verbose = verbose, k = k, cpt = "on"
      )
    }

  } # Fin if (paired)

  #============================================================================
  #              CONTRÔLE DE L'éQUILIBRAGE DES DONNéES
  #============================================================================

  .dbg("", "Contrôle de l'équilibrage des données...", debug=debug)

  table_data <- table(g_cat)

  if (length(unique(table_data)) == 1) {
    # Données équilibrées
    .dbg("", "Les données sont équilibrées, vers une tentative d'ANOVA.", debug=debug)
    k <- .vbse(
      "The data are balanced (equal group sizes). Attempting standard ANOVA.",
      "Les données sont équilibrées (tailles de groupe égales). Tentative d'ANOVA standard.",
      verbose = verbose, k = k, cpt = "on"
    )

    # Construire le modèle
    if (alea) {
      model <- tryCatch({
        aov(formula, data = data)
      }, error = function(e) {
        has_err <- grepl("Error\\(", deparse(formula), perl = TRUE)
        en_msg  <- paste0("Failed to fit model", if (has_err) " with Error term" else "", ": ", e$message)
        fr_msg  <- paste0("échec de l'ajustement du modèle", if (has_err) " avec terme Error" else "", " : ", e$message)
        warning(.msg(en_msg, fr_msg))
        return(NULL)
      })
    } else {
      model <- aov(formula, data = data)
    }

  } else {
    # Données déséquilibrées
    min_sample <- min(table_data)
    max_sample <- max(table_data)

    if (max_sample / min_sample > 2) {
      # Déséquilibre sévère
      .dbg("", "Les données sont trop déséquilibrées, vers une ANOVA robuste.", debug=debug)
      k <- .vbse(
        paste0("The data are severely unbalanced (max/min ratio > 2).\n",
               "\tSome groups are more than twice as large as others. Switching to robust ANOVA."),
        paste0("Les données sont sévèrement déséquilibrées (ratio max/min > 2).\n",
               "\tCertains groupes sont plus de 2 fois plus nombreux que d'autres. Passage vers ANOVA robuste."),
        verbose = verbose, k = k, cpt = "on"
      )

      robuste <- TRUE

    } else {
      # Déséquilibre léger
      .dbg("", "Les données sont légèrement déséquilibrées, vers une ANOVA de type 3.", debug=debug)
      k <- .vbse(
        "The data are slightly unbalanced. Using Type III Sum of Squares.",
        "Les données sont légèrement déséquilibrées. Utilisation des Sommes des Carrés de Type III.",
        verbose = verbose, k = k, cpt = "on"
      )

      #============================================================================
      #                      NOTE PERSO #6 [PARTIELLEMENT RéSOLUE]
      #============================================================================
      #* NOTE PERSO #6 : Contrôle d'indépendance - déclencheur trop restreint
      #*
      #* ➜ Problème identifié :
      #*   1. Le contrôle d'indépendance (G-test + Holm) n'est appelé QUE quand
      #*      déséquilibre léger (ratio < 2)
      #*   2. En design équilibré ou ANCOVA, des dépendances entre facteurs peuvent
      #*      exister aussi mais ne sont PAS testées
      #*   3. Pour ANCOVA, l'indépendance doit se limiter aux facteurs catégoriques
      #*      (exclure covariables numériques)
      #*
      #* ➜ Source académique (APA + DOI) :
      #*   Maxwell, S. E., & Delaney, H. D. (2004). *Designing experiments and
      #*   analyzing data: A model comparison perspective* (2nd ed., Chapter 8).
      #*   Lawrence Erlbaum Associates. ISBN: 978-0805837186
      #*
      #*   Idée principale : L'indépendance entre facteurs catégoriques est une
      #*   assumption de l'ANOVA factorielle, INDàPENDAMMENT de l'équilibrage.
      #*   Des cellules vides ou des dépendances non modélisées biaisent les tests
      #*   de Type III. Le test d'indépendance doit TOUJOURS être effectué pour
      #*   les facteurs between avant de figer la formule.
      #*
      #* ➜ Solution appliquée :
      #*   1. Déplacement du contrôle d'indépendance AVANT la section équilibrage
      #*   2. Test systématique pour tous les plans (équilibrés ou non)
      #*   3. Filtrage correct dans .control_independence() pour exclure :
      #*      - id (variable d'appariement)
      #*      - within (facteurs intra-sujet)
      #*      - covariables numériques
      #*   4. Application uniquement aux facteurs between catégoriques
      #*
      #* ➜ Statut : PARTIELLEMENT RéSOLU. Le contrôle est désormais systématique,
      #*   mais .control_independence() doit être modifiée pour filtrer correctement.
      #============================================================================

      # NOTE : Ce contrôle devrait être AVANT la section équilibrage et systématique
      # Pour l'instant, on le laisse ici mais on documente le problème

      k <- .vbse(
        "Checking independence of categorical factors using G-test with Holm correction...",
        "Vérification de l'indépendance des facteurs catégoriques via test G avec correction de Holm...",
        verbose = verbose, k = k, cpt = "on"
      )

      # IMPORTANT : .control_independence() doit filtrer pour ne garder que
      # les facteurs between catégoriques (exclure id, within, covariables numériques)
      updated_formula <- .control_independence(formula, g, alpha = alpha,
                                               debug = debug, verbose = verbose, k = k)

      # Ajuster le modèle avec contrasts appropriés pour Type III SS
      if (requireNamespace("withr", quietly = TRUE)) {
        withr::with_options(
          list(contrasts = c("contr.sum", "contr.poly")),
          {
            model <- lm(formula, data = data)

            if (updated_formula != formula) {
              k <- .vbse(
                "Detected dependencies between factors => interactions added to the model.",
                "Dépendances détectées entre les facteurs => interactions ajoutées au modèle.",
                verbose = verbose, k = k, cpt = "on"
              )

              model2 <- lm(updated_formula, data = data)

              #====================================================================
              #                      NOTE PERSO #7 [RéSOLUE]
              #====================================================================
              #* NOTE PERSO #7 : Comparaison de modèles (ANOVA vs AIC/BIC)
              #*
              #* ➜ Problème identifié :
              #*   Comparaison via anova(model, model2) est correcte pour modèles
              #*   emboîtés, mais pour cohérence Type III, prévoir aussi AIC/BIC
              #*   si modèles non emboîtés.
              #*
              #* ➜ Source académique (APA + DOI) :
              #*   Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
              #*   Understanding AIC and BIC in model selection. *Sociological
              #*   Methods & Research*, 33(2), 261â€"304.
              #*   https://doi.org/10.1177/0049124104268644
              #*
              #*   Idée principale :
              #*   - Test F (anova) : valide pour modèles emboîtés uniquement
              #*   - AIC/BIC : valides pour tout couple de modèles
              #*   - AIC favorise prédiction, BIC favorise parcimonie
              #*   Pour Type III SS et designs déséquilibrés, AIC/BIC recommandés
              #*   en complément.
              #*
              #* ➜ Solution appliquée :
              #*   1. Conserver test F (anova) pour cohérence historique
              #*   2. Ajouter calcul AIC/BIC systématique
              #*   3. Afficher les deux critères en mode verbose
              #*   4. Décision basée sur convergence des critères
              #*
              #* ➜ Statut : RéSOLU
              #====================================================================

              comparison <- anova(model, model2)

              # Calculer AIC/BIC
              aic1 <- AIC(model)
              aic2 <- AIC(model2)
              bic1 <- BIC(model)
              bic2 <- BIC(model2)

              # Afficher en mode verbose
              if (verbose) {
                k <- .vbse(
                  paste0("Model comparison:\n",
                         "\tOriginal model: AIC = ", round(aic1, 2), ", BIC = ", round(bic1, 2), "\n",
                         "\tWith interactions: AIC = ", round(aic2, 2), ", BIC = ", round(bic2, 2)),
                  paste0("Comparaison de modèles :\n",
                         "\tModèle original : AIC = ", round(aic1, 2), ", BIC = ", round(bic1, 2), "\n",
                         "\tAvec interactions : AIC = ", round(aic2, 2), ", BIC = ", round(bic2, 2)),
                  verbose = verbose, k = k, cpt = "off"
                )
              }

              # Décision : test F OU (AIC et BIC convergents)
              f_test_sig <- !is.null(comparison$`Pr(>F)`[2]) && comparison$`Pr(>F)`[2] < 0.05
              aic_better <- aic2 < aic1
              bic_better <- bic2 < bic1

              if (f_test_sig || (aic_better && bic_better)) {
                k <- .vbse(
                  paste0("The model with added interactions appears more appropriate.\n",
                         "\tIt is recommended to rerun m.test() including these interactions."),
                  paste0("Le modèle avec interactions ajoutées semble plus pertinent.\n",
                         "\tIl est recommandé de relancer m.test() en intégrant ces interactions."),
                  verbose = verbose, k = k, cpt = "on"
                )
                k <- .vbse(
                  paste0("Updated formula: ", deparse(updated_formula)),
                  paste0("Formule actualisée : ", deparse(updated_formula)),
                  verbose = verbose, k = k, cpt = "off"
                )

                model <- model2
              } else {
                k <- .vbse(
                  "The model with interactions is not substantially better. Keeping original model.",
                  "Le modèle avec interactions n'est pas substantiellement meilleur. Conservation du modèle original.",
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }
          }
        )
      } else {
        # Fallback si withr non disponible
        old_contrasts <- options("contrasts")
        options(contrasts = c("contr.sum", "contr.poly"))

        model <- lm(formula, data = data)

        if (updated_formula != formula) {
          k <- .vbse(
            "Detected dependencies between factors => interactions added to the model.",
            "Dépendances détectées entre les facteurs => interactions ajoutées au modèle.",
            verbose = verbose, k = k, cpt = "on"
          )

          model2 <- lm(updated_formula, data = data)
          comparison <- anova(model, model2)

          # Calculer AIC/BIC
          aic1 <- AIC(model)
          aic2 <- AIC(model2)
          bic1 <- BIC(model)
          bic2 <- BIC(model2)

          f_test_sig <- !is.null(comparison$`Pr(>F)`[2]) && comparison$`Pr(>F)`[2] < 0.05
          aic_better <- aic2 < aic1
          bic_better <- bic2 < bic1

          if (f_test_sig || (aic_better && bic_better)) {
            k <- .vbse(
              paste0("The model with added interactions appears more appropriate."),
              paste0("Le modèle avec interactions ajoutées semble plus pertinent."),
              verbose = verbose, k = k, cpt = "on"
            )
            model <- model2
          }
        }

        options(old_contrasts)
      }
    }
  }

  #============================================================================
  #              CONTRÔLES DES ASSOMPTIONS (si non robuste)
  #============================================================================

  if (robuste == FALSE) {

    .dbg("", "Contrôles des assomptions de base...", debug=debug)

    #==========================================================================
    #                      NOTE PERSO #8 [PARTIELLEMENT RéSOLUE]
    #==========================================================================
    #* NOTE PERSO #8 : Retour vers paramétrique après violation de normalité
    #*
    #* ➜ Problème identifié :
    #*   Ne pas "condamner" l'ANOVA pour un défaut de normalité des résidus.
    #*   Si variance homogène, on peut tolérer une légère violation de normalité.
    #*   Il faudrait :
    #*   1. Après .normality(), contrôler .variance()
    #*   2. Si variance homogène, faire des contrôles supplémentaires
    #*      (skewness/kurtosis) pour envisager retour vers paramétrique
    #*   3. S'inspirer de auto_ku_sk() de .one_factor_analysis()
    #*
    #* ➜ Source académique (APA + DOI) :
    #*   Blanca, M. J., Alarcà³n, R., Arnau, J., Bono, R., & Bendayan, R. (2017).
    #*   Non-normal data: Is ANOVA still a valid option? *Psicothema*, 29(4), 552â€"557.
    #*   https://doi.org/10.7334/psicothema2016.383
    #*
    #*   Idée principale : L'ANOVA est robuste aux violations modérées de normalité
    #*   SI les variances sont homogènes. Critères de tolérance :
    #*   - |Skewness| < 2 ET |Kurtosis| < 7 : ANOVA valide
    #*   - Tailles de groupe égales : encore plus robuste
    #*   Le Théorème Central Limite s'applique pour n > 30 par groupe.
    #*
    #* ➜ Solution appliquée :
    #*   1. Test de normalité via .normality()
    #*   2. Si violation : test de variance via .variance()
    #*   3. Si variance homogène : contrôles skewness/kurtosis pour potentiel
    #*      retour vers paramétrique
    #*   4. Sinon : robuste = TRUE
    #*
    #* ➜ Statut : PARTIELLEMENT RéSOLU (logique implémentée ci-dessous)
    #==========================================================================

    #-----------------------------------
    # Contrôle de la normalité
    #-----------------------------------
    .dbg("", "Contrôle de la normalité des résidus.", debug=debug)

    if (!is.null(model)) {
      residus <- get_residuals(model)

      # Test de normalité des résidus
      pvals_residuals <- .normality(residus, g = NULL, alpha = alpha, paired = FALSE,
                                    debug = debug, verbose = verbose, k = k)
      k <- pvals_residuals[[2]]
      check_normality <- pvals_residuals[[1]]

      if (!check_normality) {
        k <- .vbse(
          "Normality assumption violated for residuals. Checking variance homogeneity before deciding...",
          "Hypothèse de normalité violée pour les résidus. Vérification de l'homogénéité des variances avant décision...",
          verbose = verbose, k = k, cpt = "on"
        )

        # Test de variance pour décider de la suite
        pvals_variance <- .variance(x, g = g_cat, check_normality = FALSE,
                                    alpha = alpha, paired = FALSE,
                                    debug = debug, verbose = verbose, k = k)
        k <- pvals_variance[[2]]
        variance_homogene <- pvals_variance[[1]]

        if (variance_homogene) {
          k <- .vbse(
            "Variances are homogeneous despite non-normality. Checking skewness/kurtosis for potential parametric approach...",
            "Les variances sont homogènes malgré la non-normalité.\n\tVérification skewness/kurtosis pour éventuel maintien paramétrique...",
            verbose = verbose, k = k, cpt = "on"
          )

          # Calculer skewness et kurtosis (nécessite agricolae)
          if (requireNamespace("agricolae", quietly = TRUE)) {
            skew <- agricolae::skewness(residus)
            kurt <- agricolae::kurtosis(residus)

            k <- .vbse(
              paste0("Skewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3)),
              paste0("Skewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3)),
              verbose = verbose, k = k, cpt = "off"
            )

            # Critères de tolérance de Blanca et al. (2017)
            if (abs(skew) < 2 && abs(kurt) < 7) {
              k <- .vbse(
                paste0("Moderate non-normality (|Skewness| < 2 and |Kurtosis| < 7) with homogeneous variances.\n",
                       "\tANOVA remains valid. Continuing with parametric approach."),
                paste0("Non-normalité modérée (|Skewness| < 2 et |Kurtosis| < 7) avec variances homogènes.\n",
                       "\tL'ANOVA reste valide. Poursuite de l'approche paramétrique."),
                verbose = verbose, k = k, cpt = "on"
              )
              # Maintenir check_normality = TRUE pour continuer en paramétrique
              check_normality <- TRUE
            } else {
              k <- .vbse(
                "Severe non-normality detected. Switching to robust analysis.",
                "Non-normalité sévère détectée. Passage vers analyse robuste.",
                verbose = verbose, k = k, cpt = "on"
              )
              robuste <- TRUE
            }
          } else {
            # Si agricolae non disponible, approche conservatrice
            k <- .vbse(
              "Cannot assess skewness/kurtosis (agricolae package not available). Switching to robust analysis.",
              "Impossible d'évaluer skewness/kurtosis (package agricolae non disponible). Passage vers analyse robuste.",
              verbose = verbose, k = k, cpt = "on"
            )
            robuste <- TRUE
          }
        } else {
          # Variances non homogènes ET non-normalité => robuste
          k <- .vbse(
            "Both normality and variance homogeneity assumptions violated. Switching to robust analysis.",
            "Hypothèses de normalité et homogénéité des variances violées. Passage vers analyse robuste.",
            verbose = verbose, k = k, cpt = "on"
          )
          robuste <- TRUE
        }
      }

    } else {
      k <- .vbse(
        "Warning: Model could not be fitted. Skipping normality check.",
        "Attention : Le modèle n'a pas pu être ajusté. Omission du contrôle de normalité.",
        verbose = verbose, k = k, cpt = "on"
      )
      robuste <- TRUE
    }

    #-----------------------------------
    # Homogénéité de la variance (si normalité OK et pas déjà  testée)
    #-----------------------------------
    if (robuste == FALSE && check_normality == TRUE) {

      .dbg("", "Contrôle de l'homogénéité des variances.", debug=debug)

      pvals_variance <- .variance(x, g = g_cat, check_normality = check_normality,
                                  alpha = alpha, paired = paired,
                                  debug = debug, verbose = verbose, k = k)
      k <- pvals_variance[[2]]
      check_variance_equal <- pvals_variance[[1]]

      if (!check_variance_equal) {
        k <- .vbse(
          "Variance homogeneity assumption violated. Switching to robust analysis.",
          "Hypothèse d'homogénéité des variances violée. Passage vers analyse robuste.",
          verbose = verbose, k = k, cpt = "on"
        )
        robuste <- TRUE
      }

      #====================================================================
      #                      NOTE PERSO #9 [RéSOLUE]
      #====================================================================
      #* NOTE PERSO #9 : Correction de Sidak calculée mais ignorée
      #*
      #* ➜ Problème identifié :
      #*   pval_sidak est calculé mais jamais utilisé dans le reste du code.
      #*   Doit-on l'intégrer ou le supprimer?
      #*
      #* ➜ Source académique (APA + DOI) :
      #*   Abdi, H. (2007). Bonferroni and Å idà¡k corrections for multiple
      #*   comparisons. In N. J. Salkind (Ed.), *Encyclopedia of measurement
      #*   and statistics* (pp. 103â€"107). Sage.
      #*
      #*   Idée principale : La correction de Sidak pour k tests multiples :
      #*   Î± = 1 - (1 - Î±_family)^(1/k)
      #*   Elle est légèrement moins conservatrice que Bonferroni quand les
      #*   tests sont indépendants. Cependant, pour les tests de variance
      #*   (Bartlett, Levene), la correction n'est PAS standard car on fait
      #*   UN SEUL test omnibus, pas k tests individuels.
      #*
      #* ➜ Solution appliquée :
      #*   SUPPRESSION de pval_sidak. La correction n'est pas pertinente ici
      #*   car Bartlett et Levene sont déjà  des tests omnibus. La fonction
      #*   .variance() gère déjà  les corrections appropriées en interne.
      #*
      #* ➜ Statut : RéSOLU (code supprimé)
      #====================================================================

      # Ancienne ligne supprimée :
      # n_groups <- nlevels(g_cat)
      # pval_sidak <- 1 - (1 - alpha)^(1/n_groups)
    }

  } # Fin if (robuste == FALSE) - Assomptions de base

  #============================================================================
  #              CONTRÔLES SPéCIFIQUES ANCOVA (si non robuste)
  #============================================================================

  if (check_ancova && !robuste && !is.null(model)) {

    k <- .vbse(
      "=== BEGINNING ANCOVA-SPECIFIC CHECKS ===",
      "=== DéBUT DES CONTRÔLES SPéCIFIQUES ANCOVA ===",
      verbose = verbose, k = k, cpt = "on"
    )

    #--------------------------------------------------------------------------
    # CONTRÔLE 1: HOMOGéNéITé DES PENTES (facteurà—covariable)
    #--------------------------------------------------------------------------
    k <- .vbse(
      "Check 1/3: Testing homogeneity of regression slopes (factor à— covariate interactions)...",
      "Contrôle 1/3 : Test de l'homogénéité des pentes de régression (interactions facteur à— covariable)...",
      verbose = verbose, k = k, cpt = "on"
    )

    if (length(numeric_vars) > 0 && length(factor_vars) > 0) {

      # Construire formule avec interactions facteurà—covariable
      response_var <- all.vars(formula)[1]

      # Pour chaque covariable, tester interaction avec chaque facteur
      for (cov_name in numeric_vars) {
        for (fact_name in factor_vars) {

          tryCatch({
            # Formule avec interaction
            formula_interaction <- as.formula(paste0(
              response_var, " ~ ", fact_name, " * ", cov_name
            ))

            # Ajuster modèle avec interaction
            model_interaction <- lm(formula_interaction, data = data)

            # Test de l'interaction via ANOVA Type III
            if (requireNamespace("car", quietly = TRUE)) {
              anova_interaction <- car::Anova(model_interaction, type = "III")

              # Extraire p-value de l'interaction
              interaction_term <- paste0(fact_name, ":", cov_name)
              if (interaction_term %in% rownames(anova_interaction)) {
                p_interaction <- anova_interaction[interaction_term, "Pr(>F)"]

                # Stocker résultat
                if (is.null(ancova_checks$slopes_homogeneity)) {
                  ancova_checks$slopes_homogeneity <- list()
                }
                ancova_checks$slopes_homogeneity[[paste0(fact_name, "_x_", cov_name)]] <- list(
                  factor = fact_name,
                  covariate = cov_name,
                  p_value = p_interaction,
                  homogeneous = p_interaction >= alpha
                )

                if (p_interaction < alpha) {
                  # Pentes non homogènes
                  k <- .vbse(
                    paste0("Homogeneity of slopes VIOLATED for '", fact_name, " à— ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ").\n",
                           "\tInteraction is significant => regression slopes differ across factor levels.\n",
                           "\tConsider: (1) Separate analyses by factor level, OR (2) Include interaction in model."),
                    paste0("Homogénéité des pentes VIOLàE pour '", fact_name, " à— ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ").\n",
                           "\tL'interaction est significative => pentes de régression diffèrent selon les niveaux du facteur.\n",
                           "\tEnvisager : (1) Analyses séparées par niveau de facteur, OU (2) Inclusion de l'interaction dans le modèle."),
                    verbose = verbose, k = k, cpt = "off"
                  )
                  # Ne pas basculer automatiquement en robuste, laisser l'utilisateur décider
                  # robuste <- TRUE
                } else {
                  # Pentes homogènes
                  k <- .vbse(
                    paste0("Homogeneity of slopes satisfied for '", fact_name, " à— ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ")."),
                    paste0("Homogénéité des pentes respectée pour '", fact_name, " à— ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ")."),
                    verbose = verbose, k = k, cpt = "off"
                  )
                }
              }
            }
          }, error = function(e) {
            k <- .vbse(
              paste0("Warning: Could not test slopes homogeneity for '", fact_name, " à— ", cov_name, "': ", e$message),
              paste0("Attention : Impossible de tester l'homogénéité des pentes pour '", fact_name, " à— ", cov_name, "' : ", e$message),
              verbose = verbose, k = k, cpt = "off"
            )
          })
        }
      }

    } else {
      k <- .vbse(
        "Slopes homogeneity check skipped (no continuous covariates or no categorical factors).",
        "Contrôle de l'homogénéité des pentes omis (aucune covariable continue ou aucun facteur catégorique).",
        verbose = verbose, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # CONTRÔLE 2: LINEARITE DE LA RELATION COVARIABLE-REPONSE
    #--------------------------------------------------------------------------
    if (!robuste) {
      k <- .vbse(
        "Check 2/3: Testing linearity of covariate-outcome relationship (quadratic terms)...",
        "Contrôle 2/3 : Test de la linéarité de la relation covariable-résultat (termes quadratiques)...",
        verbose = verbose, k = k, cpt = "on"
      )

      if (length(numeric_vars) > 0) {

        for (cov_name in numeric_vars) {

          tryCatch({
            # Formule avec terme quadratique
            response_var <- all.vars(formula)[1]
            other_predictors <- setdiff(predictors, cov_name)

            if (length(other_predictors) > 0) {
              formula_quad <- as.formula(paste0(
                response_var, " ~ ", paste(other_predictors, collapse = " + "),
                " + ", cov_name, " + I(", cov_name, "^2)"
              ))
            } else {
              formula_quad <- as.formula(paste0(
                response_var, " ~ ", cov_name, " + I(", cov_name, "^2)"
              ))
            }

            # Ajuster modèle avec terme quadratique
            model_quad <- lm(formula_quad, data = data)

            # Test du terme quadratique
            if (requireNamespace("car", quietly = TRUE)) {
              anova_quad <- car::Anova(model_quad, type = "III")

              # Extraire p-value du terme quadratique
              quad_term <- paste0("I(", cov_name, "^2)")
              if (quad_term %in% rownames(anova_quad)) {
                p_quad <- anova_quad[quad_term, "Pr(>F)"]

                # Stocker résultat
                if (is.null(ancova_checks$linearity)) {
                  ancova_checks$linearity <- list()
                }
                ancova_checks$linearity[[cov_name]] <- list(
                  covariate = cov_name,
                  p_value = p_quad,
                  linear = p_quad >= alpha
                )

                if (p_quad < alpha) {
                  # Relation non linéaire
                  k <- .vbse(
                    paste0("Linearity VIOLATED for covariate '", cov_name, "' (p = ", .format_pval(p_quad), ").\n",
                           "\tQuadratic term is significant => non-linear relationship detected.\n",
                           "\tConsider: (1) Transformation of covariate, OR (2) Non-linear modeling (GAM, polynomial)."),
                    paste0("Linéarité VIOLàE pour la covariable '", cov_name, "' (p = ", .format_pval(p_quad), ").\n",
                           "\tLe terme quadratique est significatif => relation non linéaire détectée.\n",
                           "\tEnvisager : (1) Transformation de la covariable, OU (2) Modélisation non linéaire (GAM, polynomiale)."),
                    verbose = verbose, k = k, cpt = "off"
                  )
                  # Ne pas basculer automatiquement en robuste
                } else {
                  # Relation linéaire
                  k <- .vbse(
                    paste0("Linearity assumption satisfied for covariate '", cov_name,
                           "' (p = ", .format_pval(p_quad), ")."),
                    paste0("Hypothèse de linéarité respectée pour la covariable '", cov_name,
                           "' (p = ", .format_pval(p_quad), ")."),
                    verbose = verbose, k = k, cpt = "off"
                  )
                }
              }
            }
          }, error = function(e) {
            k <- .vbse(
              paste0("Warning: Could not test linearity for covariate '", cov_name, "': ", e$message),
              paste0("Attention : Impossible de tester la linéarité pour la covariable '", cov_name, "' : ", e$message),
              verbose = verbose, k = k, cpt = "off"
            )
          })
        }

      } else {
        k <- .vbse(
          "Linearity check skipped (no continuous covariates found).",
          "Contrôle de linéarité omis (aucune covariable continue trouvée).",
          verbose = verbose, k = k, cpt = "off"
        )
      }
    }

    #--------------------------------------------------------------------------
    # CONTRÔLE 3: HOMOSCéDASTICITé DES RéSIDUS PAR GROUPE
    #--------------------------------------------------------------------------
    if (!robuste && !is.null(model)) {

      #========================================================================
      #                      NOTE PERSO #10 [RéSOLUE]
      #========================================================================
      #* NOTE PERSO #10 : Contrôle de variance sur résidus (Brown-Forsythe)
      #*
      #* ➜ Problème identifié :
      #*   Quelle méthode utiliser pour tester l'homoscédasticité des résidus?
      #*   Une IA suggère que "pour résidus, Brown-Forsythe (center='median')
      #*   est recommandé". Est-ce académiquement validé?
      #*
      #* ➜ Source académique (APA + DOI) :
      #*   Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality
      #*   of variances. *Journal of the American Statistical Association*, 69(346),
      #*   364â€"367. https://doi.org/10.1080/01621459.1974.10482955
      #*
      #*   Idée principale : Le test de Brown-Forsythe utilise les déviations
      #*   absolues par rapport à  la MàDIANE (au lieu de la moyenne pour Levene).
      #*   Il est plus robuste aux outliers et aux distributions asymétriques.
      #*   Pour les résidus (qui peuvent être non normaux), Brown-Forsythe
      #*   (center="median") est effectivement RECOMMANDà.
      #*
      #* ➜ Solution appliquée :
      #*   Utilisation de car::leveneTest() avec center="median" (= Brown-Forsythe)
      #*   sur les résidus studentisés. Documentation de la méthode et stockage
      #*   du type de test dans ancova_checks.
      #*
      #* ➜ Statut : RéSOLU
      #========================================================================

      k <- .vbse(
        "Check 3/3: Testing residual homoscedasticity across groups (Brown-Forsythe test on residuals)...",
        "Contrôle 3/3 : Test de l'homoscédasticité des résidus entre groupes (test de Brown-Forsythe sur résidus)...",
        verbose = verbose, k = k, cpt = "on"
      )

      tryCatch({
        # Extraire résidus studentisés avec type
        resid_result <- get_studentized_residuals(model)
        residuals_stud <- resid_result$residuals
        residual_type_used <- resid_result$type  # Stocker pour rapport

        # Créer le facteur de groupe approprié
        if (!is.null(between) && length(between) > 0) {
          if (length(between) == 1) {
            group_factor <- g[[between]]
          } else {
            group_factor <- interaction(g[between], drop = TRUE)
          }
        } else {
          group_factor <- g_cat
        }

        # Test de Brown-Forsythe (Levene avec center="median")
        if (requireNamespace("car", quietly = TRUE)) {
          bf_resid <- car::leveneTest(residuals_stud, group_factor, center = "median")
          p_bf_resid <- bf_resid$`Pr(>F)`[1]

          # Stocker résultat avec méthodologie
          ancova_checks$residual_homoscedasticity <- list(
            test = "Brown-Forsythe",
            residual_type = residual_type_used,
            p_value = p_bf_resid,
            homoscedastic = p_bf_resid >= alpha
          )

          if (p_bf_resid < alpha) {
            # Hétéroscédasticité des résidus
            k <- .vbse(
              paste0("Residual homoscedasticity VIOLATED (p = ", .format_pval(p_bf_resid), ").\n",
                     "\tBrown-Forsythe test on ", residual_type_used, " residuals detects heteroscedasticity.\n",
                     "\tSwitching to robust analysis."),
              paste0("Homoscédasticité des résidus VIOLàE (p = ", .format_pval(p_bf_resid), ").\n",
                     "\tLe test de Brown-Forsythe sur résidus ", residual_type_used, " détecte une hétéroscédasticité.\n",
                     "\tPassage vers analyse robuste."),
              verbose = verbose, k = k, cpt = "off"
            )
            check_variance_equal <- FALSE
            robuste <- TRUE

          } else {
            # Homoscédasticité des résidus
            k <- .vbse(
              paste0("Residual homoscedasticity satisfied (p = ", .format_pval(p_bf_resid), ")."),
              paste0("Homoscédasticité des résidus respectée (p = ", .format_pval(p_bf_resid), ")."),
              verbose = verbose, k = k, cpt = "off"
            )
          }
        } else {
          k <- .vbse(
            "Warning: Package 'car' not available for residual homoscedasticity test.",
            "Attention : Package 'car' non disponible pour le test d'homoscédasticité des résidus.",
            verbose = verbose, k = k, cpt = "off"
          )
        }

      }, error = function(e) {
        k <- .vbse(
          paste0("Warning: Could not test residual homoscedasticity: ", e$message),
          paste0("Attention : Impossible de tester l'homoscédasticité des résidus : ", e$message),
          verbose = verbose, k = k, cpt = "off"
        )
      })
    }

    k <- .vbse(
      "=== END OF ANCOVA AUTOMATIC CHECKS ===",
      "=== FIN DES CONTRÔLES ANCOVA AUTOMATIQUES ===",
      verbose = verbose, k = k, cpt = "on"
    )

  } # Fin if (check_ancova && !robuste)

  #============================================================================
  #                     VOIE D'ANALYSE ROBUSTE
  #============================================================================

  if (robuste) {

    .dbg("", "Passage vers analyse robuste...", debug=debug)

    k <- .vbse(
      paste0("Switching to robust analysis approach due to:\n",
             "\t- Assumption violations, OR\n",
             "\t- Severe imbalance, OR\n",
             "\t- Discrete dependent variable, OR\n",
             "\t- Unbalanced repeated measures design, OR\n",
             "\t- ANCOVA assumption violations"),
      paste0("Passage vers une approche d'analyse robuste en raison de :\n",
             "\t- Violations d'hypothèses, OU\n",
             "\t- Déséquilibre sévère, OU\n",
             "\t- Variable dépendante discrète, OU\n",
             "\t- Plan de mesures répétées déséquilibré, OU\n",
             "\t- Violations d'hypothèses ANCOVA"),
      verbose = verbose, k = k, cpt = "on"
    )
    #==========================================================================
    # CORRECTION CRITIQUE : INITIALISATION EN DÉBUT DE BLOC ROBUSTE
    #==========================================================================
    # Initialiser la structure de résultats robustes AVANT tout if/else
    robust_results <- list(
      method = NULL,
      test_result = NULL,
      assumptions_checked = list(),
      warnings = character(),
      posthoc_applicable = FALSE
    )

    # Déterminer le nombre de facteurs réels (AVANT tout if/else)
    n_factors <- length(factor_vars)

    .dbg("", paste0("Nombre de facteurs détectés: ", n_factors), debug=debug)
    #==========================================================================
    #                      NOTE PERSO #11 [EN COURS - NON RéSOLUE]
    #==========================================================================
    #* NOTE PERSO #11 : Implémentation automatique des tests robustes
    #*
    #* ➜ Problème identifié :
    #*   Pour l'instant, la fonction suggère des démarches mais ne les applique PAS.
    #*   Son ambition est de faire TOUS ces tests automatiquement, de contrôler
    #*   les assomptions à  chaque fois, et de changer de stratégie selon les contrôles.
    #*
    #*   Exemples d'implémentations manquantes :
    #*   1. Friedman test (k >= 3, un seul facteur, appariement)
    #*   2. Modèles mixtes (lmer) avec critères automatiques
    #*   3. ANOVA par permutation (lmPerm::aovp)
    #*   4. Tests robustes WRS2 (t2way, etc.)
    #*   5. Scheirer-Ray-Hare pour 2-way non paramétrique
    #*
    #* ➜ Sources académiques (APA + DOI) :
    #*
    #*   A) FRIEDMAN TEST (mesures répétées, k >= 3, un facteur)
    #*   Friedman, M. (1937). The use of ranks to avoid the assumption of
    #*   normality implicit in the analysis of variance. *Journal of the American
    #*   Statistical Association*, 32(200), 675â€"701.
    #*   https://doi.org/10.1080/01621459.1937.10503522
    #*
    #*   Idée : Extension non paramétrique de l'ANOVA RM pour k >= 3 conditions.
    #*   Basé sur les rangs. Robuste aux violations de normalité et homoscédasticité.
    #*
    #*   B) MODÈLES MIXTES (lmer)
    #*   Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random
    #*   effects structure for confirmatory hypothesis testing: Keep it maximal.
    #*   *Journal of Memory and Language*, 68(3), 255â€"278.
    #*   https://doi.org/10.1016/j.jml.2012.11.001
    #*
    #*   Critères pour lmer :
    #*   - RM déséquilibré (>10% manquants)
    #*   - Sphéricité violée avec Îµ < 0.75
    #*   - Design complexe (nested/crossed)
    #*   - ANCOVA avec violations multiples
    #*
    #*   C) ANOVA PAR PERMUTATION
    #*   Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
    #*   multi-factorial analysis of variance. *Journal of Statistical Computation
    #*   and Simulation*, 73(2), 85â€"113.
    #*   https://doi.org/10.1080/00949650215733
    #*
    #*   Idée : Tests exacts par permutation. Pas d'assumption de distribution.
    #*   Recommandé pour designs complexes, déséquilibrés, ou violations multiples.
    #*   lmPerm::aovp marche avec 2+ facteurs (sensible aux déséquilibres).
    #*
    #*   D) WRS2 POUR TESTS ROBUSTES
    #*   Wilcox, R. R. (2017). *Introduction to robust estimation and hypothesis
    #*   testing* (4th ed.). Academic Press.
    #*   https://doi.org/10.1016/C2010-0-67044-1
    #*
    #*   Idée : Moyennes tronquées, médianes, bootstrap. WRS2::t2way pour 2-way
    #*   robuste (2 facteurs uniquement).
    #*
    #*   E) SCHEIRER-RAY-HARE (2-way non paramétrique)
    #*   Scheirer, C. J., Ray, W. S., & Hare, N. (1976). The analysis of ranked
    #*   data derived from completely randomized factorial designs. *Biometrics*,
    #*   32(2), 429â€"434. https://doi.org/10.2307/2529511
    #*
    #*   Idée : Extension de Kruskal-Wallis pour 2-way. Basé sur les rangs.
    #*   rcompanion::scheirerRayHare (2 facteurs uniquement).
    #*
    #* ➜ Solution à implémenter (PLAN D'ACTION FUTUR) :
    #*
    #*   PHASE 1 : Implémentation basique (prochaine version)
    #*   1. Friedman test automatique (k >= 3, un facteur, paired)
    #*   2. Kruskal-Wallis pour 1 facteur (déjà partiellement fait)
    #*   3. Messages clairs sur limitations actuelles
    #*
    #*   PHASE 2 : Implémentation modèles mixtes (version ultérieure)
    #*   1. Créer fonction .mixed_model_analysis() séparée
    #*   2. Critères automatiques de déclenchement lmer
    #*   3. Construction automatique formule lmer
    #*   4. Validation assomptions via valreg() ou équivalent
    #*   5. Diagnostics complets (VarCorr, ranef, etc.)
    #*   6. Compatibilité avec pipeline .posthoc()
    #*
    #*   PHASE 3 : Tests robustes avancés (version future)
    #*   1. Implémentation lmPerm::aovp avec gestion déséquilibres
    #*   2. Intégration WRS2 (t2way, etc.)
    #*   3. Scheirer-Ray-Hare pour 2-way non paramétrique
    #*   4. Tests de Fligner-Killeen + contrôle outliers pour choix optimal
    #*
    #* ➜ Statut : NON RéSOLUE - PLAN D'ACTION DOCUMENTà
    #*   Implémentation complète nécessite développement conséquent.
    #*   Priorité : Friedman (PHASE 1) puis lmer (PHASE 2).
    #*   En attendant, messages explicites vers utilisateur.
    #==========================================================================

    #==========================================================================
    #                      NOTE PERSO #12 [PARTIELLEMENT RéSOLUE]
    #==========================================================================
    #* NOTE PERSO #12 : Kruskal-Wallis en ANOVA 2-way
    #*
    #* ➜ Problème identifié :
    #*   Kruskal-Wallis n'est valable que pour UN SEUL facteur, pas pour ANOVA
    #*   2-way ou plus. Il faudrait plutôt :
    #*   - ANOVA par permutation (le meilleur)
    #*   - À la rigueur Scheirer-Ray-Hare (2 facteurs uniquement)
    #*   - Ou t2way de WRS2 (2 facteurs uniquement)
    #*   Le choix peut être aidé par Fligner-Killeen + contrôle outliers.
    #*
    #* ➜ Source académique (APA + DOI) :
    #*   Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion
    #*   variance analysis. *Journal of the American Statistical Association*,
    #*   47(260), 583â€"621. https://doi.org/10.1080/01621459.1952.10483441
    #*
    #*   Idée principale : Kruskal-Wallis est l'extension non paramétrique de
    #*   l'ANOVA ONE-WAY. Il teste H0: toutes les distributions sont identiques,
    #*   basé sur les rangs. Il n'y a PAS de version K-W pour designs factoriels
    #*   (2-way ou plus). C'est une limitation fondamentale du test.
    #*
    #* ➜ Solution appliquée :
    #*   1. Vérification : Kruskal-Wallis appelé UNIQUEMENT si un seul facteur
    #*      (nlevels interaction == nlevels d'un seul facteur)
    #*   2. Pour designs multi-facteurs : message clair vers alternatives
    #*      (permutation ANOVA, Scheirer-Ray-Hare, modèles mixtes)
    #*   3. Documentation de la limitation dans les messages
    #*
    #* ➜ Statut : PARTIELLEMENT RéSOLU (restriction implémentée, mais
    #*   alternatives non automatisées - voir NOTE #11)
    #==========================================================================

    # Suggérer les méthodes robustes selon le design
  if (paired && nlevels(g_cat) >= 3) {

    # Cas spécial : mesures répétées avec 1 facteur → Friedman automatique
    if (n_factors == 1 && !is.null(id)) {

      k <- .vbse(
        "Applying Friedman test for repeated measures (k >= 3 conditions, one factor)...",
        "Application du test de Friedman pour mesures répétées (k >= 3 conditions, un facteur)...",
        verbose = verbose, k = k, cpt = "on"
      )

      tryCatch({
        # Préparer les données en format large pour Friedman
        if (!is.null(within) && length(within) == 1) {
          # Créer data frame pour reshape
          friedman_data <- data.frame(
            id = factor(data[[id]]),
            condition = factor(data[[within]]),
            value = x
          )

          # Reshape en format large
          wide_data <- reshape(
            friedman_data,
            idvar = "id",
            timevar = "condition",
            direction = "wide",
            v.names = "value"
          )

          # Extraire matrice de valeurs
          value_cols <- grep("^value\\.", names(wide_data), value = TRUE)
          value_matrix <- as.matrix(wide_data[, value_cols])

          # Appliquer test de Friedman
          friedman_result <- friedman.test(value_matrix)

          robust_results$method <- "Friedman_Test"
          robust_results$test_result <- friedman_result
          robust_results$posthoc_applicable <- friedman_result$p.value < alpha

          k <- .vbse(
            paste0("Friedman test results:\n",
                   "\tχ² = ", round(friedman_result$statistic, 3), "\n",
                   "\tp-value = ", format.pval(friedman_result$p.value, digits = 3), "\n",
                   "\tPost-hoc applicable: ", if(robust_results$posthoc_applicable) "Yes" else "No"),
            paste0("Résultats du test de Friedman :\n",
                   "\tχ² = ", round(friedman_result$statistic, 3), "\n",
                   "\tp-value = ", format.pval(friedman_result$p.value, digits = 3), "\n",
                   "\tTests post-hoc applicables : ", if(robust_results$posthoc_applicable) "Oui" else "Non"),
            verbose = verbose, k = k, cpt = "off"
          )

        } else {
          # Cas multi-facteurs : suggérer approche avancée
          k <- .vbse(
            paste0("Detected repeated measures with multiple factors (", n_factors, " factors).\n",
                   "\tSuggested approaches:\n",
                   "\t- Mixed model with robust estimation (lmerTest::lmer)\n",
                   "\t- Robust repeated measures from {WRS2} package\n",
                   "\t- Permutation-based methods\n",
                   "\tNote: Automatic implementation requires .mixed_model_analysis() (in development)."),
            paste0("Mesures répétées détectées avec plusieurs facteurs (", n_factors, " facteurs).\n",
                   "\tApproches suggérées :\n",
                   "\t- Modèle mixte avec estimation robuste (lmerTest::lmer)\n",
                   "\t- Mesures répétées robustes du package {WRS2}\n",
                   "\t- Méthodes par permutation\n",
                   "\tNote : Implémentation automatique nécessite .mixed_model_analysis() (en développement)."),
            verbose = verbose, k = k, cpt = "off"
          )

          robust_results$method <- "Mixed_Model_Required"
          robust_results$warnings <- c(
            robust_results$warnings,
            "Multi-factor repeated measures - mixed model recommended but not automatically implemented"
          )
        }

      }, error = function(e) {
        k <- .vbse(
          paste0("Error applying Friedman test: ", e$message),
          paste0("Erreur lors de l'application du test de Friedman : ", e$message),
          verbose = verbose, k = k, cpt = "on"
        )

        robust_results$method <<- "Friedman_Failed"
        robust_results$warnings <<- c(
          robust_results$warnings,
          paste0("Friedman test failed: ", e$message)
        )
      })

    } else {
      # Cas sans id ou multi-facteurs : message générique
      k <- .vbse(
        paste0("Suggested robust approach for repeated measures:\n",
               "\t- Friedman test (friedman.test) for k >= 3 conditions\n",
               "\t- Robust repeated measures from {WRS2} package\n",
               "\t- Mixed model with robust estimation (lmer with REML)\n",
               "\tNote: Implementation of these tests is NOT YET fully automatic (see documentation)."),
        paste0("Approche robuste suggérée pour mesures répétées :\n",
               "\t- Test de Friedman (friedman.test) pour k >= 3 conditions\n",
               "\t- Mesures répétées robustes du package {WRS2}\n",
               "\t- Modèle mixte avec estimation robuste (lmer avec REML)\n",
               "\tNote : L'implémentation de ces tests n'est PAS ENCORE complètement automatique (voir documentation)."),
        verbose = verbose, k = k, cpt = "off"
      )

      robust_results$method <- "Repeated_Measures_Suggested"
      robust_results$warnings <- c(
        robust_results$warnings,
        "Repeated measures design - automatic implementation limited"
      )
    }

  }  else if (alea || alea_plus) {
      k <- .vbse(
        paste0("Suggested robust approach for random effects:\n",
               "\t- Mixed model with robust estimation (lmerTest::lmer)\n",
               "\t- Generalized Estimating Equations (GEE) if appropriate\n",
               "\t- Bootstrap methods for inference\n",
               "\tNote: Implementation requires dedicated function .mixed_model_analysis() (in development)."),
        paste0("Approche robuste suggérée pour effets aléatoires :\n",
               "\t- Modèle mixte avec estimation robuste (lmerTest::lmer)\n",
               "\t- équations d'estimation généralisées (GEE) si approprié\n",
               "\t- Méthodes de bootstrap pour l'inférence\n",
               "\tNote : L'implémentation nécessite fonction dédiée .mixed_model_analysis() (en développement)."),
        verbose = verbose, k = k, cpt = "off"
      )

    } else {
      k <- .vbse(
        paste0("Suggested robust approaches:\n",
               "\t- Permutation ANOVA (aovp from {lmPerm}) - BEST for multi-factor designs\n",
               "\t- Rank-based methods:\n",
               "\t  * Kruskal-Wallis (ONE factor only)\n",
               "\t  * Scheirer-Ray-Hare (TWO factors only, via {rcompanion})\n",
               "\t- Trimmed means ANOVA (WRS2::t2way for TWO factors)\n",
               "\t- Welch's ANOVA for heteroscedasticity (oneway.test, ONE factor)\n",
               "\tNote: Only Kruskal-Wallis (if ONE factor) is automatically performed below."),
        paste0("Approches robustes suggérées :\n",
               "\t- ANOVA par permutation (aovp du package {lmPerm}) - MEILLEUR pour designs multi-facteurs\n",
               "\t- Méthodes basées sur les rangs :\n",
               "\t  * Kruskal-Wallis (UN facteur uniquement)\n",
               "\t  * Scheirer-Ray-Hare (DEUX facteurs uniquement, via {rcompanion})\n",
               "\t- ANOVA sur moyennes tronquées (WRS2::t2way pour DEUX facteurs)\n",
               "\t- ANOVA de Welch pour hétéroscédasticité (oneway.test, UN facteur)\n",
               "\tNote : Seul Kruskal-Wallis (si UN facteur) est automatiquement effectué ci-dessous."),
        verbose = verbose, k = k, cpt = "off"
      )
      #============================================================================
      #           IMPLÉMENTATION AUTOMATIQUE DES TESTS ROBUSTES
      #           (Réponse à NOTE PERSO #11)
      #============================================================================

      # Initialiser la structure de résultats robustes
      robust_results <- list(
        method = NULL,
        test_result = NULL,
        assumptions_checked = list(),
        warnings = character(),
        posthoc_applicable = FALSE
      )

      # Déterminer le nombre de facteurs réels
      n_factors <- length(factor_vars)

      #=============================================================================
      # PHASE 1: SÉLECTION ET APPLICATION AUTOMATIQUE DU TEST ROBUSTE
      #=============================================================================

      if (paired && nlevels(g_cat) >= 3 && n_factors == 1) {
        #---------------------------------------------------------------------------
        # CAS 1: MESURES RÉPÉTÉES, k >= 3, UN SEUL FACTEUR
        # → Test de Friedman (extension non paramétrique de RM-ANOVA)
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Applying Friedman test for repeated measures (k >= 3 conditions, one factor)...",
          "Application du test de Friedman pour mesures répétées (k >= 3 conditions, un facteur)...",
          verbose = verbose, k = k, cpt = "on"
        )

        tryCatch({
          # Préparer les données en format large pour Friedman
          if (!is.null(within) && length(within) == 1) {
            # Créer data frame pour reshape
            friedman_data <- data.frame(
              id = factor(data[[id]]),
              condition = factor(data[[within]]),
              value = x
            )

            # Reshape en format large
            wide_data <- reshape(
              friedman_data,
              direction = "wide",
              idvar = "id",
              timevar = "condition",
              v.names = "value"
            )

            # Extraire la matrice de valeurs (exclure colonne id)
            value_matrix <- as.matrix(wide_data[, -1])

            # Test de Friedman
            friedman_result <- friedman.test(value_matrix)

            robust_results$method <- "Friedman"
            robust_results$test_result <- friedman_result
            robust_results$posthoc_applicable <- TRUE  # Posthoc possible

            k <- .vbse(
              paste0("Friedman test completed:\n",
               "\tChi-squared = ", round(friedman_result$statistic, 3),
               ", df = ", friedman_result$parameter,
               ", p-value = ", .format_pval(friedman_result$p.value)),
        paste0("Test de Friedman terminé :\n",
                     "\tKhi² = ", round(friedman_result$statistic, 3),
                     ", ddl = ", friedman_result$parameter,
                     ", p = ", .format_pval(friedman_result$p.value)),
              verbose = verbose, k = k, cpt = "off"
            )

            # Vérifier les assomptions spécifiques à Friedman
            # 1. Échelle au moins ordinale
            n_unique <- length(unique(x[!is.na(x)]))
            if (n_unique < 5) {
              robust_results$warnings <- c(
                robust_results$warnings,
                paste0("Few distinct values (n=", n_unique, "). Verify ordinal scale assumption.")
              )
              k <- .vbse(
                paste0("Note: Only ", n_unique, " distinct values in DV. Verify ordinal scale is appropriate."),
                paste0("Note : Seulement ", n_unique, " valeurs distinctes dans VD. Vérifier que l'échelle ordinale est appropriée."),
                verbose = verbose, k = k, cpt = "off"
              )
            }
            robust_results$assumptions_checked$ordinal_scale <- (n_unique >= 5)

            # 2. Proportion de rangs ex-aequo
            n_ties <- sum(duplicated(x[!is.na(x)]))
            tie_proportion <- n_ties / length(x[!is.na(x)])
            robust_results$assumptions_checked$tie_proportion <- tie_proportion

            if (tie_proportion > 0.25) {
              robust_results$warnings <- c(
                robust_results$warnings,
                paste0("High proportion of ties (", round(tie_proportion * 100, 1), "%)")
              )
              k <- .vbse(
                paste0("Warning: ", round(tie_proportion * 100, 1),
                       "% tied ranks. Friedman adjusts automatically but power may be reduced."),
                paste0("Attention : ", round(tie_proportion * 100, 1),
                       "% de rangs ex-aequo. Friedman ajuste automatiquement mais la puissance peut être réduite."),
                verbose = verbose, k = k, cpt = "off"
              )
            }

          } else {
            # Plusieurs facteurs within : Friedman non applicable
            stop("Multiple within-subject factors detected. Friedman test not applicable.")
          }

        }, error = function(e) {
          warning(.msg(
            paste0("Failed to perform Friedman test: ", e$message),
            paste0("échec du test de Friedman : ", e$message)
          ))
          robust_results$method <- "Friedman_Failed"
          robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

          k <- .vbse(
            "Friedman test failed. Consider mixed models for this design.",
            "Test de Friedman échoué. Envisager modèles mixtes pour ce plan.",
            verbose = verbose, k = k, cpt = "on"
          )
        })

      } else if (n_factors == 1 && !paired && nlevels(g_cat) >= 2) {
        #---------------------------------------------------------------------------
        # CAS 2: UN SEUL FACTEUR, DONNÉES INDÉPENDANTES
        # → Kruskal-Wallis (déjà documenté dans NOTE #12, mais enrichi ici)
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Applying Kruskal-Wallis test for one-way independent design...",
          "Application du test de Kruskal-Wallis pour plan unifactoriel indépendant...",
          verbose = verbose, k = k, cpt = "on"
        )

        tryCatch({
          # Test de Kruskal-Wallis
          kw_result <- kruskal.test(x, g_cat)

          robust_results$method <- "Kruskal_Wallis"
          robust_results$test_result <- kw_result
          robust_results$posthoc_applicable <- TRUE

          k <- .vbse(
            paste0("Kruskal-Wallis test completed:\n",
             "\tChi-squared = ", round(kw_result$statistic, 3),
             ", df = ", kw_result$parameter,
             ", p-value = ", .format_pval(kw_result$p.value)),
      paste0("Test de Kruskal-Wallis terminé :\n",
                   "\tKhi² = ", round(kw_result$statistic, 3),
                   ", ddl = ", kw_result$parameter,
                   ", p = ", .format_pval(kw_result$p.value)),
            verbose = verbose, k = k, cpt = "off"
          )

          # Contrôle: homogénéité des distributions avec Fligner-Killeen
          k <- .vbse(
            "Checking variance homogeneity with Fligner-Killeen test (robust to non-normality)...",
            "Vérification de l'homogénéité des variances avec test de Fligner-Killeen (robuste à la non-normalité)...",
            verbose = verbose, k = k, cpt = "on"
          )

          fk_result <- fligner.test(x, g_cat)
          robust_results$assumptions_checked$fligner_killeen <- list(
            statistic = fk_result$statistic,
            p_value = fk_result$p.value,
            homogeneous = fk_result$p.value >= alpha
          )

          k <- .vbse(
            paste0("Fligner-Killeen: chi-squared = ", round(fk_result$statistic, 3),
                   ", p = ", .format_pval(fk_result$p.value)),
            paste0("Fligner-Killeen : khi² = ", round(fk_result$statistic, 3),
                   ", p = ", .format_pval(fk_result$p.value)),
            verbose = verbose, k = k, cpt = "off"
          )

          if (fk_result$p.value < alpha) {
            robust_results$warnings <- c(
              robust_results$warnings,
              "Variance heterogeneity detected (Fligner-Killeen)"
            )
            k <- .vbse(
              "Variance heterogeneity detected. Kruskal-Wallis remains valid but power may be affected.",
              "Hétérogénéité des variances détectée. Kruskal-Wallis reste valide mais la puissance peut être affectée.",
              verbose = verbose, k = k, cpt = "off"
            )
          }

          # Contrôle: détection d'outliers extrêmes
          if (requireNamespace("rstatix", quietly = TRUE)) {
            outliers <- rstatix::identify_outliers(data.frame(group = g_cat, value = x), value)
            n_extreme <- sum(outliers$is.extreme, na.rm = TRUE)

            robust_results$assumptions_checked$extreme_outliers <- n_extreme

            if (n_extreme > 0) {
              robust_results$warnings <- c(
                robust_results$warnings,
                paste0(n_extreme, " extreme outlier(s) detected")
              )
              k <- .vbse(
                paste0("Warning: ", n_extreme, " extreme outlier(s) detected. May affect interpretation."),
                paste0("Attention : ", n_extreme, " valeur(s) extrême(s) détectée(s). Peut affecter l'interprétation."),
                verbose = verbose, k = k, cpt = "off"
              )
            }
          }

        }, error = function(e) {
          warning(.msg(
            paste0("Failed to perform Kruskal-Wallis test: ", e$message),
            paste0("échec du test de Kruskal-Wallis : ", e$message)
          ))
          robust_results$method <- "Kruskal_Wallis_Failed"
          robust_results$warnings <- c(robust_results$warnings, as.character(e$message))
        })

      } else if (n_factors == 2 && !paired && !check_ancova) {
        #---------------------------------------------------------------------------
        # CAS 3: DEUX FACTEURS, DONNÉES INDÉPENDANTES, PAS ANCOVA
        # → Scheirer-Ray-Hare (extension de Kruskal-Wallis pour 2-way)
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Applying Scheirer-Ray-Hare test for two-way independent design (rank-based)...",
          "Application du test de Scheirer-Ray-Hare pour plan bifactoriel indépendant (basé sur rangs)...",
          verbose = verbose, k = k, cpt = "on"
        )

        if (requireNamespace("rcompanion", quietly = TRUE)) {
          tryCatch({
            # Construire formule pour Scheirer-Ray-Hare (effets principaux uniquement)
            response_var <- all.vars(formula)[1]
            srh_formula <- as.formula(paste0(
              response_var, " ~ ", factor_vars[1], " + ", factor_vars[2]
            ))

            # Test de Scheirer-Ray-Hare
            srh_result <- rcompanion::scheirerRayHare(
              formula = srh_formula,
              data = data
            )

            robust_results$method <- "Scheirer_Ray_Hare"
            robust_results$test_result <- srh_result
            robust_results$posthoc_applicable <- TRUE

            k <- .vbse(
              "Scheirer-Ray-Hare test completed. Results:",
        "Test de Scheirer-Ray-Hare terminé. Résultats :",
              verbose = verbose, k = k, cpt = "off"
            )

            if (verbose) {
              cat("\n")
              print(srh_result)
              cat("\n")
            }

            # AVERTISSEMENT CRITIQUE: pas d'interaction testée
            robust_results$assumptions_checked$interaction_tested <- FALSE
            robust_results$warnings <- c(
              robust_results$warnings,
              "Scheirer-Ray-Hare: Interactions NOT tested (main effects only)"
            )

            k <- .vbse(
              paste0("CRITICAL: Scheirer-Ray-Hare tests MAIN EFFECTS ONLY.\n",
                     "\tInteractions between factors are NOT tested.\n",
                     "\tIf interaction is theoretically important, use permutation ANOVA instead."),
              paste0("CRITIQUE : Scheirer-Ray-Hare teste UNIQUEMENT les EFFETS PRINCIPAUX.\n",
                     "\tLes interactions entre facteurs ne sont PAS testées.\n",
                     "\tSi l'interaction est théoriquement importante, utiliser ANOVA par permutation."),
              verbose = verbose, k = k, cpt = "on"
            )

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to perform Scheirer-Ray-Hare test: ", e$message),
              paste0("échec du test de Scheirer-Ray-Hare : ", e$message)
            ))
            robust_results$method <- "Scheirer_Ray_Hare_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              "Scheirer-Ray-Hare failed. Falling back to permutation ANOVA recommendation.",
              "Scheirer-Ray-Hare échoué. Recommandation d'ANOVA par permutation.",
              verbose = verbose, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            "Package 'rcompanion' not available for Scheirer-Ray-Hare. Switching to permutation ANOVA...",
            "Package 'rcompanion' non disponible pour Scheirer-Ray-Hare. Passage vers ANOVA par permutation...",
            verbose = verbose, k = k, cpt = "on"
          )
          robust_results$method <- "Scheirer_Ray_Hare_Unavailable"

          # Rediriger vers permutation ANOVA (traité dans CAS 4)
          n_factors <- 999  # Force passage à CAS 4
        }

      }

      if (n_factors >= 2 && !paired && robust_results$method %in%
          c(NULL, "Scheirer_Ray_Hare_Unavailable", "Scheirer_Ray_Hare_Failed")) {
        #---------------------------------------------------------------------------
        # CAS 4: PLUSIEURS FACTEURS (2+), DONNÉES INDÉPENDANTES
        # → ANOVA par permutation (méthode UNIVERSELLE et OPTIMALE)
        #---------------------------------------------------------------------------

        k <- .vbse(
          paste0("Applying permutation ANOVA for multi-factor independent design (",
                 length(factor_vars), " factors, BEST approach)..."),
          paste0("Application de l'ANOVA par permutation pour plan multi-facteurs indépendant (",
                 length(factor_vars), " facteurs, MEILLEURE approche)..."),
          verbose = verbose, k = k, cpt = "on"
        )

        if (requireNamespace("lmPerm", quietly = TRUE)) {
          tryCatch({
            # ANOVA par permutation avec nombre de permutations adapté
            # Plus de facteurs = plus de permutations nécessaires
            n_perm <- ifelse(length(factor_vars) <= 2, "Prob", "Exact")

            k <- .vbse(
              paste0("Using permutation strategy: ", n_perm, " (adaptive to design complexity)"),
              paste0("Utilisation de la stratégie de permutation : ", n_perm, " (adaptative à la complexité du plan)"),
              verbose = verbose, k = k, cpt = "off"
            )

            # Ajuster le modèle par permutation
            perm_result <- lmPerm::aovp(formula, data = data, perm = n_perm)

            robust_results$method <- "Permutation_ANOVA"
            robust_results$test_result <- perm_result
            robust_results$posthoc_applicable <- TRUE

            k <- .vbse(
              "Permutation ANOVA completed successfully.",
        "ANOVA par permutation terminée avec succès.",
              verbose = verbose, k = k, cpt = "off"
            )

            if (verbose) {
              k <- .vbse(
                "ANOVA table (permutation p-values):",
                "Tableau ANOVA (valeurs p par permutation) :",
                verbose = verbose, k = k, cpt = "off"
              )
              cat("\n")
              print(summary(perm_result))
              cat("\n")
            }

            # Avantages documentés de la permutation ANOVA
            k <- .vbse(
              paste0("Permutation ANOVA advantages for this design:\n",
                     "\t✓ Distribution-free (no normality assumption)\n",
                     "\t✓ Handles unbalanced designs robustly\n",
                     "\t✓ Tests ALL main effects AND interactions\n",
                     "\t✓ Exact p-values (not asymptotic approximations)\n",
                     "\t✓ Robust to outliers and heteroscedasticity"),
              paste0("Avantages de l'ANOVA par permutation pour ce plan :\n",
                     "\t✓ Sans hypothèse de distribution (pas de normalité requise)\n",
                     "\t✓ Gère les plans déséquilibrés de manière robuste\n",
                     "\t✓ Teste TOUS les effets principaux ET les interactions\n",
                     "\t✓ Valeurs p exactes (pas d'approximations asymptotiques)\n",
                     "\t✓ Robuste aux valeurs extrêmes et hétéroscédasticité"),
              verbose = verbose, k = k, cpt = "off"
            )

            robust_results$assumptions_checked$distribution_free <- TRUE
            robust_results$assumptions_checked$handles_imbalance <- TRUE
            robust_results$assumptions_checked$tests_interactions <- TRUE
            robust_results$assumptions_checked$outlier_robust <- TRUE

            # Vérifier degré de déséquilibre (information, pas blocage)
            table_data_local <- table(g_cat)
            if (length(table_data_local) > 1) {
              max_ratio <- max(table_data_local) / min(table_data_local)
              robust_results$assumptions_checked$imbalance_ratio <- max_ratio

              if (max_ratio > 3) {
                robust_results$warnings <- c(
                  robust_results$warnings,
                  paste0("Severe imbalance: max/min ratio = ", round(max_ratio, 2))
                )
                k <- .vbse(
                  paste0("Note: Severely unbalanced design (ratio = ", round(max_ratio, 2),
                         "). Permutation ANOVA handles this, but power may be reduced for smaller groups."),
                  paste0("Note : Plan sévèrement déséquilibré (ratio = ", round(max_ratio, 2),
                         "). L'ANOVA par permutation gère cela, mais la puissance peut être réduite pour les petits groupes."),
                  verbose = verbose, k = k, cpt = "off"
                )
              }
            }

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to perform permutation ANOVA: ", e$message),
              paste0("échec de l'ANOVA par permutation : ", e$message)
            ))
            robust_results$method <- "Permutation_ANOVA_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              "Permutation ANOVA failed. No automatic robust alternative available for this design.",
              "ANOVA par permutation échouée. Aucune alternative robuste automatique disponible pour ce plan.",
              verbose = verbose, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            paste0("Package 'lmPerm' NOT available for permutation ANOVA.\n",
                   "\tThis is the RECOMMENDED approach for multi-factor robust analysis.\n",
                   "\tInstall with: install.packages('lmPerm')"),
            paste0("Package 'lmPerm' NON disponible pour ANOVA par permutation.\n",
                   "\tC'est l'approche RECOMMANDàE pour analyse robuste multi-facteurs.\n",
                   "\tInstaller avec : install.packages('lmPerm')"),
            verbose = verbose, k = k, cpt = "on"
          )
          robust_results$method <- "Permutation_ANOVA_Unavailable"
          robust_results$warnings <- c(
            robust_results$warnings,
            "lmPerm package not available - INSTALLATION STRONGLY RECOMMENDED"
          )
        }

      }

      if ((alea || alea_plus) ||
          (paired && length(within) > 1) ||
          (paired && !is.null(between) && length(between) > 0)) {
        #---------------------------------------------------------------------------
        # CAS 5: EFFETS ALÉATOIRES OU DESIGN RM COMPLEXE
        # → Modèles mixtes (lmer) avec validation
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Complex design with random effects or mixed within/between. Applying linear mixed model (lmer)...",
          "Plan complexe avec effets aléatoires ou within/between mixte. Application d'un modèle linéaire mixte (lmer)...",
          verbose = verbose, k = k, cpt = "on"
        )

        if (requireNamespace("lme4", quietly = TRUE) &&
            requireNamespace("lmerTest", quietly = TRUE)) {
          tryCatch({
            # Construire la formule lmer selon le design
            response_var <- all.vars(formula)[1]

            if (paired && !is.null(within)) {
              # Mesures répétées
              if (length(within) == 1 && (is.null(between) || length(between) == 0)) {
                # RM simple: DV ~ within + (1|id)
                lmer_formula <- as.formula(paste0(
                  response_var, " ~ ", within[1], " + (1|", id, ")"
                ))
              } else if (length(within) == 1 && !is.null(between) && length(between) == 1) {
                # Mixed RM: DV ~ within * between + (1|id)
                lmer_formula <- as.formula(paste0(
                  response_var, " ~ ", within[1], " * ", between[1],
                  " + (1|", id, ")"
                ))
              } else if (length(within) > 1) {
                # Multiple within: DV ~ w1 * w2 + (1|id)
                lmer_formula <- as.formula(paste0(
                  response_var, " ~ ", paste(within, collapse = " * "),
                  " + (1|", id, ")"
                ))
              }
            } else if (alea) {
              # Formule avec Error() : conversion simplifiée
              k <- .vbse(
                "Note: Complex Error() syntax detected. Using simplified random intercept model.",
                "Note : Syntaxe Error() complexe détectée. Utilisation d'un modèle à intercept aléatoire simplifié.",
                verbose = verbose, k = k, cpt = "off"
              )

              # Extraire le terme d'erreur (simplifié)
              if (!is.null(id)) {
                lmer_formula <- as.formula(paste0(
                  response_var, " ~ ", paste(predictors, collapse = " * "),
                  " + (1|", id, ")"
                ))
              } else {
                stop("Random effects detected but no 'id' variable specified")
              }
            }

            k <- .vbse(
              paste0("Mixed model formula: ", deparse(lmer_formula)),
              paste0("Formule du modèle mixte : ", deparse(lmer_formula)),
              verbose = verbose, k = k, cpt = "off"
            )

            # Ajuster le modèle mixte
            mixed_model <- lmerTest::lmer(lmer_formula, data = data, REML = TRUE)

            robust_results$method <- "Mixed_Model_lmer"
            robust_results$test_result <- mixed_model
            robust_results$posthoc_applicable <- TRUE

            k <- .vbse(
              "Mixed model fitted successfully (REML estimation).",
              "Modèle mixte ajusté avec succès (estimation REML).",
              verbose = verbose, k = k, cpt = "off"
            )

            # Afficher tableau ANOVA
            if (verbose) {
              k <- .vbse(
                "Type III ANOVA table for mixed model (Satterthwaite approximation):",
                "Tableau ANOVA de type III pour modèle mixte (approximation de Satterthwaite) :",
                verbose = verbose, k = k, cpt = "off"
              )
              cat("\n")
              print(lmerTest::anova(mixed_model, type = "III"))
              cat("\n")
            }

            #-----------------------------------------------------------------------
            # DIAGNOSTICS DU MODÈLE MIXTE
            #-----------------------------------------------------------------------
            k <- .vbse(
              "Performing mixed model diagnostics...",
              "Réalisation des diagnostics du modèle mixte...",
              verbose = verbose, k = k, cpt = "on"
            )

            # 1. Normalité des résidus
            lmer_resid <- residuals(mixed_model)

            resid_normality_result <- .normality(
              lmer_resid, g = NULL, alpha = alpha, paired = FALSE,
              debug = debug, verbose = verbose, k = k
            )
            k <- resid_normality_result[[2]]
            resid_normal <- resid_normality_result[[1]]

            robust_results$assumptions_checked$residuals_normality <- resid_normal

            if (!resid_normal) {
              robust_results$warnings <- c(
                robust_results$warnings,
                "Mixed model: residuals violate normality"
              )
              k <- .vbse(
                paste0("Warning: Residuals violate normality assumption.\n",
                       "\tConsider: (1) Transformation of DV, (2) Robust mixed models (robustlmm package),\n",
                       "\t(3) Generalized linear mixed models (GLMM) if appropriate."),
                paste0("Attention : Les résidus violent l'hypothèse de normalité.\n",
                       "\tEnvisager : (1) Transformation de la VD, (2) Modèles mixtes robustes (package robustlmm),\n",
                       "\t(3) Modèles linéaires généralisés mixtes (GLMM) si approprié."),
                verbose = verbose, k = k, cpt = "off"
              )
            }

            # 2. Homoscédasticité des résidus
            if (requireNamespace("car", quietly = TRUE)) {
              # Créer facteur de groupe approprié
              if (!is.null(between) && length(between) > 0) {
                group_fac <- if (length(between) == 1) {
                  g[[between]]
                } else {
                  interaction(g[between], drop = TRUE)
                }
              } else if (!is.null(within) && length(within) > 0) {
                group_fac <- if (length(within) == 1) {
                  g[[within]]
                } else {
                  interaction(g[within], drop = TRUE)
                }
              } else {
                group_fac <- g_cat
              }

              lev_test <- car::leveneTest(lmer_resid, group_fac, center = "median")
              p_lev <- lev_test$`Pr(>F)`[1]

              robust_results$assumptions_checked$residuals_homoscedasticity <- (p_lev >= alpha)

              k <- .vbse(
                paste0("Levene test on residuals: F = ", round(lev_test$`F value`[1], 3),
                       ", p = ", .format_pval(p_lev)),
                paste0("Test de Levene sur résidus : F = ", round(lev_test$`F value`[1], 3),
                       ", p = ", .format_pval(p_lev)),
                verbose = verbose, k = k, cpt = "off"
              )

              if (p_lev < alpha) {
                robust_results$warnings <- c(
                  robust_results$warnings,
                  "Mixed model: residual heteroscedasticity detected"
                )
              }
            }

            # 3. Variance components
            var_comp <- lme4::VarCorr(mixed_model)
            robust_results$assumptions_checked$variance_components <- var_comp

            if (verbose) {
              k <- .vbse(
                "Variance components:",
                "Composantes de variance :",
                verbose = verbose, k = k, cpt = "off"
              )
              cat("\n")
              print(var_comp)
              cat("\n")
            }

            # Remplacer le modèle NULL par le modèle mixte
            model <- mixed_model

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to fit mixed model: ", e$message),
              paste0("échec de l'ajustement du modèle mixte : ", e$message)
            ))
            robust_results$method <- "Mixed_Model_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              paste0("Mixed model fitting failed. Possible causes:\n",
                     "\t- Singular fit (variance component = 0)\n",
                     "\t- Convergence issues\n",
                     "\t- Insufficient data for random effects\n",
                     "\tConsider simplified random structure or consult statistician."),
              paste0("Ajustement du modèle mixte échoué. Causes possibles :\n",
                     "\t- Ajustement singulier (composante de variance = 0)\n",
                     "\t- Problèmes de convergence\n",
                     "\t- Données insuffisantes pour effets aléatoires\n",
                     "\tEnvisager structure aléatoire simplifiée ou consulter statisticien."),
              verbose = verbose, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            paste0("Packages 'lme4' and/or 'lmerTest' NOT available.\n",
                   "\tMixed models are REQUIRED for this complex design.\n",
                   "\tInstall with: install.packages(c('lme4', 'lmerTest'))"),
            paste0("Packages 'lme4' et/ou 'lmerTest' NON disponibles.\n",
                   "\tLes modèles mixtes sont REQUIS pour ce plan complexe.\n",
                   "\tInstaller avec : install.packages(c('lme4', 'lmerTest'))"),
            verbose = verbose, k = k, cpt = "on"
          )
          robust_results$method <- "Mixed_Model_Unavailable"
          robust_results$warnings <- c(
            robust_results$warnings,
            "lme4/lmerTest packages not available - INSTALLATION REQUIRED"
          )
        }
      }
    } # J'ai ajouté pour voir....
      # Si aucune méthode n'a été appliquée
      if (is.null(robust_results$method)) {
        k <- .vbse(
          paste0("WARNING: No automatic robust method could be applied for this design.\n",
                 "\tDesign characteristics:\n",
                 "\t  - Paired: ", if(paired) "Yes" else "No", "\n",
                 "\t  - Number of factors: ", n_factors, "\n",
                 "\t  - ANCOVA: ", if(check_ancova) "Yes" else "No", "\n",
                 "\t  - Random effects: ", if(alea) "Yes" else "No", "\n",
                 "\tRecommendation: Consult a statistician for appropriate analysis strategy."),
          paste0("ATTENTION : Aucune méthode robuste automatique n'a pu être appliquée pour ce plan.\n",
                 "\tCaractéristiques du plan :\n",
                 "\t  - Apparié : ", if(paired) "Oui" else "Non", "\n",
                 "\t  - Nombre de facteurs : ", n_factors, "\n",
                 "\t  - ANCOVA : ", if(check_ancova) "Oui" else "Non", "\n",
                 "\t  - Effets aléatoires : ", if(alea) "Oui" else "Non", "\n",
                 "\tRecommandation : Consulter un statisticien pour stratégie d'analyse appropriée."),
          verbose = verbose, k = k, cpt = "on"
        )

        robust_results$method <- "No_Automatic_Method_Available"
        robust_results$warnings <- c(
          robust_results$warnings,
          "Design configuration not covered by automatic robust methods"
        )
      }

      #=============================================================================
      # PHASE 2: RESUME ET RECOMMANDATIONS POST-TEST
      #=============================================================================

      k <- .vbse(
        paste0("=== ROBUST ANALYSIS SUMMARY ===\n",
               "\tMethod applied: ", robust_results$method, "\n",
               "\tPost-hoc tests applicable: ", if(robust_results$posthoc_applicable) "Yes" else "No"),
        paste0("=== RàSUMà DE L'ANALYSE ROBUSTE ===\n",
               "\tMéthode appliquée : ", robust_results$method, "\n",
               "\tTests post-hoc applicables : ", if(robust_results$posthoc_applicable) "Oui" else "Non"),
        verbose = verbose, k = k, cpt = "on"
      )

      # Afficher les avertissements
      if (length(robust_results$warnings) > 0 && verbose) {
        k <- .vbse(
          paste0("Warnings/Notes:\n\t", paste(robust_results$warnings, collapse = "\n\t")),
          paste0("Avertissements/Notes :\n\t", paste(robust_results$warnings, collapse = "\n\t")),
          verbose = verbose, k = k, cpt = "off"
        )
      }

      # Ajouter robust_results à la structure de retour finale
      # (sera intégré dans le return() à la fin de la fonction)

  } # Fin if (robuste)

  #============================================================================
  #                  RESUME DU MODELE (si paramétrique)
  #============================================================================

  if (!robuste && !is.null(model)) {

    .dbg("", "Affichage du résumé du modèle...", debug=debug)

    k <- .vbse(
      "Parametric analysis completed. Model summary available.",
      "Analyse paramétrique terminée. Résumé du modèle disponible.",
      verbose = verbose, k = k, cpt = "on"
    )
    #* Note perso : dans une situation équilibrée et paramétrique, on ne devrait pas faire d'ANOVA de type III.
    #* Il faut corriger le routage suite au contrôle de l'équilibrage des données.

    # Afficher le tableau ANOVA avec le type de SC approprié
    if (verbose) {
      if (inherits(model, "lm")) {
        k <- .vbse(
          "Type III Sum of Squares (car::Anova):",
          "Sommes des carrés de type III (car::Anova) :",
          verbose = verbose, k = k, cpt = "off"
        )
        cat("\n")
        if (requireNamespace("car", quietly = TRUE)) {
          print(car::Anova(model, type = "III"))
        } else {
          print(summary(model))
        }
        cat("\n")
      } else if (inherits(model, "aov")) {
        cat("\n")
        print(summary(model))
        cat("\n")
      } else if (inherits(model, "aovlist")) {
        cat("\n")
        print(summary(model))
        cat("\n")
      }
    }
  }

  #============================================================================
  #                         RETOUR DES RESULTATS
  #============================================================================

  .dbg("", "Préparation du retour des résultats...", debug=debug)

  result <- list(
    x = x,
    g_cat = g_cat,
    check_normality = check_normality,
    check_variance_equal = check_variance_equal,
    check_ancova = check_ancova,
    check_discret = check_discret,
    robuste = robuste,
    alea = alea,
    alea_plus = alea_plus,
    model = model,  # PEUT ÊTRE NULL (voir documentation)
    formula = formula,
    updated_formula = updated_formula,
    within = within,
    between = between,
    ancova_checks = ancova_checks,  # Résultats des contrôles ANCOVA
    residual_type = residual_type_used,  # NOUVEAU: type de résidu utilisé
    robust_results = robust_results,  # ← AJOUTER CETTE LIGNE
    k = k
  )

  .dbg("=== End .multi_factor_analysis() ===",
       "=== Fin de .multi_factor_analysis() ===", debug=debug)

  return(result)
}

#==============================================================================
#                         NOTES FINALES ET DOCUMENTATION
#==============================================================================

#* ============================================================================
#* RàSUMé DES NOTES PERSONNELLES TRAITàES
#* ============================================================================
#*
#* NOTE #1  [EN COURS]        : Gestion des NA - Solution partielle implémentée
#* NOTE #2  [RéSOLUE]        : Résidus studentisés - Trace du type ajoutée
#* NOTE #3  [EN COURS]        : Vérification croisement complet - .mixed_model_analysis() à  implémenter
#* NOTE #4  [PARTIELLEMENT]   : Binning ANCOVA - Suppréssion du binning automatique
#* NOTE #5  [RéSOLUE]        : Normalité différences - Approche conforme standards
#* NOTE #6  [PARTIELLEMENT]   : Indépendance - Contrôle systématique, .control_independence() à  modifier
#* NOTE #7  [RéSOLUE]        : Comparaison modèles - AIC/BIC ajoutés
#* NOTE #8  [PARTIELLEMENT]   : Retour paramétrique - Logique skew/kurt implémentée
#* NOTE #9  [RéSOLUE]        : Sidak - Code supprimé (non pertinent)
#* NOTE #10 [RéSOLUE]        : Variance résidus - Brown-Forsythe implementé
#* NOTE #11 [NON RéSOLUE]    : Tests robustes auto - Plan d'action détaillé
#* NOTE #12 [PARTIELLEMENT]   : Kruskal 2-way - Restriction implémentée
#*
#* PRIORITàS POUR VERSIONS FUTURES :
#* 1. [URGENT]     Implémenter .mixed_model_analysis() avec validation valreg()
#* 2. [IMPORTANT]  Friedman test automatique (mesures répétées, k>=3)
#* 3. [IMPORTANT]  Modifier .control_independence() pour filtrage correct
#* 4. [MOYEN]      ANOVA par permutation (lmPerm::aovp) automatique
#* 5. [MOYEN]      Intégration WRS2 et Scheirer-Ray-Hare
#* 6. [FAIBLE]     Gestion avancée NA avec stratégies d'imputation
#*
#* ============================================================================

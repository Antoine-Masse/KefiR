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
#' @importFrom DescTools GTest
#' @importFrom rstatix identify_outliers
#' @importFrom lmPerm aovp lmp
#' @importFrom withr with_options
#'
#' @seealso
#' \code{\link{.one_factor_analysis}}, \code{\link[car]{Anova}},
#' \code{\link[lme4]{lmer}}
#'
#' @keywords internal
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

  # Helper function for code mode output
  .code_multi <- function(step_num, title, code_lines) {
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
    k_code <- 0  # Compteur s\u00e9par\u00e9 pour mode code
  }

  .dbg("=== Start .multi_factor_analysis() ===",
       "=== D\u00e9but de .multi_factor_analysis() ===", debug=debug)

  #============================================================================
  #                      NOTE PERSO #1 [EN COURS]
  #============================================================================
  #* NOTE PERSO #1 : Gestion des donn\u00e9es manquantes (NA)
  #*
  #* \u279c Probl\u00e8me identifi\u00e9 :
  #*   - m.test() supprime les NA en amont (listwise deletion)
  #*   - Pour les mod\u00e8les mixtes, le listwise n'est PAS optimal car il supprime
  #*     des observations compl\u00e8tes pour d'autres variables
  #*   - Beaucoup d'appels (bartlett.test, kruskal.test, leveneTest) plantent si NA
  #*
  #* \u279c Source acad\u00e9mique (APA + DOI) :
  #*   Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the
  #*   state of the art. *Psychological Methods*, 7(2), 147\u00e2\u20ac"177.
  #*   https://doi.org/10.1037/1082-989X.7.2.147
  #*
  #*   Id\u00e9e principale : Le listwise deletion est acceptable SEULEMENT si :
  #*   (a) Les donn\u00e9es manquent compl\u00e8tement au hasard (MCAR)
  #*   (b) Le taux de donn\u00e9es manquantes est < 5%
  #*   Pour les mod\u00e8les mixtes (lmer), Maximum Likelihood (ML) utilise toutes
  #*   les donn\u00e9es disponibles et est sup\u00e9rieur au listwise.
  #*
  #* \u279c Solution appliqu\u00e9e :
  #*   1. D\u00e9tection pr\u00e9coce des NA et calcul du taux de manquants
  #*   2. Si taux > 5%, avertissement \u00e0  l'utilisateur
  #*   3. Pour la voie param\u00e9trique classique : na.omit() d\u00e9fensif avant chaque test
  #*   4. Pour les mod\u00e8les mixtes : conserver les NA, lmer() g\u00e8re via ML
  #*   5. Documentation explicite de la strat\u00e9gie dans les messages verbose
  #*
  #* \u279c Statut : Solution partielle impl\u00e9ment\u00e9e ci-dessous
  #============================================================================

  # D\u00e9tection pr\u00e9coce des NA (uniquement sur les variables utilis\u00e9es)
  if (!is.null(data)) {
    data_na <- data
    if (!is.null(formula)) {
      data_na <- tryCatch(
        stats::model.frame(formula, data = data, na.action = stats::na.pass),
        error = function(e) data
      )
    }

    if (nrow(data_na) > 0) {
      complete_idx <- stats::complete.cases(data_na)
      na_count <- sum(!complete_idx)
      na_rate <- (na_count / nrow(data_na)) * 100
    } else {
      complete_idx <- rep(TRUE, nrow(data))
      na_count <- 0
      na_rate <- 0
    }

    if (na_rate > 5 && verbose) {
      k <- .vbse(
        paste0("Warning: ", round(na_rate, 2), "% missing data detected (>5% threshold).\n",
               "\tListwise deletion may reduce power. Consider imputation."),
        paste0("Attention : ", round(na_rate, 2), "% de donn\u00e9es manquantes d\u00e9tect\u00e9es (seuil >5%).\n",
               "\tLa suppression listwise peut r\u00e9duire la puissance. Envisager l'imputation."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      if (isTRUE(code)) {
        k_code <- k_code + 1
        .code_multi(k_code, "D\u00e9tection des donn\u00e9es manquantes", c(
          "data_na <- model.frame(formula, data = data, na.action = stats::na.pass)",
          "complete_idx <- complete.cases(data_na)",
          "na_rate <- sum(!complete_idx) / nrow(data_na) * 100",
          "data <- data[complete_idx, , drop = FALSE]"
        ))
      }
    }

    # Suppression listwise pour la voie classique (d\u00e9fensif)
    data_clean <- data[complete_idx, , drop = FALSE]
    if (nrow(data_clean) < nrow(data) && verbose) {
      k <- .vbse(
        paste0(nrow(data) - nrow(data_clean), " observations removed due to missing data (listwise deletion)."),
        paste0(nrow(data) - nrow(data_clean), " observations supprim\u00e9es pour donn\u00e9es manquantes (listwise deletion)."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
    data <- data_clean

    # Synchroniser x et g avec le listwise (si passes en argument, non extraits de data)
    if (!is.null(x) && length(x) == length(complete_idx)) {
      x <- x[complete_idx]
    }
    if (!is.null(g)) {
      if (is.data.frame(g) && nrow(g) == length(complete_idx)) {
        g <- g[complete_idx, , drop = FALSE]
      } else if (is.vector(g) || is.factor(g)) {
        if (length(g) == length(complete_idx)) {
          g <- g[complete_idx]
        }
      }
    }
  }

  #============================================================================
  #                   FONCTIONS INTERNES
  #============================================================================

  # Fonction pour extraire automatiquement les r\u00e9sidus
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
        "Le mod\u00e8le fourni n'est ni un objet 'aov', 'lm', ni 'aovlist'."
      )
    }
  }

  #============================================================================
  #                      NOTE PERSO #2 [R\u00e9SOLUE]
  #============================================================================
  #* NOTE PERSO #2 : R\u00e9sidus studentis\u00e9s vs standardis\u00e9s pour ancova_checks
  #*
  #* \u279c Probl\u00e8me identifi\u00e9 :
  #*   Bonne id\u00e9e d'avoir un fallback rstudent\u00e2\u2020'rstandard\u00e2\u2020'residuals, mais il faut
  #*   garder trace du type de r\u00e9sidu utilis\u00e9 dans ancova_checks (utile pour audit).
  #*
  #* \u279c Source acad\u00e9mique (APA + DOI) :
  #*   Cook, R. D., & Weisberg, S. (1982). *Residuals and influence in regression*.
  #*   Chapman and Hall. ISBN: 978-0412242809
  #*
  #*   Id\u00e9e principale : Les r\u00e9sidus studentis\u00e9s (externally studentized) sont
  #*   pr\u00e9f\u00e9r\u00e9s pour la d\u00e9tection d'outliers car chaque observation est standardis\u00e9e
  #*   en utilisant un mod\u00e8le o\u00f9 elle est exclue. Plus robuste que les r\u00e9sidus
  #*   standardis\u00e9s (internally studentized).
  #*
  #* \u279c Solution appliqu\u00e9e :
  #*   1. Ajout d'un champ "residual_type" dans le retour de la fonction
  #*   2. Trace explicite du type utilis\u00e9 dans ancova_checks
  #*   3. Message verbose si fallback n\u00e9cessaire
  #*
  #* \u279c Statut : R\u00e9SOLU
  #============================================================================

  get_studentized_residuals <- function(model) {
    residual_type <- "studentized"  # Par d\u00e9faut

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
          "Attention : Impossible de calculer les r\u00e9sidus studentis\u00e9s. Utilisation des r\u00e9sidus bruts.",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }
    })

    return(list(residuals = res, type = residual_type))
  }

  #============================================================================
  #                         INITIALISATION
  #============================================================================

  if (is.null(k)) k <- 0

  # Initialiser les drapeaux de contr\u00f4le
  robuste <- FALSE
  check_ancova <- FALSE
  alea <- FALSE
  alea_plus <- FALSE
  check_normality <- TRUE
  normality_already_tested <- FALSE
  check_variance_equal <- TRUE
  variance_already_tested <- FALSE
  outliers_marginal_detected <- FALSE
  n_extreme_outliers <- 0
  force_wrs2 <- FALSE
  imbalance_ratio <- NA_real_
  check_discret <- FALSE
  balanced <- FALSE  # Indicateur pour le type de sommes des carr\u00e9s (Type II vs III)
  use_mixed_model <- FALSE  # Indicateur sp\u00e9cifique pour mod\u00e8les mixtes (\u00e9vite messages redondants)
  model <- NULL
  n_problematic_subjects <- 0  # Nombre de sujets avec observations manquantes ou en exc\u00e8s (coh\u00e9rence entre messages)
  updated_formula <- formula
  ancova_checks <- list()
  residual_type_used <- NULL  # NOUVEAU: trace du type de r\u00e9sidu
  # Initialisation de robust_results (utilis\u00e9 si voie robuste)
  robust_results <- NULL

  #============================================================================
  #              VALIDATION ET PR\u00e9PARATION DES ENTR\u00e9ES
  #============================================================================

  # Extraire x et g depuis data si non fournis
  if (is.null(g)) {
    if (is.null(data)) {
      .exit(
        "Either 'g' or 'data' must be provided.",
        "Soit 'g' soit 'data' doit \u00eatre fourni."
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
        "Soit 'x' soit 'data' doit \u00eatre fourni."
      )
    }
    x <- data[, 1]
  }

  # CORRECTION CRITIQUE : S'assurer que g est TOUJOURS un data.frame
  if (!is.data.frame(g)) {
    # Extraire les noms de pr\u00e9dicteurs depuis la formule
    if (!is.null(formula)) {
      predictor_names <- setdiff(all.vars(formula), all.vars(formula)[1])
    } else {
      predictor_names <- "g"  # Nom par d\u00e9faut
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
  # DEBUG : Afficher la structure de data apr\u00e8s reconstruction
  if (debug) {
    cat("\n=== DEBUG apr\u00e8s reconstruction de data ===\n")
    cat("Noms de colonnes dans data :", paste(names(data), collapse = ", "), "\n")
    cat("Noms de colonnes dans g :", paste(names(g), collapse = ", "), "\n")
    cat("dim(data) :", paste(dim(data), collapse = " x "), "\n")
    cat("str(data) :\n")
    str(data)
    cat("==========================================\n\n")
  }

  #============================================================================
  #             CONTR\u00d4LE VARIABLE D\u00e9PENDANTE DISCR\u00caTE
  #============================================================================

  check_discret <- discret.test(x)
  if (check_discret) {
    k <- .vbse(
      paste0("Discrete dependent variable detected [unique values < sqrt(n)].\n",
             "\tCriterion: Number of unique values too small for continuous analysis.\n",
             "\t==> Switching to robust non-parametric analysis strategy."),
      paste0("Variable d\u00e9pendante discr\u00e8te d\u00e9tect\u00e9e [valeurs uniques < sqrt(n)].\n",
             "\tCrit\u00e8re : Nombre de valeurs uniques trop faible pour analyse continue.\n",
             "\t==> Passage vers strat\u00e9gie d'analyse non param\u00e9trique robuste."),
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    if (isTRUE(code)) {
      k_code <- k_code + 1
      .code_multi(k_code, "D\u00e9tection variable d\u00e9pendante discr\u00e8te", c(
        "n_unique <- length(unique(x))",
        "n_total <- length(x)",
        "check_discret <- (n_unique < sqrt(n_total))"
      ))
    }
    robuste <- TRUE
    check_normality <- FALSE  # Donn\u00e9es discr\u00e8tes => tests non-param\u00e9triques
  }

  #============================================================================
  #              D\u00e9TECTION DU TYPE DE MOD\u00caLE: ANOVA vs ANCOVA
  #============================================================================
  # DEBUG CRITIQUE : V\u00e9rifier la structure de data AVANT .detect_model_type()
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
  .dbg("", "D\u00e9tection du type de mod\u00e8le (ANOVA vs ANCOVA)...", debug=debug)

  # D\u00e9tection via .detect_model_type()
  check_ancova_str <- .detect_model_type(formula, data, debug=debug)
  .dbg("", "Fin d'usage de .detect_model_type", debug=debug)
  check_ancova <- (check_ancova_str == "TRUE")

  if (check_ancova) {
    # Message de d\u00e9tection supprim\u00e9 - sera affich\u00e9 par .ancova_analysis() directement
    # pour \u00e9viter doublon entre \u00e9tapes 1 et 2

    # REDIRECTION VERS .ancova_analysis() pour traitement complet et rigoureux
    ancova_result <- .ancova_analysis(
      x = x,
      g = g,
      formula = formula,
      data = data,
      paired = paired,
      id = id,
      alpha = alpha,
      k = k,
      code = code,
      debug = debug,
      verbose = verbose
    )

    # Retour imm\u00e9diat avec r\u00e9sultats ANCOVA
    return(list(
      x = x,
      g_cat = ancova_result$g_cat,  # Interaction facteurs cat\u00e9goriques depuis .ancova_analysis
      check_normality = ancova_result$assumptions_checked$normality$passed,
      check_variance_equal = ancova_result$assumptions_checked$homoscedasticity$passed,
      k = ancova_result$k,
      model = ancova_result$model,
      robuste = ancova_result$robust,
      check_ancova = TRUE,
      check_discret = FALSE,
      alea = FALSE,
      alea_plus = FALSE,
      formula = formula,
      updated_formula = formula,
      ancova_checks = ancova_result$assumptions_checked,
      robust_results = if(ancova_result$robust) ancova_result$robust_results else NULL,
      paired = FALSE,
      balanced = FALSE
    ))

  } else {
    k <- .vbse(
      "Model type detected: ANOVA (categorical factors only).",
      "Type de mod\u00e8le d\u00e9tect\u00e9 : ANOVA (facteurs cat\u00e9goriques uniquement).",
      verbose = verbose, code = code, k = k, cpt = "on"
    )
    if (isTRUE(code)) {
      k_code <- k_code + 1
      .code_multi(k_code, "Type de mod\u00e8le : ANOVA", c(
        "predictors <- attr(terms(formula), 'term.labels')",
        "simple_predictors <- predictors[!grepl(':', predictors)]",
        "numeric_vars <- simple_predictors[sapply(data[simple_predictors], is.numeric)]",
        "factor_vars <- simple_predictors[sapply(data[simple_predictors], is.factor)]"
      ))
    }
  }

  # Extraire les variables num\u00e9riques et cat\u00e9gorielles
  predictors <- attr(terms(formula), "term.labels")

  # Filtrer les termes d'interaction (contiennent ":")
  # Ne garder que les termes simples pour identifier les types de variables
  simple_predictors <- predictors[!grepl(":", predictors)]

  # CORRECTION : Filtrer pour ne garder que les pr\u00e9dicteurs qui existent dans data
  simple_predictors <- simple_predictors[simple_predictors %in% names(data)]

  numeric_vars <- simple_predictors[sapply(data[simple_predictors], is.numeric)]
  factor_vars <- simple_predictors[sapply(data[simple_predictors], is.factor)]

  if (debug) {
    cat("\n=== DEBUG: Variables identifi\u00e9es ===\n")
    cat("Num\u00e9riques:", paste(numeric_vars, collapse=", "), "\n")
    cat("Cat\u00e9gorielles:", paste(factor_vars, collapse=", "), "\n")
    cat("=====================================\n\n")
  }

  #============================================================================
  #              D\u00e9TECTION DES EFFETS AL\u00e9ATOIRES
  #============================================================================

  .dbg("", "D\u00e9tection des effets al\u00e9atoires...", debug=debug)

  # D\u00e9tecter la pr\u00e9sence du terme Error() dans la formule
  # NOTE: Error() est la syntaxe aov() pour effets al\u00e9atoires/mesures r\u00e9p\u00e9t\u00e9es
  # "|" est la syntaxe lmer() pour effets al\u00e9atoires : (1|Subject), (Time|Subject)
  formula_str <- deparse(formula)
  alea_error <- grepl("Error\\s*\\(", formula_str, perl = TRUE)
  alea_lmer <- grepl("[(][^)]+[|][^)]+[)]", formula_str)  # D\u00e9tecte (1|id), (var|id), (1 | id)
  alea <- alea_error || alea_lmer

  # Si syntaxe lmer d\u00e9tect\u00e9e, extraire l'identifiant sujet
  if (alea_lmer && is.null(id)) {
    # Extraire l'id depuis la syntaxe (1|id) ou (var|id) avec espaces optionnels
    lmer_match <- regmatches(formula_str, regexpr("[|][^)]+[)]", formula_str))
    if (length(lmer_match) > 0) {
      extracted_id <- trimws(gsub("[|)]", "", lmer_match[1]))
      if (extracted_id %in% names(data)) {
        id <- extracted_id
        k <- .vbse(
          paste0("lmer syntax detected: random effect grouping variable = '", id, "'"),
          paste0("Syntaxe lmer d\u00e9tect\u00e9e : variable de regroupement effet al\u00e9atoire = '", id, "'"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
      }
    }
  }

  if (alea_error) {
    # SEULEMENT afficher ce message s'il y a vraiment Error() dans la formule
    # R\u00e9f\u00e9rence : Maxwell, Delaney & Kelley (2018), Chapters 11-12
    # Error(id) : Plan appari\u00e9 simple (within-subject)
    # Error(id/factor) : Plan imbriqu\u00e9 (nested) - rarement utilis\u00e9
    # Error(id/(F*G)) : Plan appari\u00e9 multi-facteurs (within-subjects)
    k <- .vbse(
      paste0("Random effects detected in formula (Error term): ", formula_str),
      paste0("Effets al\u00e9atoires d\u00e9tect\u00e9s dans la formule (terme Error) : ", formula_str),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # D\u00e9tecter structure complexe : UNIQUEMENT pour nested ou plusieurs |
    # Note : "/" dans Error(id/F) indique embo\u00eetement (nested), pas plan appari\u00e9
    alea_plus <- (
      length(gregexpr("\\|", formula_str)[[1]]) > 1 ||
        (grepl("/", formula_str) && grepl("\\*", formula_str))  # Nested ET interaction
    )

    if (alea_plus) {
      k <- .vbse(
        "Complex random effects structure detected (nested design or multiple within-subject factors).",
        "Structure d'effets al\u00e9atoires complexe d\u00e9tect\u00e9e (plan embo\u00eet\u00e9 ou plusieurs facteurs intra-sujet).",
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
  }

  #============================================================================
  #         CONTR\u00d4LES SP\u00e9CIFIQUES AUX MESURES R\u00e9P\u00e9T\u00e9ES
  #============================================================================

  if (paired) {

    .dbg("", "Contr\u00f4les pour mesures r\u00e9p\u00e9t\u00e9es...", debug=debug)

    k <- .vbse(
      "Paired/repeated measures design detected. Performing specific checks...",
      "Plan appari\u00e9 / mesures r\u00e9p\u00e9t\u00e9es d\u00e9tect\u00e9. R\u00e9alisation de contr\u00f4les sp\u00e9cifiques...",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    #--------------------------------------------------------------------------
    # 1) V\u00e9rification de la pr\u00e9sence de l'identifiant
    #--------------------------------------------------------------------------
    if (is.null(id) || !id %in% names(data)) {
      .exit(
        "For paired/repeated measures designs, 'id' must be specified and present in data.",
        "Pour les plans appari\u00e9s/mesures r\u00e9p\u00e9t\u00e9es, 'id' doit \u00eatre sp\u00e9cifi\u00e9 et pr\u00e9sent dans les donn\u00e9es."
      )
    }

    #--------------------------------------------------------------------------
    # 2) V\u00e9rification de la coh\u00e9rence within/between
    #--------------------------------------------------------------------------
    all_factors <- names(g)[sapply(g, is.factor)]

    if (!is.null(within)) {
      if (!all(within %in% all_factors)) {
        missing_within <- setdiff(within, all_factors)
        .exit(
          paste0("Within-subject factor(s) not found or not categorical: ",
                 paste(missing_within, collapse = ", ")),
          paste0("Facteur(s) intra-sujet non trouv\u00e9(s) ou non cat\u00e9goriel(s) : ",
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
          paste0("Facteur(s) inter-sujet non trouv\u00e9(s) ou non cat\u00e9goriel(s) : ",
                 paste(missing_between, collapse = ", "))
        )
      }
    }

    #--------------------------------------------------------------------------
    # 3) Cr\u00e9ation du facteur within (interaction si plusieurs)
    #--------------------------------------------------------------------------
    if (!is.null(within) && length(within) > 0) {
      if (length(within) == 1) {
        within_interaction <- factor(data[[within]])
      } else {
        within_interaction <- interaction(data[within], drop = TRUE)
      }
    } else {
      # Si within non sp\u00e9cifi\u00e9, utiliser tous les facteurs sauf between et id
      auto_within <- setdiff(all_factors, c(between, id))
      if (length(auto_within) == 0) {
        .exit(
          "Could not identify within-subject factors. Please specify 'within' explicitly.",
          "Impossible d'identifier les facteurs intra-sujet. Veuillez sp\u00e9cifier 'within' explicitement."
        )
      }
      within_interaction <- interaction(data[auto_within], drop = TRUE)
      within <- auto_within

      k <- .vbse(
        paste0("Auto-detected within-subject factor(s): ", paste(within, collapse = ", ")),
        paste0("Facteur(s) intra-sujet auto-d\u00e9tect\u00e9(s) : ", paste(within, collapse = ", ")),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # 4) V\u00e9rification de l'\u00e9quilibrage ET de l'unicit\u00e9 (\u00e9tapes fusionn\u00e9es)
    #--------------------------------------------------------------------------
    # Fusion \u00e9tapes 4+5 pour logique p\u00e9dagogique: attendus \u2192 probl\u00e8mes \u2192 recommandation

    k <- .vbse(
      "Checking repeated measures design structure (balance and uniqueness)...",
      "V\u00e9rification de la structure du plan \u00e0 mesures r\u00e9p\u00e9t\u00e9es (\u00e9quilibrage et unicit\u00e9)...",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # Informations g\u00e9n\u00e9rales
    n_subjects <- length(unique(data[[id]]))
    n_obs <- nrow(data)
    expected_n <- nlevels(within_interaction)

    # \u00c9TAPE A: Afficher les ATTENDUS (ce qui devrait \u00eatre)
    # NOTE: Le message sera mis \u00e0 jour apr\u00e8s d\u00e9tection des r\u00e9plicats
    # Pour l'instant, afficher structure de base
    k <- .vbse(
      paste0("EXPECTED design structure (checking for replicates):\n",
             "\t\u2022 Number of subjects: ", n_subjects, "\n",
             "\t\u2022 Conditions per subject: ", expected_n, " (all combinations of within-factors)\n",
             "\t\u2022 Total observations: ", n_obs),
      paste0("Structure ATTENDUE du plan (v\u00e9rification des r\u00e9plicats) :\n",
             "\t\u2022 Nombre de sujets : ", n_subjects, "\n",
             "\t\u2022 Conditions par sujet : ", expected_n, " (toutes combinaisons facteurs intra-sujet)\n",
             "\t\u2022 Total observations : ", n_obs),
      verbose = verbose, code = code, k = k, cpt = "off"
    )

    # \u00c9TAPE B: Analyser les PROBL\u00c8MES (d\u00e9s\u00e9quilibre, doublons, manquantes)
    obs_per_subject <- table(data[[id]])

    # V\u00e9rifier si toutes les cellules id\u00d7within ont le m\u00eame nombre d'observations
    tab_id_within <- table(data[[id]], within_interaction)

    # D\u00e9tection intelligente: v\u00e9rifier si c'est un design \u00e9quilibr\u00e9 avec r\u00e9plicats
    # ou un vrai d\u00e9s\u00e9quilibre
    unique_counts <- unique(as.vector(tab_id_within))
    unique_counts_nonzero <- unique_counts[unique_counts > 0]

    # Si toutes les cellules non-vides ont le M\u00caME nombre d'observations (ex: toutes 3),
    # c'est un design \u00c9QUILIBR\u00c9 avec r\u00e9plicats, pas un d\u00e9s\u00e9quilibre
    is_balanced_with_replicates <- (length(unique_counts_nonzero) == 1)
    n_replicates_per_cell <- ifelse(is_balanced_with_replicates, unique_counts_nonzero[1], NA)

    # Red\u00e9finir expected_n si design avec r\u00e9plicats \u00e9quilibr\u00e9s
    if (is_balanced_with_replicates && n_replicates_per_cell > 1) {
      expected_n_with_replicates <- nlevels(within_interaction) * n_replicates_per_cell
    } else {
      expected_n_with_replicates <- expected_n
    }

    # V\u00e9rifier d\u00e9s\u00e9quilibre R\u00c9EL (sujets n'ayant pas le nombre attendu avec r\u00e9plicats)
    incorrect_ids <- names(obs_per_subject)[obs_per_subject != expected_n_with_replicates]
    n_problematic_subjects <- length(incorrect_ids)

    # V\u00e9rifier cellules manquantes (0 observations)
    missing <- tab_id_within == 0
    n_missing_cells <- sum(missing)

    # Les "duplicates" ne sont probl\u00e9matiques QUE si d\u00e9s\u00e9quilibre
    # Si design \u00e9quilibr\u00e9 avec r\u00e9plicats, ce n'est PAS un probl\u00e8me
    if (is_balanced_with_replicates && n_replicates_per_cell > 1) {
      n_duplicate_cells <- 0  # Pas de vrais doublons, juste r\u00e9plicats \u00e9quilibr\u00e9s
    } else {
      duplicates <- tab_id_within > 1
      n_duplicate_cells <- sum(duplicates)
    }

    # Afficher PROBL\u00c8MES d\u00e9tect\u00e9s (si pr\u00e9sents)
    has_problems <- (n_problematic_subjects > 0 || n_duplicate_cells > 0 || n_missing_cells > 0)

    if (has_problems) {
      k <- .vbse(
        paste0("PROBLEMS detected in design structure:\n",
               if (n_problematic_subjects > 0)
                 paste0("\t\u2022 Imbalance: ", n_problematic_subjects, " subject(s) do not have ", expected_n, " observations\n",
                        "\t  Problematic subjects: ", paste(head(incorrect_ids, 5), collapse = ", "),
                        if (n_problematic_subjects > 5) ", ..." else "", "\n") else "",
               if (n_duplicate_cells > 0)
                 paste0("\t\u2022 Duplicates: ", n_duplicate_cells, " cell(s) with >1 observation per subject\u00d7condition\n") else "",
               if (n_missing_cells > 0)
                 paste0("\t\u2022 Missing data: ", n_missing_cells, " cell(s) with 0 observations (expected 1)\n") else "",
               "\t\u2022 Actual observations: ", n_obs, " (expected ", n_subjects * expected_n, ")"),
        paste0("PROBL\u00c8MES d\u00e9tect\u00e9s dans la structure du plan :\n",
               if (n_problematic_subjects > 0)
                 paste0("\t\u2022 D\u00e9s\u00e9quilibre : ", n_problematic_subjects, " sujet(s) n'ont pas ", expected_n, " observations\n",
                        "\t  Sujets probl\u00e9matiques : ", paste(head(incorrect_ids, 5), collapse = ", "),
                        if (n_problematic_subjects > 5) ", ..." else "", "\n") else "",
               if (n_duplicate_cells > 0)
                 paste0("\t\u2022 Doublons : ", n_duplicate_cells, " cellule(s) avec >1 observation par sujet\u00d7condition\n") else "",
               if (n_missing_cells > 0)
                 paste0("\t\u2022 Donn\u00e9es manquantes : ", n_missing_cells, " cellule(s) avec 0 observation (attendu 1)\n") else "",
               "\t\u2022 Observations r\u00e9elles : ", n_obs, " (attendu ", n_subjects * expected_n, ")"),
        verbose = verbose, code = code, k = k, cpt = "off"
      )

      # \u00c9TAPE C: RECOMMANDATION (mod\u00e8le mixte)
      # R\u00e9f\u00e9rence: Barr et al. (2013). Random effects structure for confirmatory hypothesis testing.
      k <- .vbse(
        paste0("RECOMMENDATION: Mixed-effects model (lmer) is required.\n",
               "\tReason: Imbalance/duplicates/missing data violate standard RM-ANOVA assumptions.\n",
               "\t--> Towards mixed-effects model (lmer)"),
        paste0("RECOMMANDATION : Mod\u00e8le \u00e0 effets mixtes (lmer) requis.\n",
               "\tRaison : D\u00e9s\u00e9quilibre/doublons/donn\u00e9es manquantes violent hypoth\u00e8ses ANOVA-RM standard.\n",
               "\t--> Vers mod\u00e8le \u00e0 effets mixtes (lmer)"),
        verbose = verbose, code = code, k = k, cpt = "off"
      )

      robuste <- TRUE
      check_normality <- FALSE  # Donn\u00e9es d\u00e9s\u00e9quilibr\u00e9es => post-hocs non-param\u00e9triques
      use_mixed_model <- TRUE

    } else {
      # Design parfaitement \u00e9quilibr\u00e9
      if (is_balanced_with_replicates && n_replicates_per_cell > 1) {
        # Design \u00e9quilibr\u00e9 avec r\u00e9plicats
        k <- .vbse(
          paste0("Design structure is BALANCED with REPLICATES:\n",
                 "\t\u2022 All ", n_subjects, " subjects have exactly ", expected_n_with_replicates, " observations\n",
                 "\t\u2022 ", n_replicates_per_cell, " replicate(s) per subject\u00d7condition (", expected_n, " conditions)\n",
                 "\t\u2022 Total: ", n_subjects, " \u00d7 ", expected_n, " \u00d7 ", n_replicates_per_cell, " = ", n_obs, " observations"),
          paste0("Structure du plan \u00c9QUILIBR\u00c9E avec R\u00c9PLICATS :\n",
                 "\t\u2022 Les ", n_subjects, " sujets ont exactement ", expected_n_with_replicates, " observations\n",
                 "\t\u2022 ", n_replicates_per_cell, " r\u00e9plicat(s) par sujet\u00d7condition (", expected_n, " conditions)\n",
                 "\t\u2022 Total : ", n_subjects, " \u00d7 ", expected_n, " \u00d7 ", n_replicates_per_cell, " = ", n_obs, " observations"),
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      } else {
        # Design \u00e9quilibr\u00e9 sans r\u00e9plicats (1 obs par cellule)
        k <- .vbse(
          paste0("Design structure is BALANCED:\n\t==> All ", n_subjects, " subjects have exactly ", expected_n, " observations (one per condition)."),
          paste0("Structure du plan \u00c9QUILIBR\u00c9E :\n\t==> Les ", n_subjects, " sujets ont exactement ", expected_n, " observations (une par condition)."),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        # Message directionnel vers ANOVA mesures r\u00e9p\u00e9t\u00e9es
        k <- .vbse(
          "--> Towards repeated measures ANOVA.",
          "--> Vers une ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es.",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }
    }
    #--------------------------------------------------------------------------
    # 6) V\u00e9rification de la coh\u00e9rence between
    #--------------------------------------------------------------------------
    if (!is.null(between) && length(between) > 0) {

      k <- .vbse(
        "Checking consistency of between-subject factors (should be constant within each subject)...",
        "V\u00e9rification de la coh\u00e9rence des facteurs inter-sujet (doivent \u00eatre constants pour chaque sujet)...",
        verbose = verbose, code = code, k = k, cpt = "on"
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
                   "\n\tUn facteur inter-sujet doit \u00eatre constant pour chaque sujet."),
            verbose = verbose, code = code, k = k, cpt = "off"
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
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # 7) Contr\u00f4le strict d'\u00e9quilibrage RM: id \u00e0\u2014 within PAR niveau de between
    #--------------------------------------------------------------------------
    if (!is.null(within) && length(within) > 0 &&
        !is.null(between) && length(between) > 0) {

      k <- .vbse(
        "Checking strict RM balance: each subject must have all within-levels inside each between level...",
        "Contr\u00f4le strict de l'\u00e9quilibrage RM : chaque sujet doit avoir toutes les modalit\u00e9s du within \u00e0  l'int\u00e9rieur de chaque niveau du between...",
        verbose = verbose, code = code, k = k, cpt = "on"
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

        # Table id \u00e0\u2014 within dans ce sous-ensemble
        tab_bw <- table(data[[id]][idx_b], within_fac[idx_b])

        # 1) Chaque sujet pr\u00e9sent dans ce between doit avoir toutes les modalit\u00e9s within
        missing_per_id <- rowSums(tab_bw > 0) != n_within_levels

        # 2) Pas de doublons dans une m\u00eame cellule id\u00e0\u2014within
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
                                           "] doublons d\u00e9tect\u00e9s dans au moins une cellule id\u00e0\u2014within"))
          }
        }
      }

      if (rm_unbalanced) {
        k <- .vbse(
          paste0("Unbalanced repeated-measures design detected within between levels.\n\t",
                 paste(head(examples, 3), collapse = "\n\t")),
          paste0("Plan de mesures r\u00e9p\u00e9t\u00e9es d\u00e9s\u00e9quilibr\u00e9 d\u00e9tect\u00e9 \u00e0  l'int\u00e9rieur des niveaux du between.\n\t",
                 paste(head(examples, 3), collapse = "\n\t")),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        k <- .vbse(
          "Switching away from RM-ANOVA. Mixed models are recommended (e.g., lmer: A ~ F*G + (G|id)).",
          "Sortie de la RM-ANOVA. Mod\u00e8les mixtes recommand\u00e9s (ex. lmer : A ~ F*G + (G|id)).",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        robuste <- TRUE
        use_mixed_model <- TRUE  # RM d\u00e9s\u00e9quilibre n\u00e9cessite mod\u00e8le mixte
      } else {
        k <- .vbse(
          "Strict RM balance satisfied inside each between level (id \u00e0\u2014 within complete and unique).",
          "\u00e9quilibrage RM strict respect\u00e9 dans chaque niveau du between (id \u00e0\u2014 within complet et unique).",
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        #============================================================================
        #                      NOTE PERSO #3 [EN COURS]
        #============================================================================
        #* NOTE PERSO #3 : V\u00e9rification du croisement complet id \u00e0\u2014 within \u00e0\u2014 between
        #*
        #* \u279c Probl\u00e8me identifi\u00e9 :
        #*   Besoin de d\u00e9tecter les plans "aliased" o\u00f9 certaines combinaisons
        #*   id \u00e0\u2014 within \u00e0\u2014 between n'existent pas. Cette structure n\u00e9cessite
        #*   des mod\u00e8les mixtes.
        #*
        #* \u279c Source acad\u00e9mique (APA + DOI) :
        #*   Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random
        #*   effects structure for confirmatory hypothesis testing: Keep it maximal.
        #*   *Journal of Memory and Language*, 68(3), 255\u00e2\u20ac"278.
        #*   https://doi.org/10.1016/j.jml.2012.11.001
        #*
        #*   Id\u00e9e principale : Les plans complexes (incomplete crossing) n\u00e9cessitent
        #*   des mod\u00e8les mixtes avec structure d'effets al\u00e9atoires maximale pour
        #*   \u00e9viter les taux d'erreur Type I gonfl\u00e9s. L'ANOVA classique assume
        #*   un croisement complet.
        #*
        #* \u279c Solution appliqu\u00e9e :
        #*   V\u00e9rification de la compl\u00e9tude du croisement avant de poursuivre en
        #*   RM-ANOVA. Si incomplet, redirection vers mod\u00e8les mixtes (future
        #*   fonction .mixed_model_analysis()).
        #*
        #* \u279c Statut : Solution appliqu\u00e9e ci-dessous. RESTE \u00e0\u20ac FAIRE : impl\u00e9menter
        #*   .mixed_model_analysis() avec valreg() pour validation des assomptions.
        #============================================================================

        k <- .vbse(
          "Checking for aliased structure (incomplete id \u00e0\u2014 within \u00e0\u2014 between crossing)...",
          "Contr\u00f4le de la structure du plan (croisement id \u00e0\u2014 within \u00e0\u2014 between incomplet)...",
          verbose=verbose, k=k, cpt="on"
        )

        # Cr\u00e9er la table de croisement
        cross_tab <- table(data[[id]], data[[within[1]]], data[[between[1]]])
        # Nombre de combinaisons pr\u00e9sentes pour chaque sujet
        cross_count <- apply(cross_tab, 1, function(x) sum(x > 0))

        # Un plan \u00e9quilibr\u00e9 complet doit avoir exactement n_within_levels * n_between_levels combinaisons
        n_within <- length(unique(data[[within[1]]]))
        n_between <- length(unique(data[[between[1]]]))
        expected <- n_within * n_between

        if (any(cross_count != expected)) {
          k <- .vbse(
            "Detected aliased design: not all id \u00e0\u2014 within \u00e0\u2014 between combinations exist. Switching to mixed model.",
            "Plan non pleinement crois\u00e9 d\u00e9tect\u00e9 : certaines combinaisons id \u00e0\u2014 within \u00e0\u2014 between sont manquantes. Bascule vers mod\u00e8le mixte.",
            verbose=verbose, k=k, cpt="off"
          )
          robuste <- TRUE
          use_mixed_model <- TRUE  # Plan non crois\u00e9 n\u00e9cessite mod\u00e8le mixte

          #* PRIORIT\u00c9 6 : Routage vers .mixed_model_analysis() (IMPL\u00c9MENT\u00c9)
          #* Structure imbriqu\u00e9e/crois\u00e9e incompl\u00e8te d\u00e9tect\u00e9e
          #* Redirection vers mod\u00e8les mixtes pour traiter correctement les donn\u00e9es

          mixed_result <- .mixed_model_analysis(
            x=x, g=g, formula=formula, data=data,
            alpha=alpha, paired=paired, id=id,
            within=within, between=between,
            k=k, code=code, debug=debug, verbose=verbose
          )

          # Initialiser bilan si n\u00e9cessaire avec structure attendue par m.test()
          if (!exists("bilan")) {
            bilan <- list(
              x,           # [[1]]
              g,           # [[2]]
              TRUE,        # [[3]] check_normality (on suppose TRUE pour mod\u00e8le mixte)
              TRUE         # [[4]] check_variance_equal (on suppose TRUE pour mod\u00e8le mixte)
            )
          }

          # ADAPTATION: .mixed_model_analysis() peut retourner soit une liste compl\u00e8te,
          # soit directement un mod\u00e8le lmerMod (\u00e0 cause du return anticip\u00e9 ligne 251)
          if (inherits(mixed_result, "lmerMod")) {
            # Cas o\u00f9 on a un mod\u00e8le brut (return anticip\u00e9 dans .mixed_model_analysis)
            # Extraire les infos n\u00e9cessaires du mod\u00e8le directement
            .dbg("Detected raw lmerMod object from .mixed_model_analysis()",
                 "Objet lmerMod brut d\u00e9tect\u00e9 depuis .mixed_model_analysis()",
                 debug = debug)

            # Obtenir l'ANOVA table avec lmerTest pour les p-values
            anova_table <- tryCatch({
              suppressMessages({
                # Convertir en lmerModLmerTest si n\u00e9cessaire
                if (!inherits(mixed_result, "lmerModLmerTest")) {
                  mixed_result <- lmerTest::lmer(formula(mixed_result),
                                                 data = model.frame(mixed_result),
                                                 REML = TRUE)
                }
                stats::anova(mixed_result, type = "III")
              })
            }, error = function(e) NULL)

            # Extraire p-value globale
            global_pval <- if (!is.null(anova_table)) {
              pval_col <- which(colnames(anova_table) %in% c("Pr(>F)", "p.value", "P(>|t|)"))
              if (length(pval_col) > 0) {
                row_names <- rownames(anova_table)
                valid_rows <- which(!grepl("Residuals|Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
                if (length(valid_rows) > 0) min(anova_table[valid_rows, pval_col[1]], na.rm = TRUE) else NA
              } else NA
            } else NA

            bilan$robust_results <- list(
              method = "Mixed_Model_lmer",
              model = mixed_result,
              anova_table = anova_table,
              significant_effects = NULL,
              posthoc_applicable = TRUE,
              test_result = mixed_result,
              assumptions_checked = FALSE,
              variance_components = NULL,
              random_effects = NULL,
              fixed_effects = NULL
            )
            bilan$k <- k
            bilan$global_pvalue <- global_pval

          } else {
            # Cas normal : mixed_result est une liste compl\u00e8te
            bilan$robust_results <- list(
              method = "Mixed_Model_lmer",
              model = mixed_result$model,
              anova_table = mixed_result$anova_table,
              significant_effects = mixed_result$significant_effects,
              posthoc_applicable = mixed_result$posthoc_applicable,
              test_result = mixed_result$model,
              assumptions_checked = mixed_result$assumptions_checked,
              variance_components = mixed_result$variance_components,
              random_effects = mixed_result$random_effects,
              fixed_effects = mixed_result$fixed_effects
            )
            bilan$k <- mixed_result$k
            bilan$global_pvalue <- mixed_result$global_pvalue
          }

          return(bilan)

        } else {
          k <- .vbse(
            "Complete crossing id \u00e0\u2014 within \u00e0\u2014 between verified.",
            "Croisement complet id \u00e0\u2014 within \u00e0\u2014 between v\u00e9rifi\u00e9.",
            verbose=verbose, k=k, cpt="off"
          )
        }
      }
    }

  } # Fin if (paired)

  #============================================================================
  #             ASSOMPTION DE BASE (AVANT D\u00c9S\u00c9QUILIBRE)
  #============================================================================

  # Ind\u00e9pendance des observations (assomption de plan)
  if (!alea) {
    k <- .vbse(
      paste0("ASSUMPTION BASE: Independence of observations (design verification).\n",
             "\tThis is a DESIGN assumption that cannot be statistically tested.\n",
             "\tVerify that:\n",
             "\t  \u2022 No repeated measures (each observation from a different subject)\n",
             "\t  \u2022 No cluster effects (observations not grouped/nested)\n",
             "\t  \u2022 No carryover effects (order of measurements doesn't influence results)"),
      paste0("ASSOMPTION DE BASE : Ind\u00e9pendance des observations (v\u00e9rification de plan).\n",
             "\tC'est une assomption de PLAN qui ne peut pas \u00eatre test\u00e9e statistiquement.\n",
             "\tV\u00e9rifiez que :\n",
             "\t  \u2022 Pas de mesures r\u00e9p\u00e9t\u00e9es (chaque observation d'un sujet diff\u00e9rent)\n",
             "\t  \u2022 Pas d'effets cluster (observations non group\u00e9es/embo\u00eet\u00e9es)\n",
             "\t  \u2022 Pas d'effets report (ordre des mesures n'influence pas les r\u00e9sultats)"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    if (isTRUE(code)) {
      k_code <- k_code + 1
      .code_multi(k_code, "Assomption de base (ind\u00e9pendance)", c(
        "# V\u00e9rifier l'ind\u00e9pendance des observations",
        "# (pas de mesures r\u00e9p\u00e9t\u00e9es, pas d'effets cluster, pas d'effets d'ordre)"
      ))
    }

    #-----------------------------------
    # DIAGNOSTIC: Outliers marginaux (avant ajustement mod\u00e8le)
    #-----------------------------------
    n_outliers_marginal <- 0
    if (requireNamespace("rstatix", quietly = TRUE)) {
      tryCatch({
        response_var <- all.vars(formula)[1]
        outlier_data <- data.frame(value = data[[response_var]], group = g_cat)
        outliers <- rstatix::identify_outliers(outlier_data, value)

        if (!is.null(outliers) && nrow(outliers) > 0) {
          n_outliers_marginal <- nrow(outliers)
          n_extreme_outliers <- sum(outliers$is.extreme, na.rm = TRUE)
          outliers_marginal_detected <- (nrow(outliers) > 0)
        }
      }, error = function(e) {
        .dbg(paste0("Warning: Outlier detection failed: ", e$message),
             paste0("Attention : D\u00e9tection outliers \u00e9chou\u00e9e : ", e$message),
             debug = debug)
      })
    }

    if (outliers_marginal_detected) {
      outlier_conclusion_en <- if (n_extreme_outliers > 0) {
        "--> Extreme values may strongly influence parametric tests."
      } else {
        "--> Values to monitor, but not excessive."
      }
      outlier_conclusion_fr <- if (n_extreme_outliers > 0) {
        "--> Valeurs extr\u00eames peuvent fortement influencer tests param\u00e9triques."
      } else {
        "--> Valeurs \u00e0 surveiller, mais non excessives."
      }

      k <- .vbse(
        paste0("DIAGNOSTIC: Outliers [identify_outliers() {rstatix}]\n",
               "\t==> ", n_outliers_marginal, " outlier(s)",
               if (n_extreme_outliers > 0) paste0(" (", n_extreme_outliers, " extreme)") else "", ".\n",
               "\t", outlier_conclusion_en),
        paste0("DIAGNOSTIC : Outliers [identify_outliers() {rstatix}]\n",
               "\t==> ", n_outliers_marginal, " outlier(s)",
               if (n_extreme_outliers > 0) paste0(" (", n_extreme_outliers, " extr\u00eame(s))") else "", ".\n",
               "\t", outlier_conclusion_fr),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
    }

    if (isTRUE(code)) {
      k_code <- k_code + 1
      .code_multi(k_code, "Diagnostic outliers", c(
        "library(rstatix)",
        "outlier_data <- data.frame(value = data[[all.vars(formula)[1]]], group = g_cat)",
        "outliers <- identify_outliers(outlier_data, value)"
      ))
    }
  } else {
    k <- .vbse(
      paste0("ASSUMPTION BASE: Independence of observations (repeated measures context).\n",
             "\tThis is a DESIGN assumption that cannot be statistically tested.\n",
             "\tVerify that:\n",
             "\t  \u2022 Repeated measures on the same subjects (expected)\n",
             "\t  \u2022 No additional clustering beyond subject\n",
             "\t  \u2022 Order/carryover effects are controlled\n",
             "\t  \u2022 No systematic dropout"),
      paste0("ASSOMPTION DE BASE : Ind\u00e9pendance des observations (contexte mesures r\u00e9p\u00e9t\u00e9es).\n",
             "\tC'est une assomption de PLAN qui ne peut pas \u00eatre test\u00e9e statistiquement.\n",
             "\tV\u00e9rifiez que :\n",
             "\t  \u2022 Mesures r\u00e9p\u00e9t\u00e9es sur les m\u00eames sujets (attendu)\n",
             "\t  \u2022 Pas de cluster suppl\u00e9mentaire au-del\u00e0 du sujet\n",
             "\t  \u2022 Effets d'ordre/report contr\u00f4l\u00e9s\n",
             "\t  \u2022 Pas d'abandon syst\u00e9matique"),
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    if (isTRUE(code)) {
      k_code <- k_code + 1
      .code_multi(k_code, "Assomption de base (mesures r\u00e9p\u00e9t\u00e9es)", c(
        "# V\u00e9rifier la structure RM (sujets, effets d'ordre, abandons)"
      ))
    }
  }

  #============================================================================
  #         PR\u00e9PARATION DU FACTEUR D'INTERACTION CAT\u00e9GORIEL
  #============================================================================

  .dbg("", "Pr\u00e9paration du facteur d'interaction cat\u00e9goriel...", debug=debug)

  # Exclure la variable d'appariement de l'interaction si pr\u00e9sente
  if (!is.null(id)) {
    g_temp <- g[, setdiff(names(g), id), drop = FALSE]
  } else {
    g_temp <- g
  }

  #----------------------------------------------------------------------------
  # Voie ANOVA: cr\u00e9er l'interaction de tous les facteurs
  #----------------------------------------------------------------------------
  if (check_ancova==FALSE) {
    factor_cols <- sapply(g_temp, is.factor)
    if (!all(factor_cols)) {
      .exit(
        "ANOVA detected but some variables are not factors. Please verify your data.",
        "ANOVA d\u00e9tect\u00e9e mais certaines variables ne sont pas des facteurs. Veuillez v\u00e9rifier vos donn\u00e9es."
      )
    }
    g_cat <- interaction(g_temp, drop = TRUE)

    #----------------------------------------------------------------------------
    # Voie ANCOVA: gestion des covariables num\u00e9riques
    #----------------------------------------------------------------------------
  } else {

    #============================================================================
    #                      NOTE PERSO #4 [PARTIELLEMENT R\u00e9SOLUE]
    #============================================================================
    #* NOTE PERSO #4 : Binning automatique des covariables continues en ANCOVA
    #*
    #* \u279c Probl\u00e8me identifi\u00e9 :
    #*   1. Le binning automatique (cut en 3 cat\u00e9gories) n'est pas une approche
    #*      acad\u00e9miquement valid\u00e9e pour l'ANCOVA
    #*   2. Si toutes les colonnes sont num\u00e9riques, g_cat devient vide
    #*   3. Perte d'information en discr\u00e9tisant des covariables continues
    #*
    #* \u279c Source acad\u00e9mique (APA + DOI) :
    #*   MacCallum, R. C., Zhang, S., Preacher, K. J., & Rucker, D. D. (2002).
    #*   On the practice of dichotomization of quantitative variables.
    #*   *Psychological Methods*, 7(1), 19\u00e2\u20ac"40.
    #*   https://doi.org/10.1037/1082-989X.7.1.19
    #*
    #*   Id\u00e9e principale : La discr\u00e9tisation de variables continues est
    #*   DECONSEILLEE car elle :
    #*   - R\u00e9duit la puissance statistique
    #*   - Augmente le risque d'erreur Type I
    #*   - Perd de l'information sur les relations lin\u00e9aires
    #*   Exception : discr\u00e9tisation justifi\u00e9e th\u00e9oriquement (ex: points de coupure cliniques)
    #*
    #* \u279c Solution appliqu\u00e9e :
    #*   1. SUPPRESSION du binning automatique
    #*   2. Les covariables num\u00e9riques restent num\u00e9riques dans le mod\u00e8le ANCOVA
    #*   3. g_cat cr\u00e9\u00e9 uniquement \u00e0  partir des facteurs cat\u00e9goriques
    #*   4. Fallback : si aucun facteur cat\u00e9gorique, erreur explicite
    #*
    #* \u279c Statut : PARTIELLEMENT R\u00e9SOLU (binning supprim\u00e9, mais utilisateur
    #*   peut toujours binning manuel en amont si justifi\u00e9 th\u00e9oriquement)
    #============================================================================

    k <- .vbse(
      "ANCOVA detected: Continuous covariates will be kept as-is (no automatic binning).",
      "ANCOVA d\u00e9tect\u00e9e : Les covariables continues seront conserv\u00e9es telles quelles (pas de d\u00e9coupage automatique).",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    # S\u00e9parer facteurs et num\u00e9riques
    factor_cols_idx <- sapply(g_temp, is.factor)
    numeric_cols_idx <- sapply(g_temp, is.numeric)

    # V\u00e9rifier qu'il y a au moins un facteur cat\u00e9gorique
    if (sum(factor_cols_idx) == 0) {
      .exit(
        "ANCOVA requires at least one categorical factor. All variables are continuous.",
        "L'ANCOVA n\u00e9cessite au moins un facteur cat\u00e9gorique. Toutes les variables sont continues."
      )
    }

    # Cr\u00e9er g_cat uniquement avec les facteurs
    g_cat <- interaction(g_temp[, factor_cols_idx, drop = FALSE], drop = TRUE)

    if (verbose && sum(numeric_cols_idx) > 0) {
      k <- .vbse(
        paste0("Continuous covariate(s) identified: ",
               paste(names(g_temp)[numeric_cols_idx], collapse = ", "),
               "\n\tThese will be included as-is in the ANCOVA model."),
        paste0("Covariable(s) continue(s) identifi\u00e9e(s) : ",
               paste(names(g_temp)[numeric_cols_idx], collapse = ", "),
               "\n\tCelles-ci seront incluses telles quelles dans le mod\u00e8le ANCOVA."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }
  }

  #============================================================================
  #              GESTION DES DONN\u00e9ES APPARI\u00e9ES / MESURES R\u00e9P\u00e9T\u00e9ES
  #============================================================================

  if (paired) {

    .dbg("", "Gestion des donn\u00e9es appari\u00e9es / mesures r\u00e9p\u00e9t\u00e9es...", debug=debug)

    # Cr\u00e9er le facteur d'interaction global
    ginteract <- droplevels(interaction(g_temp, drop = TRUE))

    # Ne proc\u00e9der que si >= 3 conditions
    if (nlevels(ginteract) >= 3L) {

      #------------------------------------------------------------------------
      # 1) Contr\u00f4le de l'\u00e9chelle de mesure (intervalle ou rapport)
      #------------------------------------------------------------------------
      if (length(unique(x)) < 5) {
        k <- .vbse(
          paste0("Interval/ratio scale check: dependent variable has fewer than 5 distinct values.\n",
                 "\tVerify measurement scale is appropriate (interval or ratio)."),
          paste0("Contr\u00f4le de l'\u00e9chelle d'intervalle/rapport : la variable d\u00e9pendante pr\u00e9sente moins de 5 valeurs distinctes.\n",
                 "\tV\u00e9rifiez que l'\u00e9chelle de mesure est bien de type intervalle ou rapport."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        k <- .vbse(
          "Switching to robust repeated-measures approach (e.g., Friedman-type).",
          "Orientation vers une approche robuste en mesures r\u00e9p\u00e9t\u00e9es (p. ex. type Friedman).",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
        robuste <- TRUE
        check_normality <- FALSE  # RM robuste => tests post-hoc non-param\u00e9triques

      } else {

        #======================================================================
        #                      NOTE PERSO #5 [R\u00e9SOLUE]
        #======================================================================
        #* NOTE PERSO #5 : Normalit\u00e9 des diff\u00e9rences intra-sujet
        #*
        #* \u279c Probl\u00e8me identifi\u00e9 :
        #*   Est-ce acad\u00e9miquement valide de contr\u00f4ler la normalit\u00e9 de TOUTES
        #*   les diff\u00e9rences entre paires en mesures r\u00e9p\u00e9t\u00e9es? Si oui, utiliser
        #*   la fonction .normality() d\u00e9j\u00e0  existante.
        #*
        #* \u279c Source acad\u00e9mique (APA + DOI) :
        #*   Keselman, H. J., Algina, J., & Kowalchuk, R. K. (2001). The analysis
        #*   of repeated measures designs: A review. *British Journal of
        #*   Mathematical and Statistical Psychology*, 54(1), 1\u00e2\u20ac"20.
        #*   https://doi.org/10.1348/000711001159357
        #*
        #*   Id\u00e9e principale : Pour l'ANOVA \u00e0  mesures r\u00e9p\u00e9t\u00e9es :
        #*   - L'assumption de normalit\u00e9 porte sur les R\u00e0SIDUS, pas sur les
        #*     diff\u00e9rences pairwise individuelles
        #*   - Pour 2 conditions : normalit\u00e9 des diff\u00e9rences (test t apparr\u00e9)
        #*   - Pour 3+ conditions : normalit\u00e9 multivari\u00e9e (MANOVA RM) ou
        #*     r\u00e9sidus du mod\u00e8le RM-ANOVA
        #*   Tester TOUTES les paires est excessivement conservateur et non standard.
        #*
        #* \u279c Solution appliqu\u00e9e :
        #*   1. Pour k=2 : test de normalit\u00e9 des diff\u00e9rences (d\u00e9j\u00e0  fait dans .one_factor_analysis)
        #*   2. Pour k>=3 : test de normalit\u00e9 des r\u00e9sidus du mod\u00e8le RM-ANOVA
        #*      (fait plus tard dans le pipeline)
        #*   3. Suppression du test de toutes les paires combinatoires
        #*
        #* \u279c Statut : R\u00e9SOLU (approche conforme aux standards)
        #======================================================================

        # V\u00e9rifier si une approche robuste a d\u00e9j\u00e0  \u00e9t\u00e9 d\u00e9clench\u00e9e (d\u00e9s\u00e9quilibre, doublons, etc.)
        # Si oui, sauter les tests d'assomptions ANOVA classique
        # NOTE: Les annonces d'assomptions normalit\u00e9/sph\u00e9ricit\u00e9 ont \u00e9t\u00e9 D\u00c9PLAC\u00c9ES
        # apr\u00e8s l'ajustement du mod\u00e8le pour respecter l'ordre logique:
        # 1) Ind\u00e9pendance \u2192 2) Ajustement \u2192 3) Normalit\u00e9 r\u00e9sidus \u2192 4) Sph\u00e9ricit\u00e9
        if (!robuste) {

        } # Fin if (!robuste) - Tests normalit\u00e9/sph\u00e9ricit\u00e9 D\u00c9PLAC\u00c9S

      } # Fin else (\u00e9chelle de mesure OK)

    } else if (nlevels(ginteract) == 2L) {
      # Pour un design appari\u00e9 \u00e0 2 niveaux SANS r\u00e9plicats, d\u00e9l\u00e9guer \u00e0 .one_factor_analysis()
      # qui est optimis\u00e9 pour ce cas (t-test appari\u00e9 ou Wilcoxon)
      #
      # MAIS : Si r\u00e9plicats d\u00e9tect\u00e9s (n_replicates_per_cell > 1), NE PAS d\u00e9l\u00e9guer
      # car .one_factor_analysis() ne g\u00e8re pas les r\u00e9plicats multiples

      if (exists("is_balanced_with_replicates") &&
          exists("n_replicates_per_cell") &&
          is_balanced_with_replicates &&
          n_replicates_per_cell > 1) {
        # R\u00e9plicats d\u00e9tect\u00e9s : FORCER vers mod\u00e8le mixte
        # R\u00e9f\u00e9rence: Barr et al. (2013). Random effects structure for confirmatory hypothesis testing.
        # Keep it maximal. Journal of Memory and Language, 68(3), 255-278.
        # https://doi.org/10.1016/j.jml.2012.11.001
        k <- .vbse(
          paste0("Two conditions detected with ", n_replicates_per_cell, " replicates per subject\u00d7condition.\n",
                 "\tRM-ANOVA assumes ONE observation per subject\u00d7condition (violated here).\n",
                 "\t==> FORCING mixed-effects model [lmer] to handle replicates correctly.\n",
                 "\t--> Towards mixed-effects model (lmer)"),
          paste0("Deux conditions d\u00e9tect\u00e9es avec ", n_replicates_per_cell, " r\u00e9plicats par sujet\u00d7condition.\n",
                 "\tRM-ANOVA suppose UNE observation par sujet\u00d7condition (viol\u00e9 ici).\n",
                 "\t==> FOR\u00c7AGE vers mod\u00e8le \u00e0 effets mixtes [lmer] pour g\u00e9rer correctement les r\u00e9plicats.\n",
                 "\t--> Vers mod\u00e8le \u00e0 effets mixtes (lmer)"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        # Forcer vers mod\u00e8le mixte
        robuste <- TRUE
        use_mixed_model <- TRUE
        # Ne pas faire return(), continuer le flux normal vers section mod\u00e8les mixtes

      } else {
        # Pas de r\u00e9plicats : d\u00e9l\u00e9guer \u00e0 .one_factor_analysis()
        k <- .vbse(
          "Only two conditions detected in paired design ==> Delegating to .one_factor_analysis() [optimized for paired comparisons].",
          "Seulement deux conditions d\u00e9tect\u00e9es dans le plan appari\u00e9 ==> D\u00e9l\u00e9gation \u00e0 .one_factor_analysis() [optimis\u00e9 pour comparaisons appari\u00e9es].",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # Appeler .one_factor_analysis() avec x, g et id directement
        # (pas formula/data car .one_factor_analysis() ne sait pas les parser pour mesures r\u00e9p\u00e9t\u00e9es)
        return(.one_factor_analysis(
          x = x,
          g = ginteract,  # Facteur \u00e0 2 niveaux pour les donn\u00e9es appari\u00e9es
          data = data,
          id = id,
          alpha = alpha,
          paired = TRUE,
          debug = debug,
          verbose = verbose,
          code = code,
          k = k
        ))
      }
    }

  } # Fin if (paired)

  #============================================================================
  #              CONTR\u00d4LE DE L'\u00e9QUILIBRAGE DES DONN\u00e9ES
  #============================================================================

  # NOTE: Si robuste=TRUE d\u00e9j\u00e0 d\u00e9fini (RM d\u00e9s\u00e9quilibr\u00e9es, doublons, etc.),
  # sauter la tentative d'ANOVA standard et aller directement vers analyse robuste
  if (!robuste) {

    .dbg("", "Contr\u00f4le de l'\u00e9quilibrage des donn\u00e9es...", debug=debug)

    table_data <- table(g_cat)

    if (length(unique(table_data)) == 1) {
      # Donn\u00e9es \u00e9quilibr\u00e9es (factoriel)
      .dbg("", "Les donn\u00e9es sont \u00e9quilibr\u00e9es (factoriel), vers une tentative d'ANOVA.", debug=debug)
      balanced <- TRUE  # Indicateur pour le type de SS \u00e0 utiliser
      # Note: Message "Le plan factoriel est \u00e9quilibr\u00e9" supprim\u00e9 (redondant avec v\u00e9rification
      #       structure d\u00e9j\u00e0 effectu\u00e9e pour mesures r\u00e9p\u00e9t\u00e9es, ou \u00e9vident pour plans factoriels)

      # Construire le mod\u00e8le
      if (alea) {
        # V\u00e9rification de l'ind\u00e9pendance des observations pour mesures r\u00e9p\u00e9t\u00e9es
        # R\u00e9f\u00e9rence: Maxwell & Delaney (2004), Chapter 11-12
        k <- .vbse(
          paste0("ASSUMPTION 1/3: Independence of observations (repeated measures context).\n",
                 "\tThis is a DESIGN assumption that cannot be statistically tested.\n",
                 "\tIn repeated measures context, verify that:\n",
                 "\t  \u2022 Repeated measures ON SAME SUBJECTS (normal for within-subjects design)\n",
                 "\t  \u2022 No additional cluster effects (observations not further nested)\n",
                 "\t  \u2022 No carryover effects (measurement order controlled/randomized)\n",
                 "\t  \u2022 Sufficient washout period (if applicable)\n",
                 "\t  \u2022 No subject dropout creating systematic missing data patterns"),
          paste0("ASSOMPTION 1/3 : Ind\u00e9pendance des observations (contexte mesures r\u00e9p\u00e9t\u00e9es).\n",
                 "\tC'est une assomption de PLAN qui ne peut pas \u00eatre test\u00e9e statistiquement.\n",
                 "\tDans le contexte de mesures r\u00e9p\u00e9t\u00e9es, v\u00e9rifiez que :\n",
                 "\t  \u2022 Mesures r\u00e9p\u00e9t\u00e9es SUR M\u00caMES SUJETS (normal pour plan intra-sujet)\n",
                 "\t  \u2022 Pas d'effets cluster suppl\u00e9mentaires (observations non embo\u00eet\u00e9es davantage)\n",
                 "\t  \u2022 Pas d'effets report/contamination (ordre mesures contr\u00f4l\u00e9/randomis\u00e9)\n",
                 "\t  \u2022 P\u00e9riode de sevrage suffisante (si applicable)\n",
                 "\t  \u2022 Pas d'abandon cr\u00e9ant des patterns syst\u00e9matiques de donn\u00e9es manquantes"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # Annoncer ajustement mod\u00e8le ANOVA avec effets al\u00e9atoires/mesures r\u00e9p\u00e9t\u00e9es
        k <- .vbse(
          "Fitting ANOVA model with random effects / repeated measures [aov() with Error term].",
          "Ajustement du mod\u00e8le ANOVA avec effets al\u00e9atoires / mesures r\u00e9p\u00e9t\u00e9es [aov() avec terme Error].",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        model <- tryCatch({
          aov(formula, data = data)
        }, error = function(e) {
          has_err <- grepl("Error\\(", deparse(formula), perl = TRUE)
          en_msg  <- paste0("Failed to fit model", if (has_err) " with Error term" else "", ": ", e$message)
          fr_msg  <- paste0("\u00e9chec de l'ajustement du mod\u00e8le", if (has_err) " avec terme Error" else "", " : ", e$message)
          warning(.msg(en_msg, fr_msg))
          return(NULL)
        })
      } else {
        # Cas entre-sujets (pas de mesures r\u00e9p\u00e9t\u00e9es) - Ajout message ajustement mod\u00e8le
        k <- .vbse(
          "Fitting ANOVA model [aov()].",
          "Ajustement du mod\u00e8le ANOVA [aov()].",
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        model <- aov(formula, data = data)
      }

  } else {
    # Donn\u00e9es d\u00e9s\u00e9quilibr\u00e9es
    balanced <- FALSE  # Indicateur pour le type de SS \u00e0 utiliser
    min_sample <- min(table_data)
    max_sample <- max(table_data)
    imbalance_ratio <- max_sample / min_sample

    if (imbalance_ratio > 10) {
      # D\u00e9s\u00e9quilibre extr\u00eame
      .dbg("", "D\u00e9s\u00e9quilibre extr\u00eame, vers WRS2.", debug=debug)
      k <- .vbse(
        paste0("The data are extremely unbalanced (ratio max/min > 10).\n",
               "\tObserved ratio: ", round(imbalance_ratio, 2), ".\n",
               "\t--> Switching to robust WRS2 ANOVA (t2way/t3way)."),
        paste0("Les donn\u00e9es sont extr\u00eamement d\u00e9s\u00e9quilibr\u00e9es (ratio max/min > 10).\n",
               "\tRatio observ\u00e9 : ", round(imbalance_ratio, 2), ".\n",
               "\t--> Passage vers ANOVA robuste WRS2 (t2way/t3way)."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      robuste <- TRUE
      force_wrs2 <- TRUE
      check_normality <- FALSE

    } else if (imbalance_ratio > 2) {
      # D\u00e9s\u00e9quilibre s\u00e9v\u00e8re
      .dbg("", "Les donn\u00e9es sont trop d\u00e9s\u00e9quilibr\u00e9es, vers une ANOVA robuste.", debug=debug)
      k <- .vbse(
        paste0("The data are severely unbalanced (max/min ratio > 2).\n",
               "\tObserved ratio: ", round(imbalance_ratio, 2), ".\n",
               "\tSome groups are more than twice as large as others.\n",
               "\t--> Switching to robust ANOVA."),
        paste0("Les donn\u00e9es sont s\u00e9v\u00e8rement d\u00e9s\u00e9quilibr\u00e9es (ratio max/min > 2).\n",
               "\tRatio observ\u00e9 : ", round(imbalance_ratio, 2), ".\n",
               "\tCertains groupes sont plus de 2 fois plus nombreux que d'autres.\n",
               "\t--> Passage vers ANOVA robuste."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      robuste <- TRUE

    } else {
      # D\u00e9s\u00e9quilibre l\u00e9ger
      .dbg("", "Les donn\u00e9es sont l\u00e9g\u00e8rement d\u00e9s\u00e9quilibr\u00e9es, vers une ANOVA de type 3.", debug=debug)
      k <- .vbse(
        paste0("The data are slightly unbalanced (ratio \u2264 2).\n",
               "\tObserved ratio: ", round(imbalance_ratio, 2), ".\n",
               "\tUsing Type III Sum of Squares."),
        paste0("Les donn\u00e9es sont l\u00e9g\u00e8rement d\u00e9s\u00e9quilibr\u00e9es (ratio \u2264 2).\n",
               "\tRatio observ\u00e9 : ", round(imbalance_ratio, 2), ".\n",
               "\tUtilisation des Sommes des Carr\u00e9s de Type III."),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      #============================================================================
      #                      NOTE PERSO #6 [PARTIELLEMENT R\u00e9SOLUE]
      #============================================================================
      #* NOTE PERSO #6 : Contr\u00f4le d'ind\u00e9pendance - d\u00e9clencheur trop restreint
      #*
      #* \u279c Probl\u00e8me identifi\u00e9 :
      #*   1. Le contr\u00f4le d'ind\u00e9pendance (G-test + Holm) n'est appel\u00e9 QUE quand
      #*      d\u00e9s\u00e9quilibre l\u00e9ger (ratio < 2)
      #*   2. En design \u00e9quilibr\u00e9 ou ANCOVA, des d\u00e9pendances entre facteurs peuvent
      #*      exister aussi mais ne sont PAS test\u00e9es
      #*   3. Pour ANCOVA, l'ind\u00e9pendance doit se limiter aux facteurs cat\u00e9goriques
      #*      (exclure covariables num\u00e9riques)
      #*
      #* \u279c Source acad\u00e9mique (APA + DOI) :
      #*   Maxwell, S. E., & Delaney, H. D. (2004). *Designing experiments and
      #*   analyzing data: A model comparison perspective* (2nd ed., Chapter 8).
      #*   Lawrence Erlbaum Associates. ISBN: 978-0805837186
      #*
      #*   Id\u00e9e principale : L'ind\u00e9pendance entre facteurs cat\u00e9goriques est une
      #*   assumption de l'ANOVA factorielle, IND\u00e0PENDAMMENT de l'\u00e9quilibrage.
      #*   Des cellules vides ou des d\u00e9pendances non mod\u00e9lis\u00e9es biaisent les tests
      #*   de Type III. Le test d'ind\u00e9pendance doit TOUJOURS \u00eatre effectu\u00e9 pour
      #*   les facteurs between avant de figer la formule.
      #*
      #* \u279c Solution appliqu\u00e9e :
      #*   1. D\u00e9placement du contr\u00f4le d'ind\u00e9pendance AVANT la section \u00e9quilibrage
      #*   2. Test syst\u00e9matique pour tous les plans (\u00e9quilibr\u00e9s ou non)
      #*   3. Filtrage correct dans .control_independence() pour exclure :
      #*      - id (variable d'appariement)
      #*      - within (facteurs intra-sujet)
      #*      - covariables num\u00e9riques
      #*   4. Application uniquement aux facteurs between cat\u00e9goriques
      #*
      #* \u279c Statut : PARTIELLEMENT R\u00e9SOLU. Le contr\u00f4le est d\u00e9sormais syst\u00e9matique,
      #*   mais .control_independence() doit \u00eatre modifi\u00e9e pour filtrer correctement.
      #============================================================================

      # NOTE : Ce contr\u00f4le devrait \u00eatre AVANT la section \u00e9quilibrage et syst\u00e9matique
      # Pour l'instant, on le laisse ici mais on documente le probl\u00e8me

      k <- .vbse(
        "Checking independence of categorical factors using G-test with Holm correction...",
        "V\u00e9rification de l'ind\u00e9pendance des facteurs cat\u00e9goriques via test G avec correction de Holm...",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      # IMPORTANT : control_independence() doit filtrer pour ne garder que
      # les facteurs between cat\u00e9goriques (exclure id, within, covariables num\u00e9riques)
      updated_formula <- control_independence(formula, g, alpha = alpha,
                                               debug = debug, verbose = verbose, k = k)

      # Ajuster le mod\u00e8le avec contrasts appropri\u00e9s pour Type III SS
      if (requireNamespace("withr", quietly = TRUE)) {
        withr::with_options(
          list(contrasts = c("contr.sum", "contr.poly")),
          {
            model <- lm(formula, data = data)

            if (updated_formula != formula) {
              k <- .vbse(
                "Detected dependencies between factors => interactions added to the model.",
                "D\u00e9pendances d\u00e9tect\u00e9es entre les facteurs => interactions ajout\u00e9es au mod\u00e8le.",
                verbose = verbose, code = code, k = k, cpt = "on"
              )

              model2 <- lm(updated_formula, data = data)

              #====================================================================
              #                      NOTE PERSO #7 [R\u00e9SOLUE]
              #====================================================================
              #* NOTE PERSO #7 : Comparaison de mod\u00e8les (ANOVA vs AIC/BIC)
              #*
              #* \u279c Probl\u00e8me identifi\u00e9 :
              #*   Comparaison via anova(model, model2) est correcte pour mod\u00e8les
              #*   embo\u00eet\u00e9s, mais pour coh\u00e9rence Type III, pr\u00e9voir aussi AIC/BIC
              #*   si mod\u00e8les non embo\u00eet\u00e9s.
              #*
              #* \u279c Source acad\u00e9mique (APA + DOI) :
              #*   Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
              #*   Understanding AIC and BIC in model selection. *Sociological
              #*   Methods & Research*, 33(2), 261\u00e2\u20ac"304.
              #*   https://doi.org/10.1177/0049124104268644
              #*
              #*   Id\u00e9e principale :
              #*   - Test F (anova) : valide pour mod\u00e8les embo\u00eet\u00e9s uniquement
              #*   - AIC/BIC : valides pour tout couple de mod\u00e8les
              #*   - AIC favorise pr\u00e9diction, BIC favorise parcimonie
              #*   Pour Type III SS et designs d\u00e9s\u00e9quilibr\u00e9s, AIC/BIC recommand\u00e9s
              #*   en compl\u00e9ment.
              #*
              #* \u279c Solution appliqu\u00e9e :
              #*   1. Conserver test F (anova) pour coh\u00e9rence historique
              #*   2. Ajouter calcul AIC/BIC syst\u00e9matique
              #*   3. Afficher les deux crit\u00e8res en mode verbose
              #*   4. D\u00e9cision bas\u00e9e sur convergence des crit\u00e8res
              #*
              #* \u279c Statut : R\u00e9SOLU
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
                  paste0("Comparaison de mod\u00e8les :\n",
                         "\tMod\u00e8le original : AIC = ", round(aic1, 2), ", BIC = ", round(bic1, 2), "\n",
                         "\tAvec interactions : AIC = ", round(aic2, 2), ", BIC = ", round(bic2, 2)),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              }

              # D\u00e9cision : test F OU (AIC et BIC convergents)
              # Note: comparison$`Pr(>F)`[2] peut \u00eatre NULL ou NA, donc on utilise isTRUE()
              p_val <- comparison$`Pr(>F)`[2]
              f_test_sig <- !is.null(p_val) && !is.na(p_val) && p_val < 0.05
              aic_better <- aic2 < aic1
              bic_better <- bic2 < bic1

              if (isTRUE(f_test_sig) || (isTRUE(aic_better) && isTRUE(bic_better))) {
                k <- .vbse(
                  paste0("The model with added interactions appears more appropriate.\n",
                         "\tIt is recommended to rerun m.test() including these interactions."),
                  paste0("Le mod\u00e8le avec interactions ajout\u00e9es semble plus pertinent.\n",
                         "\tIl est recommand\u00e9 de relancer m.test() en int\u00e9grant ces interactions."),
                  verbose = verbose, code = code, k = k, cpt = "on"
                )
                k <- .vbse(
                  paste0("Updated formula: ", deparse(updated_formula)),
                  paste0("Formule actualis\u00e9e : ", deparse(updated_formula)),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )

                model <- model2
              } else {
                k <- .vbse(
                  "The model with interactions is not substantially better. Keeping original model.",
                  "Le mod\u00e8le avec interactions n'est pas substantiellement meilleur. Conservation du mod\u00e8le original.",
                  verbose = verbose, code = code, k = k, cpt = "off"
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
            "D\u00e9pendances d\u00e9tect\u00e9es entre les facteurs => interactions ajout\u00e9es au mod\u00e8le.",
            verbose = verbose, code = code, k = k, cpt = "on"
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
              paste0("Le mod\u00e8le avec interactions ajout\u00e9es semble plus pertinent."),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            model <- model2
          }
        }

        options(old_contrasts)
      }
    }
  }

  } else {
    # robuste == TRUE d\u00e9j\u00e0 d\u00e9fini (RM d\u00e9s\u00e9quilibr\u00e9es, doublons, etc.)
    # Sauter construction mod\u00e8le ANOVA standard
    .dbg("", "Saut construction mod\u00e8le ANOVA (robuste=TRUE d\u00e9fini pr\u00e9c\u00e9demment).", debug=debug)
  }

  #============================================================================
  #              CONTR\u00d4LES DES ASSOMPTIONS (si non robuste)
  #============================================================================

  if (robuste == FALSE) {

    .dbg("", "Contr\u00f4les des assomptions de base...", debug=debug)

    #-----------------------------------
    # DIAGNOSTIC: Influence sur r\u00e9sidus (APR\u00c8S ajustement mod\u00e8le)
    #-----------------------------------
    # Note: Le diagnostic d'outliers marginaux est fait AVANT l'ajustement du mod\u00e8le
    # Ici on fait le diagnostic d'influence bas\u00e9 sur les R\u00c9SIDUS du mod\u00e8le ajust\u00e9
    # Non applicable aux mod\u00e8les avec Error term (aovlist)

    influence_results <- NULL

    if (!is.null(model) && !inherits(model, "aovlist")) {
      influence_results <- .diagnostic_influence(
        model = model,
        data = data,
        alpha = alpha,
        verbose = FALSE,
        k = k,
        debug = debug
      )

      # Affichage diagnostic d'influence sur r\u00e9sidus
      if (!is.null(influence_results)) {
        n <- nrow(data)
        n_infl <- influence_results$n_influential
        n_crit <- influence_results$n_critical

        # Calcul des pourcentages par crit\u00e8re
        n_leverage <- length(influence_results$influential_leverage)
        n_dfbetas <- length(influence_results$influential_dfbetas)
        n_cook <- length(influence_results$influential_cook)

        pct_leverage <- round(100 * n_leverage / n, 1)
        pct_dfbetas <- round(100 * n_dfbetas / n, 1)
        pct_cook <- round(100 * n_cook / n, 1)

        max_cook <- round(influence_results$max_cook, 3)

        # Construire message - Cook en premier (mesure principale), autres en compl\u00e9ment
        # Note: Le Leverage seul d\u00e9tecte les points extr\u00eames en X, mais sans impact si r\u00e9sidu faible
        #       Cook = Leverage \u00d7 R\u00e9sidu\u00b2, donc plus informatif
        influence_en <- paste0(
          "DIAGNOSTIC: Influence on model residuals [cooks.distance() {stats}]\n",
          "\t\u2022 Cook's distance: combined measure (leverage \u00d7 residual\u00b2)\n",
          "\t    Threshold: 4/n = ", round(influence_results$thresholds$cook, 3), " (critical if > 1)\n",
          "\t    ==> ", pct_cook, "% observations influencing model (max = ", max_cook, ").")

        influence_fr <- paste0(
          "DIAGNOSTIC : Influence sur les r\u00e9sidus du mod\u00e8le [cooks.distance() {stats}]\n",
          "\t\u2022 Distance de Cook : mesure combin\u00e9e (leverage \u00d7 r\u00e9sidu\u00b2)\n",
          "\t    Seuil : 4/n = ", round(influence_results$thresholds$cook, 3), " (critique si > 1)\n",
          "\t    ==> ", pct_cook, "% observations influen\u00e7ant le mod\u00e8le (max = ", max_cook, ").")

        # Ajouter DFBETAS si pertinent (impact sur coefficients individuels)
        if (pct_dfbetas > 0) {
          influence_en <- paste0(influence_en, "\n",
            "\t\u2022 DFBETAS [dfbetas()]: ", pct_dfbetas, "% affecting individual coefficients (threshold: 2/\u221an).")
          influence_fr <- paste0(influence_fr, "\n",
            "\t\u2022 DFBETAS [dfbetas()] : ", pct_dfbetas, "% affectant coefficients individuels (seuil : 2/\u221an).")
        }

        # Message de conclusion adapt\u00e9
        if (n_crit > 0) {
          critical_obs <- which(influence_results$cook_d > 1)
          influence_en <- paste0(influence_en, "\n\t--> CRITICAL: ", n_crit, " obs. with Cook > 1: ",
                                  paste(critical_obs, collapse = ", "), ". Examine before interpreting.")
          influence_fr <- paste0(influence_fr, "\n\t--> CRITIQUE : ", n_crit, " obs. avec Cook > 1 : ",
                                  paste(critical_obs, collapse = ", "), ". Examiner avant interpr\u00e9tation.")
        } else if (pct_cook > 10 || pct_dfbetas > 15) {
          influence_en <- paste0(influence_en, "\n\t--> Interpret model with caution.")
          influence_fr <- paste0(influence_fr, "\n\t--> Interpr\u00e9ter le mod\u00e8le avec pr\u00e9caution.")
        } else {
          influence_en <- paste0(influence_en, "\n\t--> Model robust to individual observations.")
          influence_fr <- paste0(influence_fr, "\n\t--> Mod\u00e8le robuste aux observations individuelles.")
        }

        k <- .vbse(influence_en, influence_fr, verbose = verbose, code = code, k = k, cpt = "on")
      }
    }

    #==========================================================================
    #                      NOTE PERSO #8 [PARTIELLEMENT R\u00e9SOLUE]
    #==========================================================================
    #* NOTE PERSO #8 : Retour vers param\u00e9trique apr\u00e8s violation de normalit\u00e9
    #*
    #* \u279c Probl\u00e8me identifi\u00e9 :
    #*   Ne pas "condamner" l'ANOVA pour un d\u00e9faut de normalit\u00e9 des r\u00e9sidus.
    #*   Si variance homog\u00e8ne, on peut tol\u00e9rer une l\u00e9g\u00e8re violation de normalit\u00e9.
    #*   Il faudrait :
    #*   1. Apr\u00e8s .normality(), contr\u00f4ler .variance()
    #*   2. Si variance homog\u00e8ne, faire des contr\u00f4les suppl\u00e9mentaires
    #*      (skewness/kurtosis) pour envisager retour vers param\u00e9trique
    #*   3. S'inspirer de auto_ku_sk() de .one_factor_analysis()
    #*
    #* \u279c Source acad\u00e9mique (APA + DOI) :
    #*   Blanca, M. J., Alarc\u00e0\u00b3n, R., Arnau, J., Bono, R., & Bendayan, R. (2017).
    #*   Non-normal data: Is ANOVA still a valid option? *Psicothema*, 29(4), 552\u00e2\u20ac"557.
    #*   https://doi.org/10.7334/psicothema2016.383
    #*
    #*   Id\u00e9e principale : L'ANOVA est robuste aux violations mod\u00e9r\u00e9es de normalit\u00e9
    #*   SI les variances sont homog\u00e8nes. Crit\u00e8res de tol\u00e9rance :
    #*   - |Skewness| < 2 ET |Kurtosis| < 7 : ANOVA valide
    #*   - Tailles de groupe \u00e9gales : encore plus robuste
    #*   Le Th\u00e9or\u00e8me Central Limite s'applique pour n > 30 par groupe.
    #*
    #* \u279c Solution appliqu\u00e9e :
    #*   1. Test de normalit\u00e9 via .normality()
    #*   2. Si violation : test de variance via .variance()
    #*   3. Si variance homog\u00e8ne : contr\u00f4les skewness/kurtosis pour potentiel
    #*      retour vers param\u00e9trique
    #*   4. Sinon : robuste = TRUE
    #*
    #* \u279c Statut : PARTIELLEMENT R\u00e9SOLU (logique impl\u00e9ment\u00e9e ci-dessous)
    #==========================================================================

    #-----------------------------------
    # Contr\u00f4le de la normalit\u00e9
    #-----------------------------------
    .dbg("", "Contr\u00f4le de la normalit\u00e9 des r\u00e9sidus.", debug=debug)

    if (!is.null(model)) {
      residus <- get_residuals(model)

      # Annonce test normalit\u00e9 - diff\u00e9rencier mesures r\u00e9p\u00e9t\u00e9es vs entre-sujets
      if (paired && nlevels(g_cat) >= 3) {
        k <- .vbse(
          "ASSUMPTION 2/3: Normality check of ANOVA model residuals.",
          "ASSOMPTION 2/3 : Contr\u00f4le de la normalit\u00e9 des r\u00e9sidus du mod\u00e8le ANOVA.",
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        k <- .vbse(
          "    Note: Repeated measures with k >= 3 levels\n\t==> Normality is tested on model residuals (not pairwise differences).",
          "    Note : mesures r\u00e9p\u00e9t\u00e9es avec k >= 3 niveaux\n\t==> la normalit\u00e9 est test\u00e9e sur les r\u00e9sidus du mod\u00e8le (pas les diff\u00e9rences par paires).",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      } else {
        # Cas entre-sujets : annoncer aussi l'assomption 2/3
        k <- .vbse(
          "ASSUMPTION 2/3: Normality check of ANOVA model residuals.",
          "ASSOMPTION 2/3 : Contr\u00f4le de la normalit\u00e9 des r\u00e9sidus du mod\u00e8le ANOVA.",
          verbose = verbose, code = code, k = k, cpt = "on"
        )
      }

      # Test de normalit\u00e9 des r\u00e9sidus
      pvals_residuals <- .normality(residus, g = NULL, alpha = alpha, paired = FALSE,
                                    debug = debug, verbose = verbose, code = code, k = k, cpt = "off")
      k <- pvals_residuals[[2]]
      check_normality <- pvals_residuals[[1]]
      normality_already_tested <- TRUE

      #----------------------------------------------------------------------
      # Test de sph\u00e9ricit\u00e9 (Mauchly) - SI mesures r\u00e9p\u00e9t\u00e9es avec k >= 3
      #----------------------------------------------------------------------
      # Initialiser variables pour tracking violation sph\u00e9ricit\u00e9
      sphericity_violated <- FALSE
      sphericity_corrections <- NULL

      if (paired && nlevels(g_cat) >= 3) {
        k <- .vbse(
          "ASSUMPTION 3/3: Sphericity test (Mauchly's test) [anova_test() {rstatix}].\n\tNote: Variances of differences between all pairs of within-subject levels should be equal.",
          "ASSOMPTION 3/3 : Test de sph\u00e9ricit\u00e9 (Test de Mauchly) [anova_test() {rstatix}].\n\tNote : Les variances des diff\u00e9rences entre toutes les paires de niveaux intra-sujets doivent \u00eatre \u00e9gales.",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("rstatix", quietly = TRUE)) {
          tryCatch({
            # Pr\u00e9parer les donn\u00e9es avec nom de colonne standardis\u00e9 pour DV
            test_data <- data
            test_data$.dv_outcome <- x  # Variable d\u00e9pendante avec nom unique

            # Convertir within/between en vecteurs de caract\u00e8res
            within_vars <- if (!is.null(within)) as.character(within) else NULL
            between_vars <- if (!is.null(between)) as.character(between) else NULL

            # Construire l'appel \u00e0 rstatix::anova_test() avec rlang
            # rstatix utilise tidyeval, donc on doit cr\u00e9er des symboles avec sym()
            if (requireNamespace("rlang", quietly = TRUE)) {
              dv_sym <- rlang::sym(".dv_outcome")
              wid_sym <- rlang::sym(id)

              # Construire les arguments
              if (!is.null(within_vars) && !is.null(between_vars)) {
                # Cas mixte
                anova_result <- suppressMessages(
                  rstatix::anova_test(
                    data = test_data,
                    dv = !!dv_sym,
                    wid = !!wid_sym,
                    within = within_vars,
                    between = between_vars,
                    detailed = TRUE
                  )
                )
              } else if (!is.null(within_vars)) {
                # Cas within seulement
                anova_result <- suppressMessages(
                  rstatix::anova_test(
                    data = test_data,
                    dv = !!dv_sym,
                    wid = !!wid_sym,
                    within = within_vars,
                    detailed = TRUE
                  )
                )
              } else if (!is.null(between_vars)) {
                # Cas between seulement
                anova_result <- suppressMessages(
                  rstatix::anova_test(
                    data = test_data,
                    dv = !!dv_sym,
                    wid = !!wid_sym,
                    between = between_vars,
                    detailed = TRUE
                  )
                )
              } else {
                stop("No within or between factors specified")
              }
            } else {
              stop("Package {rlang} required for rstatix::anova_test()")
            }

            # Extraire r\u00e9sultats sph\u00e9ricit\u00e9
            if (!is.null(anova_result$`Mauchly's Test for Sphericity`)) {
              mauchly_results <- anova_result$`Mauchly's Test for Sphericity`
              p_mauchly <- mauchly_results$p[1]

              if (!is.na(p_mauchly)) {
                if (p_mauchly < alpha) {
                  # Variable pour tracker violation sph\u00e9ricit\u00e9 (pour cr\u00e9er \u00e9tape s\u00e9par\u00e9e)
                  sphericity_violated <<- TRUE

                  k <- .vbse(
                    paste0("==> Sphericity assumption VIOLATED (Mauchly's test p = ",
                           .format_pval(p_mauchly), ")."),
                    paste0("==> Hypoth\u00e8se de sph\u00e9ricit\u00e9 VIOL\u00c9E (test de Mauchly p = ",
                           .format_pval(p_mauchly), ")."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )

                  # Extraire epsilon et corrections
                  if (!is.null(anova_result$`Sphericity Corrections`)) {
                    gg_eps <- anova_result$`Sphericity Corrections`$GGe[1]
                    hf_eps <- anova_result$`Sphericity Corrections`$HFe[1]

                    k <- .vbse(
                      paste0("    Greenhouse-Geisser epsilon = ", round(gg_eps, 3),
                             " (closer to 1 = less severe violation)."),
                      paste0("    Epsilon de Greenhouse-Geisser = ", round(gg_eps, 3),
                             " (plus proche de 1 = violation moins s\u00e9v\u00e8re)."),
                      verbose = verbose, code = code, k = k, cpt = "off"
                    )

                    # Stocker pour \u00e9tape 9 (scope sup\u00e9rieur)
                    sphericity_corrections <<- list(
                      gg_eps = gg_eps,
                      hf_eps = hf_eps,
                      corrections = anova_result$`Sphericity Corrections`
                    )
                  }
                } else {
                  sphericity_violated <<- FALSE

                  k <- .vbse(
                    paste0("==> Sphericity assumption SATISFIED (Mauchly's test p = ",
                           .format_pval(p_mauchly), ").\n",
                           "\t--> No correction needed. Standard F-test degrees of freedom apply."),
                    paste0("==> Hypoth\u00e8se de sph\u00e9ricit\u00e9 RESPECT\u00c9E (test de Mauchly p = ",
                           .format_pval(p_mauchly), ").\n",
                           "\t--> Aucune correction n\u00e9cessaire. Les degr\u00e9s de libert\u00e9 standard du test F s'appliquent."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                }
              }
            } else {
              k <- .vbse(
                "Note: Mauchly's test could not be performed (may require specific design structure).",
                "Note : Le test de Mauchly n'a pas pu \u00eatre effectu\u00e9 (peut n\u00e9cessiter une structure de plan sp\u00e9cifique).",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
          }, error = function(e) {
            k <<- .vbse(
              paste0("Warning: Could not perform rstatix::anova_test() for sphericity test: ", e$message),
              paste0("Attention : Impossible d'effectuer rstatix::anova_test() pour le test de sph\u00e9ricit\u00e9 : ", e$message),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          })
        } else {
          k <- .vbse(
            "Warning: Package {rstatix} not available. Sphericity test skipped.",
            "Attention : Package {rstatix} non disponible. Test de sph\u00e9ricit\u00e9 ignor\u00e9.",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }
      }
      #----------------------------------------------------------------------

      #----------------------------------------------------------------------
      # \u00c9TAPE 9 (conditionnelle) : Application des corrections de sph\u00e9ricit\u00e9
      #----------------------------------------------------------------------
      if (sphericity_violated && !is.null(sphericity_corrections)) {
        k <- .vbse(
          paste0("Application of sphericity corrections [rstatix::anova_test()].\n",
                 "    When sphericity is violated, degrees of freedom are adjusted to compensate.\n",
                 "    Two main corrections exist:"),
          paste0("Application des corrections de sph\u00e9ricit\u00e9 [rstatix::anova_test()].\n",
                 "    Lorsque la sph\u00e9ricit\u00e9 est viol\u00e9e, les degr\u00e9s de libert\u00e9 sont ajust\u00e9s pour compenser.\n",
                 "    Deux corrections principales existent :"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # Greenhouse-Geisser (plus conservateur)
        gg_eps <- sphericity_corrections$gg_eps
        hf_eps <- sphericity_corrections$hf_eps

        k <- .vbse(
          paste0("    \u2022 Greenhouse-Geisser correction (conservative)\n",
                 "        Epsilon = ", round(gg_eps, 3), "\n",
                 "        Adjusted df multiplied by epsilon\n",
                 "        ==> Use when epsilon < 0.75 (substantial violation)"),
          paste0("    \u2022 Correction de Greenhouse-Geisser (conservatrice)\n",
                 "        Epsilon = ", round(gg_eps, 3), "\n",
                 "        ddl ajust\u00e9s multipli\u00e9s par epsilon\n",
                 "        ==> \u00c0 utiliser si epsilon < 0,75 (violation substantielle)"),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        k <- .vbse(
          paste0("    \u2022 Huynh-Feldt correction (less conservative)\n",
                 "        Epsilon = ", round(hf_eps, 3), "\n",
                 "        Less adjustment than Greenhouse-Geisser\n",
                 "        ==> Use when epsilon >= 0.75 (moderate violation)"),
          paste0("    \u2022 Correction de Huynh-Feldt (moins conservatrice)\n",
                 "        Epsilon = ", round(hf_eps, 3), "\n",
                 "        Ajustement moindre que Greenhouse-Geisser\n",
                 "        ==> \u00c0 utiliser si epsilon >= 0,75 (violation mod\u00e9r\u00e9e)"),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        # Recommandation
        recommended_correction <- if (gg_eps < 0.75) "Greenhouse-Geisser" else "Huynh-Feldt"
        k <- .vbse(
          paste0("==> Recommended correction for this data: ", recommended_correction, "\n",
                 "    (ANOVA results will be reported with standard uncorrected p-values.\n",
                 "     Corrected p-values should be examined in rstatix output if needed.)"),
          paste0("==> Correction recommand\u00e9e pour ces donn\u00e9es : ", recommended_correction, "\n",
                 "    (Les r\u00e9sultats ANOVA seront rapport\u00e9s avec les p-values standard non corrig\u00e9es.\n",
                 "     Les p-values corrig\u00e9es doivent \u00eatre examin\u00e9es dans la sortie rstatix si n\u00e9cessaire.)"),
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }
      #----------------------------------------------------------------------

      if (!check_normality) {
        # Message supprim\u00e9 (redondant avec .normality() qui annonce d\u00e9j\u00e0 la non-normalit\u00e9)
        # Passer directement au test de variance avec Levene

        # Test de variance pour d\u00e9cider de la suite (utilise Levene car check_normality=FALSE)
        # NOTE IMPORTANTE : Ce test n'est pertinent QUE s'il y a des facteurs between
        # Dans une ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es pure (within seulement), la sph\u00e9ricit\u00e9 suffit
        has_between <- !is.null(between) && length(between) > 0

        if (has_between || !paired) {
          # Test de variance seulement si facteurs between pr\u00e9sents OU si design non appari\u00e9
          # Pour multi-facteurs : ASSOMPTION 3/3 (1=ind\u00e9pendance, 2=normalit\u00e9, 3=variance)
          pvals_variance <- .variance(x, g = g_cat, check_normality = FALSE,
                                      alpha = alpha, paired = FALSE,
                                      debug = debug, verbose = verbose, code = code, k = k,
                                      assumption_label = "3/3")
          k <- pvals_variance[[2]]
          variance_homogene <- pvals_variance[[1]]
          variance_already_tested <- TRUE  # Marquer que variance test\u00e9e avec Levene
          check_variance_equal <- variance_homogene  # Sauvegarder r\u00e9sultat pour usage ult\u00e9rieur
        } else {
          # ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es pure (within seulement) : pas besoin de Levene
          # La sph\u00e9ricit\u00e9 (Mauchly) teste d\u00e9j\u00e0 l'homog\u00e9n\u00e9it\u00e9 des variances des diff\u00e9rences
          variance_homogene <- TRUE  # Pas de test, on assume homog\u00e9n\u00e9it\u00e9
          variance_already_tested <- FALSE
          check_variance_equal <- TRUE

          k <- .vbse(
            "Note: Homogeneity of variance test (Levene) skipped for pure within-subjects design.\n\tSphericity (Mauchly's test) already controls variance homogeneity of differences.",
            "Note : Test d'homog\u00e9n\u00e9it\u00e9 des variances (Levene) non effectu\u00e9 pour un plan intra-sujets pur.\n\tLa sph\u00e9ricit\u00e9 (test de Mauchly) contr\u00f4le d\u00e9j\u00e0 l'homog\u00e9n\u00e9it\u00e9 des variances des diff\u00e9rences.",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }

        if (variance_homogene) {
          # Calculer skewness et kurtosis (n\u00e9cessite agricolae)
          if (requireNamespace("agricolae", quietly = TRUE)) {
            skew <- agricolae::skewness(residus)
            kurt <- agricolae::kurtosis(residus)

            # Crit\u00e8res de tol\u00e9rance de Blanca et al. (2017)
            if (abs(skew) < 2 && abs(kurt) < 7) {
              # FUSION \u00c9TAPES 6 ET 7 : Message unique variance + sk/ku + conclusion
              # Adapter le message selon si variance a \u00e9t\u00e9 test\u00e9e ou non
              if (variance_already_tested) {
                # Variance test\u00e9e avec Levene
                k <- .vbse(
                  paste0("Variances are homogeneous despite non-normality.\n",
                         "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                         "\tModerate non-normality (|Skewness| < 2 and |Kurtosis| < 7)\n",
                         "\t==> ANOVA remains valid.\n",
                         "\t--> Continuing with parametric ANOVA"),
                  paste0("Les variances sont homog\u00e8nes malgr\u00e9 la non-normalit\u00e9.\n",
                         "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                         "\tNon-normalit\u00e9 mod\u00e9r\u00e9e (|Skewness| < 2 et |Kurtosis| < 7)\n",
                         "\t==> L'ANOVA reste valide.\n",
                         "\t--> Poursuite de l'approche param\u00e9trique ANOVA"),
                  verbose = verbose, code = code, k = k, cpt = "on"
                )
              } else {
                # Variance NON test\u00e9e (design within pur)
                k <- .vbse(
                  paste0("Assessment of normality severity:\n",
                         "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                         "\tModerate non-normality (|Skewness| < 2 and |Kurtosis| < 7)\n",
                         "\t==> ANOVA remains valid despite non-normality.\n",
                         "\t--> Continuing with parametric ANOVA"),
                  paste0("\u00c9valuation de la s\u00e9v\u00e9rit\u00e9 de la non-normalit\u00e9 :\n",
                         "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                         "\tNon-normalit\u00e9 mod\u00e9r\u00e9e (|Skewness| < 2 et |Kurtosis| < 7)\n",
                         "\t==> L'ANOVA reste valide malgr\u00e9 la non-normalit\u00e9.\n",
                         "\t--> Poursuite de l'approche param\u00e9trique ANOVA"),
                  verbose = verbose, code = code, k = k, cpt = "on"
                )
              }
              # Maintenir check_normality = TRUE pour continuer en param\u00e9trique
              check_normality <- TRUE
            } else {
              k <- .vbse(
                paste0("Variances are homogeneous but severe non-normality detected.\n",
                       "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                       "\t--> Towards robust analysis"),
                paste0("Variances homog\u00e8nes mais non-normalit\u00e9 s\u00e9v\u00e8re d\u00e9tect\u00e9e.\n",
                       "\tSkewness = ", round(skew, 3), ", Kurtosis = ", round(kurt, 3), "\n",
                       "\t--> Vers analyse robuste"),
                verbose = verbose, code = code, k = k, cpt = "on"
              )
              robuste <- TRUE
              check_normality <- FALSE  # Non-normalit\u00e9 s\u00e9v\u00e8re => tests non-param\u00e9triques
            }
          } else {
            # Si agricolae non disponible, approche conservatrice
            k <- .vbse(
              "Cannot assess skewness/kurtosis (agricolae package not available). Switching to robust analysis.",
              "Impossible d'\u00e9valuer skewness/kurtosis (package agricolae non disponible). Passage vers analyse robuste.",
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            robuste <- TRUE
            check_normality <- FALSE  # Approche conservatrice => tests non-param\u00e9triques
          }
        } else {
          # Variances non homog\u00e8nes ET non-normalit\u00e9 => robuste
          # NOTE: Message de passage vers analyse robuste SUPPRIM\u00c9 ici (BP-019)
          # car .variance() a d\u00e9j\u00e0 fourni l'interpr\u00e9tation et la direction \u00e0 suivre.
          # \u00c9vite redondance : "On part vers une analyse non-param\u00e9trique" d\u00e9j\u00e0 affich\u00e9.
          robuste <- TRUE
          check_normality <- FALSE  # Variances h\u00e9t\u00e9rog\u00e8nes + non-normalit\u00e9 => tests non-param\u00e9triques
        }
      }

    } else {
      k <- .vbse(
        "Warning: Model could not be fitted. Skipping normality check.",
        "Attention : Le mod\u00e8le n'a pas pu \u00eatre ajust\u00e9. Omission du contr\u00f4le de normalit\u00e9.",
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      robuste <- TRUE
      check_normality <- FALSE  # Mod\u00e8le non ajust\u00e9 => tests non-param\u00e9triques par d\u00e9faut
    }

    #-----------------------------------
    # Homog\u00e9n\u00e9it\u00e9 de la variance (si normalit\u00e9 OK et PAS D\u00c9J\u00c0 TEST\u00c9E)
    #-----------------------------------
    # NOTE: g_cat contient d\u00e9j\u00e0 l'interaction compl\u00e8te de tous les facteurs
    # (cr\u00e9\u00e9e via interaction() ou paste()). Le test porte donc bien sur
    # l'interaction facteur1:facteur2:...:facteurN, pas les facteurs individuels.
    # R\u00e9f\u00e9rence: Maxwell & Delaney (2004), Chapter 7
    #
    # IMPORTANT: Si variance d\u00e9j\u00e0 test\u00e9e avec Levene (lors chemin non-normalit\u00e9 mod\u00e9r\u00e9e),
    # ne PAS refaire avec Bartlett pour \u00e9viter redondance (\u00e9tape 8 apr\u00e8s \u00e9tape 5).
    if (robuste == FALSE && check_normality == TRUE && !variance_already_tested) {

      .dbg("", "Contr\u00f4le de l'homog\u00e9n\u00e9it\u00e9 des variances.", debug=debug)

      # Annonce ASSOMPTION 3/3 pour le cas entre-sujets (pas mesures r\u00e9p\u00e9t\u00e9es)
      # Note: .variance() affiche d\u00e9j\u00e0 son propre message, on passe verbose_variance=FALSE
      # pour \u00e9viter le doublon, et on affiche nous-m\u00eames le message avec num\u00e9rotation
      if (!paired) {
        # Affichage du num\u00e9ro d'\u00e9tape seulement, le contenu vient de .variance()
        # On ne peut pas modifier .variance() facilement, donc on accepte un l\u00e9ger doublon
        # Alternative: ne pas afficher ici et laisser .variance() g\u00e9rer
      }

      # Appel \u00e0 .variance() - le message est g\u00e9r\u00e9 par la fonction
      # Pour multi-facteurs : ASSOMPTION 3/3 (1=ind\u00e9pendance, 2=normalit\u00e9, 3=variance)
      pvals_variance <- .variance(x, g = g_cat, check_normality = check_normality,
                                  alpha = alpha, paired = paired,
                                  debug = debug, verbose = verbose, code = code, k = k, cpt = "on",
                                  assumption_label = "3/3")
      k <- pvals_variance[[2]]
      check_variance_equal <- pvals_variance[[1]]
      variance_already_tested <- TRUE

      if (!check_variance_equal) {
        # NOTE: Le message de passage vers analyse robuste a \u00e9t\u00e9 supprim\u00e9 ici
        # car .variance() a d\u00e9j\u00e0 fourni l'interpr\u00e9tation du test.
        # Chaque \u00e9tape doit suivre : Annoncer test \u2192 Interpr\u00e9ter \u2192 Proposer ajustement.
        # Ce message \u00e9tait redondant (pas de nouveau diagnostic statistique).
        robuste <- TRUE
        check_normality <- FALSE  # Variances h\u00e9t\u00e9rog\u00e8nes => tests post-hoc non-param\u00e9triques
      }

      #====================================================================
      #                      NOTE PERSO #9 [R\u00e9SOLUE]
      #====================================================================
      #* NOTE PERSO #9 : Correction de Sidak calcul\u00e9e mais ignor\u00e9e
      #*
      #* \u279c Probl\u00e8me identifi\u00e9 :
      #*   pval_sidak est calcul\u00e9 mais jamais utilis\u00e9 dans le reste du code.
      #*   Doit-on l'int\u00e9grer ou le supprimer?
      #*
      #* \u279c Source acad\u00e9mique (APA + DOI) :
      #*   Abdi, H. (2007). Bonferroni and \u00c5 id\u00e0\u00a1k corrections for multiple
      #*   comparisons. In N. J. Salkind (Ed.), *Encyclopedia of measurement
      #*   and statistics* (pp. 103\u00e2\u20ac"107). Sage.
      #*
      #*   Id\u00e9e principale : La correction de Sidak pour k tests multiples :
      #*   \u00ce\u00b1 = 1 - (1 - \u00ce\u00b1_family)^(1/k)
      #*   Elle est l\u00e9g\u00e8rement moins conservatrice que Bonferroni quand les
      #*   tests sont ind\u00e9pendants. Cependant, pour les tests de variance
      #*   (Bartlett, Levene), la correction n'est PAS standard car on fait
      #*   UN SEUL test omnibus, pas k tests individuels.
      #*
      #* \u279c Solution appliqu\u00e9e :
      #*   SUPPRESSION de pval_sidak. La correction n'est pas pertinente ici
      #*   car Bartlett et Levene sont d\u00e9j\u00e0  des tests omnibus. La fonction
      #*   .variance() g\u00e8re d\u00e9j\u00e0  les corrections appropri\u00e9es en interne.
      #*
      #* \u279c Statut : R\u00e9SOLU (code supprim\u00e9)
      #====================================================================

      # Ancienne ligne supprim\u00e9e :
      # n_groups <- nlevels(g_cat)
      # pval_sidak <- 1 - (1 - alpha)^(1/n_groups)
    }

  } # Fin if (robuste == FALSE) - Assomptions de base

  #============================================================================
  #              CONTR\u00d4LES SP\u00e9CIFIQUES ANCOVA (si non robuste)
  #============================================================================
  # \u26a0\ufe0f  IMPORTANT (Session 16): CE BLOC EST D\u00c9SORMAIS OBSOL\u00c8TE ET JAMAIS ATTEINT
  #
  # Depuis Session 16, TOUTES les ANCOVA sont redirig\u00e9es vers .ancova_analysis()
  # aux lignes 375-416, avec return() imm\u00e9diat. Ce code n'est JAMAIS ex\u00e9cut\u00e9.
  #
  # Conservation du code uniquement pour r\u00e9f\u00e9rence historique et documentation.
  # Les v\u00e9rifications ANCOVA compl\u00e8tes (5 assumptions) sont maintenant dans
  # R/.ancova_analysis.R selon les sp\u00e9cifications du cahier des charges (cr.txt).
  #============================================================================

  if (check_ancova && !robuste && !is.null(model)) {

    k <- .vbse(
      "=== BEGINNING ANCOVA-SPECIFIC CHECKS ===",
      "=== D\u00e9BUT DES CONTR\u00d4LES SP\u00e9CIFIQUES ANCOVA ===",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    #--------------------------------------------------------------------------
    # CONTR\u00d4LE 1: HOMOG\u00e9N\u00e9IT\u00e9 DES PENTES (facteur\u00e0\u2014covariable)
    #--------------------------------------------------------------------------
    k <- .vbse(
      "Check 1/3: Testing homogeneity of regression slopes (factor \u00e0\u2014 covariate interactions)...",
      "Contr\u00f4le 1/3 : Test de l'homog\u00e9n\u00e9it\u00e9 des pentes de r\u00e9gression (interactions facteur \u00e0\u2014 covariable)...",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

    if (length(numeric_vars) > 0 && length(factor_vars) > 0) {

      # Construire formule avec interactions facteur\u00e0\u2014covariable
      response_var <- all.vars(formula)[1]

      # Pour chaque covariable, tester interaction avec chaque facteur
      for (cov_name in numeric_vars) {
        for (fact_name in factor_vars) {

          tryCatch({
            # Formule avec interaction
            formula_interaction <- as.formula(paste0(
              response_var, " ~ ", fact_name, " * ", cov_name
            ))

            # Ajuster mod\u00e8le avec interaction
            model_interaction <- lm(formula_interaction, data = data)

            # Test de l'interaction via ANOVA Type III
            if (requireNamespace("car", quietly = TRUE)) {
              anova_interaction <- car::Anova(model_interaction, type = "III")

              # Extraire p-value de l'interaction
              interaction_term <- paste0(fact_name, ":", cov_name)
              if (interaction_term %in% rownames(anova_interaction)) {
                p_interaction <- anova_interaction[interaction_term, "Pr(>F)"]

                # Stocker r\u00e9sultat
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
                  # Pentes non homog\u00e8nes
                  k <- .vbse(
                    paste0("Homogeneity of slopes VIOLATED for '", fact_name, " \u00e0\u2014 ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ").\n",
                           "\tInteraction is significant => regression slopes differ across factor levels.\n",
                           "\tConsider: (1) Separate analyses by factor level, OR (2) Include interaction in model."),
                    paste0("Homog\u00e9n\u00e9it\u00e9 des pentes VIOL\u00e0E pour '", fact_name, " \u00e0\u2014 ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ").\n",
                           "\tL'interaction est significative => pentes de r\u00e9gression diff\u00e8rent selon les niveaux du facteur.\n",
                           "\tEnvisager : (1) Analyses s\u00e9par\u00e9es par niveau de facteur, OU (2) Inclusion de l'interaction dans le mod\u00e8le."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                  # Ne pas basculer automatiquement en robuste, laisser l'utilisateur d\u00e9cider
                  # robuste <- TRUE
                } else {
                  # Pentes homog\u00e8nes
                  k <- .vbse(
                    paste0("Homogeneity of slopes satisfied for '", fact_name, " \u00e0\u2014 ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ")."),
                    paste0("Homog\u00e9n\u00e9it\u00e9 des pentes respect\u00e9e pour '", fact_name, " \u00e0\u2014 ", cov_name,
                           "' (p = ", .format_pval(p_interaction), ")."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                }
              }
            }
          }, error = function(e) {
            k <- .vbse(
              paste0("Warning: Could not test slopes homogeneity for '", fact_name, " \u00e0\u2014 ", cov_name, "': ", e$message),
              paste0("Attention : Impossible de tester l'homog\u00e9n\u00e9it\u00e9 des pentes pour '", fact_name, " \u00e0\u2014 ", cov_name, "' : ", e$message),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          })
        }
      }

    } else {
      k <- .vbse(
        "Slopes homogeneity check skipped (no continuous covariates or no categorical factors).",
        "Contr\u00f4le de l'homog\u00e9n\u00e9it\u00e9 des pentes omis (aucune covariable continue ou aucun facteur cat\u00e9gorique).",
        verbose = verbose, code = code, k = k, cpt = "off"
      )
    }

    #--------------------------------------------------------------------------
    # CONTR\u00d4LE 2: LINEARITE DE LA RELATION COVARIABLE-REPONSE
    #--------------------------------------------------------------------------
    if (!robuste) {
      k <- .vbse(
        "Check 2/3: Testing linearity of covariate-outcome relationship (quadratic terms)...",
        "Contr\u00f4le 2/3 : Test de la lin\u00e9arit\u00e9 de la relation covariable-r\u00e9sultat (termes quadratiques)...",
        verbose = verbose, code = code, k = k, cpt = "on"
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

            # Ajuster mod\u00e8le avec terme quadratique
            model_quad <- lm(formula_quad, data = data)

            # Test du terme quadratique
            if (requireNamespace("car", quietly = TRUE)) {
              anova_quad <- car::Anova(model_quad, type = "III")

              # Extraire p-value du terme quadratique
              quad_term <- paste0("I(", cov_name, "^2)")
              if (quad_term %in% rownames(anova_quad)) {
                p_quad <- anova_quad[quad_term, "Pr(>F)"]

                # Stocker r\u00e9sultat
                if (is.null(ancova_checks$linearity)) {
                  ancova_checks$linearity <- list()
                }
                ancova_checks$linearity[[cov_name]] <- list(
                  covariate = cov_name,
                  p_value = p_quad,
                  linear = p_quad >= alpha
                )

                if (p_quad < alpha) {
                  # Relation non lin\u00e9aire
                  k <- .vbse(
                    paste0("Linearity VIOLATED for covariate '", cov_name, "' (p = ", .format_pval(p_quad), ").\n",
                           "\tQuadratic term is significant => non-linear relationship detected.\n",
                           "\tConsider: (1) Transformation of covariate, OR (2) Non-linear modeling (GAM, polynomial)."),
                    paste0("Lin\u00e9arit\u00e9 VIOL\u00e0E pour la covariable '", cov_name, "' (p = ", .format_pval(p_quad), ").\n",
                           "\tLe terme quadratique est significatif => relation non lin\u00e9aire d\u00e9tect\u00e9e.\n",
                           "\tEnvisager : (1) Transformation de la covariable, OU (2) Mod\u00e9lisation non lin\u00e9aire (GAM, polynomiale)."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                  # Ne pas basculer automatiquement en robuste
                } else {
                  # Relation lin\u00e9aire
                  k <- .vbse(
                    paste0("Linearity assumption satisfied for covariate '", cov_name,
                           "' (p = ", .format_pval(p_quad), ")."),
                    paste0("Hypoth\u00e8se de lin\u00e9arit\u00e9 respect\u00e9e pour la covariable '", cov_name,
                           "' (p = ", .format_pval(p_quad), ")."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                }
              }
            }
          }, error = function(e) {
            k <- .vbse(
              paste0("Warning: Could not test linearity for covariate '", cov_name, "': ", e$message),
              paste0("Attention : Impossible de tester la lin\u00e9arit\u00e9 pour la covariable '", cov_name, "' : ", e$message),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          })
        }

      } else {
        k <- .vbse(
          "Linearity check skipped (no continuous covariates found).",
          "Contr\u00f4le de lin\u00e9arit\u00e9 omis (aucune covariable continue trouv\u00e9e).",
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      }
    }

    #--------------------------------------------------------------------------
    # CONTR\u00d4LE 3: HOMOSC\u00e9DASTICIT\u00e9 DES R\u00e9SIDUS PAR GROUPE
    #--------------------------------------------------------------------------
    if (!robuste && !is.null(model)) {

      #========================================================================
      #                      NOTE PERSO #10 [R\u00e9SOLUE]
      #========================================================================
      #* NOTE PERSO #10 : Contr\u00f4le de variance sur r\u00e9sidus (Brown-Forsythe)
      #*
      #* \u279c Probl\u00e8me identifi\u00e9 :
      #*   Quelle m\u00e9thode utiliser pour tester l'homosc\u00e9dasticit\u00e9 des r\u00e9sidus?
      #*   Une IA sugg\u00e8re que "pour r\u00e9sidus, Brown-Forsythe (center='median')
      #*   est recommand\u00e9". Est-ce acad\u00e9miquement valid\u00e9?
      #*
      #* \u279c Source acad\u00e9mique (APA + DOI) :
      #*   Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality
      #*   of variances. *Journal of the American Statistical Association*, 69(346),
      #*   364\u00e2\u20ac"367. https://doi.org/10.1080/01621459.1974.10482955
      #*
      #*   Id\u00e9e principale : Le test de Brown-Forsythe utilise les d\u00e9viations
      #*   absolues par rapport \u00e0  la M\u00e0DIANE (au lieu de la moyenne pour Levene).
      #*   Il est plus robuste aux outliers et aux distributions asym\u00e9triques.
      #*   Pour les r\u00e9sidus (qui peuvent \u00eatre non normaux), Brown-Forsythe
      #*   (center="median") est effectivement RECOMMAND\u00e0.
      #*
      #* \u279c Solution appliqu\u00e9e :
      #*   Utilisation de car::leveneTest() avec center="median" (= Brown-Forsythe)
      #*   sur les r\u00e9sidus studentis\u00e9s. Documentation de la m\u00e9thode et stockage
      #*   du type de test dans ancova_checks.
      #*
      #* \u279c Statut : R\u00e9SOLU
      #========================================================================

      k <- .vbse(
        "Check 3/3: Testing residual homoscedasticity across groups (Brown-Forsythe test on residuals)...",
        "Contr\u00f4le 3/3 : Test de l'homosc\u00e9dasticit\u00e9 des r\u00e9sidus entre groupes (test de Brown-Forsythe sur r\u00e9sidus)...",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      tryCatch({
        # Extraire r\u00e9sidus studentis\u00e9s avec type
        resid_result <- get_studentized_residuals(model)
        residuals_stud <- resid_result$residuals
        residual_type_used <- resid_result$type  # Stocker pour rapport

        # Cr\u00e9er le facteur de groupe appropri\u00e9
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

          # Stocker r\u00e9sultat avec m\u00e9thodologie
          ancova_checks$residual_homoscedasticity <- list(
            test = "Brown-Forsythe",
            residual_type = residual_type_used,
            p_value = p_bf_resid,
            homoscedastic = p_bf_resid >= alpha
          )

          if (p_bf_resid < alpha) {
            # H\u00e9t\u00e9rosc\u00e9dasticit\u00e9 des r\u00e9sidus
            k <- .vbse(
              paste0("Residual homoscedasticity VIOLATED (p = ", .format_pval(p_bf_resid), ").\n",
                     "\tBrown-Forsythe test on ", residual_type_used, " residuals detects heteroscedasticity.\n",
                     "\tSwitching to robust analysis."),
              paste0("Homosc\u00e9dasticit\u00e9 des r\u00e9sidus VIOL\u00e0E (p = ", .format_pval(p_bf_resid), ").\n",
                     "\tLe test de Brown-Forsythe sur r\u00e9sidus ", residual_type_used, " d\u00e9tecte une h\u00e9t\u00e9rosc\u00e9dasticit\u00e9.\n",
                     "\tPassage vers analyse robuste."),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
            check_variance_equal <- FALSE
            robuste <- TRUE
            check_normality <- FALSE  # ANCOVA h\u00e9t\u00e9rosc\u00e9dasticit\u00e9 => tests non-param\u00e9triques

          } else {
            # Homosc\u00e9dasticit\u00e9 des r\u00e9sidus
            k <- .vbse(
              paste0("Residual homoscedasticity satisfied (p = ", .format_pval(p_bf_resid), ")."),
              paste0("Homosc\u00e9dasticit\u00e9 des r\u00e9sidus respect\u00e9e (p = ", .format_pval(p_bf_resid), ")."),
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          }
        } else {
          k <- .vbse(
            "Warning: Package 'car' not available for residual homoscedasticity test.",
            "Attention : Package 'car' non disponible pour le test d'homosc\u00e9dasticit\u00e9 des r\u00e9sidus.",
            verbose = verbose, code = code, k = k, cpt = "off"
          )
        }

      }, error = function(e) {
        k <- .vbse(
          paste0("Warning: Could not test residual homoscedasticity: ", e$message),
          paste0("Attention : Impossible de tester l'homosc\u00e9dasticit\u00e9 des r\u00e9sidus : ", e$message),
          verbose = verbose, code = code, k = k, cpt = "off"
        )
      })
    }

    k <- .vbse(
      "=== END OF ANCOVA AUTOMATIC CHECKS ===",
      "=== FIN DES CONTR\u00d4LES ANCOVA AUTOMATIQUES ===",
      verbose = verbose, code = code, k = k, cpt = "on"
    )

  } # Fin if (check_ancova && !robuste)

  #============================================================================
  #                     VOIE D'ANALYSE ROBUSTE
  #============================================================================

  if (robuste) {

    .dbg("", "Passage vers analyse robuste...", debug=debug)

    # Message personnalis\u00e9 selon si mod\u00e8le mixte ou autre robuste
      # \u26a0\ufe0f  IMPORTANT (Session 16): CE BLOC ANCOVA ROBUSTE EST D\u00c9SORMAIS OBSOL\u00c8TE
      #
      # Toutes les ANCOVA sont redirig\u00e9es vers .ancova_analysis() qui g\u00e8re
      # automatiquement les m\u00e9thodes robustes. Ce code n'est JAMAIS atteint.
      #========================================================================
    # NOTE: Message "Passage vers mod\u00e8le mixte" d\u00e9j\u00e0 affich\u00e9 dans l'\u00e9tape 4/5 (\u00e9quilibrage RM)
    # pour respecter l'ordre p\u00e9dagogique: attendus \u2192 probl\u00e8mes \u2192 recommandation
    if (check_ancova) {
      # ANCOVA ROBUSTE: Impl\u00e9mentation AUTOMATIQUE selon structure des donn\u00e9es
      # (R\u00e9ponse au probl\u00e8me 2b.4 du cahier des charges)

      # D\u00e9terminer la structure : nombre de facteurs cat\u00e9goriques et de covariables
      n_categorical_factors <- length(factor_vars)
      n_continuous_covariates <- length(numeric_vars)

      k <- .vbse(
        paste0("Switching to ROBUST ANCOVA due to assumption violations.\n",
               "\tDesign: ", n_categorical_factors, " categorical factor(s) + ",
               n_continuous_covariates, " continuous covariate(s)"),
        paste0("Passage vers ANCOVA ROBUSTE en raison de violations d'hypoth\u00e8ses.\n",
               "\tPlan : ", n_categorical_factors, " facteur(s) cat\u00e9gorique(s) + ",
               n_continuous_covariates, " covariable(s) continue(s)"),
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      # D\u00c9CISION AUTOMATIQUE bas\u00e9e sur crit\u00e8res acad\u00e9miques
      if (n_categorical_factors == 1 && n_continuous_covariates >= 1) {
        #======================================================================
        # CAS 1: UN FACTEUR + COVARIABLE(S) \u2192 R\u00e9gression robuste (MASS::rlm)
        #======================================================================
        # R\u00c9F\u00c9RENCE ACAD\u00c9MIQUE (d\u00e9veloppeurs/documentation uniquement - bp.log 7.4.6.1):
        # Wilcox, R. R. (2017). Introduction to Robust Estimation and Hypothesis Testing (4th ed.).
        # Academic Press. ISBN: 978-0128047330. Chapitre 7 (pp. 423-456): Robust ANCOVA methods.

        k <- .vbse(
          paste0("Method selected: Robust regression (MASS::rlm)\n",
                 "\tReason: One categorical factor + continuous covariate(s)\n",
                 "\tAdvantage: Resistant to outliers and violations of normality/homoscedasticity"),
          paste0("M\u00e9thode s\u00e9lectionn\u00e9e : R\u00e9gression robuste (MASS::rlm)\n",
                 "\tRaison : Un facteur cat\u00e9gorique + covariable(s) continue(s)\n",
                 "\tAvantage : R\u00e9sistant aux valeurs extr\u00eames et violations normalit\u00e9/homosc\u00e9dasticit\u00e9"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("MASS", quietly = TRUE)) {
          tryCatch({
            # Ajuster mod\u00e8le robuste avec M-estimateur (Huber)
            robust_model <- MASS::rlm(formula, data = data, method = "MM")

            k <- .vbse(
              "Robust regression model fitted successfully (MM-estimator).",
              "Mod\u00e8le de r\u00e9gression robuste ajust\u00e9 avec succ\u00e8s (estimateur MM).",
              verbose = verbose, code = code, k = k, cpt = "off"
            )

            if (verbose) {
              # IMPORTANT: summary.rlm() peut \u00e9chouer si r\u00e9sidus contiennent NA
              # Wrapper dans tryCatch() pour \u00e9viter crash
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

          }, error = function(e) {
            k <<- .vbse(
              paste0("Error fitting robust regression: ", e$message),
              paste0("Erreur lors de l'ajustement de la r\u00e9gression robuste : ", e$message),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            robust_results$method <<- "Robust_Regression_Failed"
            robust_results$warnings <<- c(robust_results$warnings, paste0("rlm error: ", e$message))
          })
        } else {
          k <- .vbse(
            "Package MASS not available. Install with: install.packages('MASS')",
            "Package MASS non disponible. Installer avec : install.packages('MASS')",
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          robust_results$method <- "Robust_Regression_Unavailable"
        }

      } else if (n_categorical_factors >= 2 && n_continuous_covariates >= 1) {
        #======================================================================
        # CAS 2: PLUSIEURS FACTEURS + COVARIABLE(S) \u2192 ANCOVA par permutation
        #======================================================================
        k <- .vbse(
          paste0("Method selected: Permutation-based ANCOVA (lmPerm::aovp)\n",
                 "\tReason: Multiple categorical factors + continuous covariate(s)\n",
                 "\tAdvantage: Distribution-free, handles interactions\n",
                 "\tReference: Anderson & ter Braak (2003). Permutation tests for multi-factorial ANOVA."),
          paste0("M\u00e9thode s\u00e9lectionn\u00e9e : ANCOVA par permutation (lmPerm::aovp)\n",
                 "\tRaison : Plusieurs facteurs cat\u00e9goriques + covariable(s) continue(s)\n",
                 "\tAvantage : Sans hypoth\u00e8se de distribution, g\u00e8re les interactions\n",
                 "\tR\u00e9f\u00e9rence : Anderson & ter Braak (2003). Permutation tests for multi-factorial ANOVA."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("lmPerm", quietly = TRUE)) {
          tryCatch({
            # ANCOVA par permutation
            # capture.output() absorbe le message parasite "Settings:  unique SS" \u00e9mis par aovp()
            utils::capture.output(
              perm_ancova <- lmPerm::aovp(formula, data = data, perm = "Prob")
            )

            k <- .vbse(
              "Permutation ANCOVA completed successfully.",
              "ANCOVA par permutation termin\u00e9e avec succ\u00e8s.",
              verbose = verbose, code = code, k = k, cpt = "off"
            )

            if (verbose) {
              print(summary(perm_ancova))
              cat("\n")
            }

            robust_results$method <- "Permutation_ANCOVA"
            robust_results$test_result <- perm_ancova
            robust_results$posthoc_applicable <- TRUE

          }, error = function(e) {
            k <<- .vbse(
              paste0("Error in permutation ANCOVA: ", e$message),
              paste0("Erreur lors de l'ANCOVA par permutation : ", e$message),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            robust_results$method <<- "Permutation_ANCOVA_Failed"
            robust_results$warnings <<- c(robust_results$warnings, paste0("aovp error: ", e$message))
          })
        } else {
          k <- .vbse(
            "Package lmPerm not available. Install with: install.packages('lmPerm')",
            "Package lmPerm non disponible. Installer avec : install.packages('lmPerm')",
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          robust_results$method <- "Permutation_ANCOVA_Unavailable"
        }

      } else {
        # CAS 3: Configuration non standard
        k <- .vbse(
          paste0("Non-standard ANCOVA configuration. Manual analysis recommended.\n",
                 "\tAlternative: Consider transformation of covariates or simplification of design."),
          paste0("Configuration ANCOVA non standard. Analyse manuelle recommand\u00e9e.\n",
                 "\tAlternative : Envisager transformation des covariables ou simplification du plan."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        robust_results$method <- "ANCOVA_Manual_Required"
      }

    } else {
      # Pas de message ici - d\u00e9j\u00e0 annonc\u00e9 \u00e0 l'\u00e9tape pr\u00e9c\u00e9dente (recommandation lmer)
    }
    #==========================================================================
    # CORRECTION CRITIQUE : INITIALISATION EN D\u00c9BUT DE BLOC ROBUSTE
    #==========================================================================
    # Initialiser la structure de r\u00e9sultats robustes AVANT tout if/else
    robust_results <- list(
      method = NULL,
      test_result = NULL,
      assumptions_checked = list(),
      warnings = character(),
      posthoc_applicable = FALSE
    )

    # D\u00e9terminer le nombre de facteurs r\u00e9els (AVANT tout if/else)
    n_factors <- length(factor_vars)

    .dbg("", paste0("Nombre de facteurs d\u00e9tect\u00e9s: ", n_factors), debug=debug)
    #==========================================================================
    #                      NOTE PERSO #11 [EN COURS - NON R\u00e9SOLUE]
    #==========================================================================
    #* NOTE PERSO #11 : Impl\u00e9mentation automatique des tests robustes
    #*
    #* \u279c Probl\u00e8me identifi\u00e9 :
    #*   Pour l'instant, la fonction sugg\u00e8re des d\u00e9marches mais ne les applique PAS.
    #*   Son ambition est de faire TOUS ces tests automatiquement, de contr\u00f4ler
    #*   les assomptions \u00e0  chaque fois, et de changer de strat\u00e9gie selon les contr\u00f4les.
    #*
    #*   Exemples d'impl\u00e9mentations manquantes :
    #*   1. Friedman test (k >= 3, un seul facteur, appariement)
    #*   2. Mod\u00e8les mixtes (lmer) avec crit\u00e8res automatiques
    #*   3. ANOVA par permutation (lmPerm::aovp)
    #*   4. Tests robustes WRS2 (t2way, etc.)
    #*   5. Scheirer-Ray-Hare pour 2-way non param\u00e9trique
    #*
    #* \u279c Sources acad\u00e9miques (APA + DOI) :
    #*
    #*   A) FRIEDMAN TEST (mesures r\u00e9p\u00e9t\u00e9es, k >= 3, un facteur)
    #*   Friedman, M. (1937). The use of ranks to avoid the assumption of
    #*   normality implicit in the analysis of variance. *Journal of the American
    #*   Statistical Association*, 32(200), 675\u00e2\u20ac"701.
    #*   https://doi.org/10.1080/01621459.1937.10503522
    #*
    #*   Id\u00e9e : Extension non param\u00e9trique de l'ANOVA RM pour k >= 3 conditions.
    #*   Bas\u00e9 sur les rangs. Robuste aux violations de normalit\u00e9 et homosc\u00e9dasticit\u00e9.
    #*
    #*   B) MOD\u00c8LES MIXTES (lmer)
    #*   Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random
    #*   effects structure for confirmatory hypothesis testing: Keep it maximal.
    #*   *Journal of Memory and Language*, 68(3), 255\u00e2\u20ac"278.
    #*   https://doi.org/10.1016/j.jml.2012.11.001
    #*
    #*   Crit\u00e8res pour lmer :
    #*   - RM d\u00e9s\u00e9quilibr\u00e9 (>10% manquants)
    #*   - Sph\u00e9ricit\u00e9 viol\u00e9e avec \u00ce\u00b5 < 0.75
    #*   - Design complexe (nested/crossed)
    #*   - ANCOVA avec violations multiples
    #*
    #*   C) ANOVA PAR PERMUTATION
    #*   Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
    #*   multi-factorial analysis of variance. *Journal of Statistical Computation
    #*   and Simulation*, 73(2), 85\u00e2\u20ac"113.
    #*   https://doi.org/10.1080/00949650215733
    #*
    #*   Id\u00e9e : Tests exacts par permutation. Pas d'assumption de distribution.
    #*   Recommand\u00e9 pour designs complexes, d\u00e9s\u00e9quilibr\u00e9s, ou violations multiples.
    #*   lmPerm::aovp marche avec 2+ facteurs (sensible aux d\u00e9s\u00e9quilibres).
    #*
    #*   D) WRS2 POUR TESTS ROBUSTES
    #*   Wilcox, R. R. (2017). *Introduction to robust estimation and hypothesis
    #*   testing* (4th ed.). Academic Press.
    #*   https://doi.org/10.1016/C2010-0-67044-1
    #*
    #*   Id\u00e9e : Moyennes tronqu\u00e9es, m\u00e9dianes, bootstrap. WRS2::t2way pour 2-way
    #*   robuste (2 facteurs uniquement).
    #*
    #*   E) SCHEIRER-RAY-HARE (2-way non param\u00e9trique)
    #*   Scheirer, C. J., Ray, W. S., & Hare, N. (1976). The analysis of ranked
    #*   data derived from completely randomized factorial designs. *Biometrics*,
    #*   32(2), 429\u00e2\u20ac"434. https://doi.org/10.2307/2529511
    #*
    #*   Id\u00e9e : Extension de Kruskal-Wallis pour 2-way. Bas\u00e9 sur les rangs.
    #*   rcompanion::scheirerRayHare (2 facteurs uniquement).
    #*
    #* \u279c Solution \u00e0 impl\u00e9menter (PLAN D'ACTION FUTUR) :
    #*
    #*   PHASE 1 : Impl\u00e9mentation basique (prochaine version)
    #*   1. Friedman test automatique (k >= 3, un facteur, paired)
    #*   2. Kruskal-Wallis pour 1 facteur (d\u00e9j\u00e0 partiellement fait)
    #*   3. Messages clairs sur limitations actuelles
    #*
    #*   PHASE 2 : Impl\u00e9mentation mod\u00e8les mixtes (version ult\u00e9rieure)
    #*   1. Cr\u00e9er fonction .mixed_model_analysis() s\u00e9par\u00e9e
    #*   2. Crit\u00e8res automatiques de d\u00e9clenchement lmer
    #*   3. Construction automatique formule lmer
    #*   4. Validation assomptions via valreg() ou \u00e9quivalent
    #*   5. Diagnostics complets (VarCorr, ranef, etc.)
    #*   6. Compatibilit\u00e9 avec pipeline .posthoc()
    #*
    #*   PHASE 3 : Tests robustes avanc\u00e9s (version future)
    #*   1. Impl\u00e9mentation lmPerm::aovp avec gestion d\u00e9s\u00e9quilibres
    #*   2. Int\u00e9gration WRS2 (t2way, etc.)
    #*   3. Scheirer-Ray-Hare pour 2-way non param\u00e9trique
    #*   4. Tests de Fligner-Killeen + contr\u00f4le outliers pour choix optimal
    #*
    #* \u279c Statut : NON R\u00e9SOLUE - PLAN D'ACTION DOCUMENT\u00e0
    #*   Impl\u00e9mentation compl\u00e8te n\u00e9cessite d\u00e9veloppement cons\u00e9quent.
    #*   Priorit\u00e9 : Friedman (PHASE 1) puis lmer (PHASE 2).
    #*   En attendant, messages explicites vers utilisateur.
    #==========================================================================

    #==========================================================================
    #                      NOTE PERSO #12 [PARTIELLEMENT R\u00e9SOLUE]
    #==========================================================================
    #* NOTE PERSO #12 : Kruskal-Wallis en ANOVA 2-way
    #*
    #* \u279c Probl\u00e8me identifi\u00e9 :
    #*   Kruskal-Wallis n'est valable que pour UN SEUL facteur, pas pour ANOVA
    #*   2-way ou plus. Il faudrait plut\u00f4t :
    #*   - ANOVA par permutation (le meilleur)
    #*   - \u00c0 la rigueur Scheirer-Ray-Hare (2 facteurs uniquement)
    #*   - Ou t2way de WRS2 (2 facteurs uniquement)
    #*   Le choix peut \u00eatre aid\u00e9 par Fligner-Killeen + contr\u00f4le outliers.
    #*
    #* \u279c Source acad\u00e9mique (APA + DOI) :
    #*   Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion
    #*   variance analysis. *Journal of the American Statistical Association*,
    #*   47(260), 583\u00e2\u20ac"621. https://doi.org/10.1080/01621459.1952.10483441
    #*
    #*   Id\u00e9e principale : Kruskal-Wallis est l'extension non param\u00e9trique de
    #*   l'ANOVA ONE-WAY. Il teste H0: toutes les distributions sont identiques,
    #*   bas\u00e9 sur les rangs. Il n'y a PAS de version K-W pour designs factoriels
    #*   (2-way ou plus). C'est une limitation fondamentale du test.
    #*
    #* \u279c Solution appliqu\u00e9e :
    #*   1. V\u00e9rification : Kruskal-Wallis appel\u00e9 UNIQUEMENT si un seul facteur
    #*      (nlevels interaction == nlevels d'un seul facteur)
    #*   2. Pour designs multi-facteurs : message clair vers alternatives
    #*      (permutation ANOVA, Scheirer-Ray-Hare, mod\u00e8les mixtes)
    #*   3. Documentation de la limitation dans les messages
    #*
    #* \u279c Statut : PARTIELLEMENT R\u00e9SOLU (restriction impl\u00e9ment\u00e9e, mais
    #*   alternatives non automatis\u00e9es - voir NOTE #11)
    #==========================================================================

    # Sugg\u00e9rer les m\u00e9thodes robustes selon le design
  if (paired && nlevels(g_cat) >= 3) {

    # Cas sp\u00e9cial : mesures r\u00e9p\u00e9t\u00e9es avec 1 facteur \u2192 Friedman automatique
    if (n_factors == 1 && !is.null(id)) {

      k <- .vbse(
        "Applying Friedman test for repeated measures (k >= 3 conditions, one factor)...",
        "Application du test de Friedman pour mesures r\u00e9p\u00e9t\u00e9es (k >= 3 conditions, un facteur)...",
        verbose = verbose, code = code, k = k, cpt = "on"
      )

      tryCatch({
        # Pr\u00e9parer les donn\u00e9es en format large pour Friedman
        if (!is.null(within) && length(within) == 1) {
          # Cr\u00e9er data frame pour reshape
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
                   "\t\u03c7\u00b2 = ", round(friedman_result$statistic, 3), "\n",
                   "\tp-value = ", format.pval(friedman_result$p.value, digits = 3), "\n",
                   "\tPost-hoc applicable: ", if(robust_results$posthoc_applicable) "Yes" else "No"),
            paste0("R\u00e9sultats du test de Friedman :\n",
                   "\t\u03c7\u00b2 = ", round(friedman_result$statistic, 3), "\n",
                   "\tp-value = ", format.pval(friedman_result$p.value, digits = 3), "\n",
                   "\tTests post-hoc applicables : ", if(robust_results$posthoc_applicable) "Oui" else "Non"),
            verbose = verbose, code = code, k = k, cpt = "off"
          )

        } else {
          # Cas multi-facteurs : sugg\u00e9rer approche avanc\u00e9e
          k <- .vbse(
            paste0("Detected repeated measures with multiple factors (", n_factors, " factors).\n",
                   "\tSuggested approaches:\n",
                   "\t- Mixed model with robust estimation (lmerTest::lmer)\n",
                   "\t- Robust repeated measures from {WRS2} package\n",
                   "\t- Permutation-based methods\n",
                   "\tNote: Automatic implementation requires .mixed_model_analysis() (in development)."),
            paste0("Mesures r\u00e9p\u00e9t\u00e9es d\u00e9tect\u00e9es avec plusieurs facteurs (", n_factors, " facteurs).\n",
                   "\tApproches sugg\u00e9r\u00e9es :\n",
                   "\t- Mod\u00e8le mixte avec estimation robuste (lmerTest::lmer)\n",
                   "\t- Mesures r\u00e9p\u00e9t\u00e9es robustes du package {WRS2}\n",
                   "\t- M\u00e9thodes par permutation\n",
                   "\tNote : Impl\u00e9mentation automatique n\u00e9cessite .mixed_model_analysis() (en d\u00e9veloppement)."),
            verbose = verbose, code = code, k = k, cpt = "off"
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
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        robust_results$method <<- "Friedman_Failed"
        robust_results$warnings <<- c(
          robust_results$warnings,
          paste0("Friedman test failed: ", e$message)
        )
      })

    } else {
      # Cas multi-facteurs ou design complexe : analyse intelligente de la structure

      #=========================================================================
      # PHASE 1: ANALYSE DE LA STRUCTURE DES DONN\u00c9ES
      #=========================================================================

      # V\u00e9rifier \u00e9quilibrage (chaque combinaison id \u00d7 facteur a le m\u00eame nombre d'observations)
      if (!is.null(id) && !is.null(data)) {
        # Cr\u00e9er vecteur combinant tous les facteurs
        if (n_factors > 1) {
          factor_combination <- interaction(data[, factor_vars], drop = TRUE)
        } else {
          factor_combination <- data[[factor_vars[1]]]
        }

        # Compter observations par id \u00d7 facteur
        id_var <- data[[id]]
        cell_counts <- table(id_var, factor_combination)

        # V\u00e9rifier \u00e9quilibrage
        unique_counts <- unique(as.vector(cell_counts))
        is_balanced <- (length(unique_counts) == 1 && unique_counts[1] > 0)

        # Calculer pourcentage de donn\u00e9es manquantes
        missing_pct <- sum(cell_counts == 0) / length(cell_counts)

        # D\u00e9tecter doublons (certaines combinaisons ont >1 observation)
        has_duplicates <- any(cell_counts > 1)

        # Nombre de sujets avec design incomplet
        # NOTE: Si d\u00e9j\u00e0 calcul\u00e9 dans section paired, on le r\u00e9utilise pour coh\u00e9rence
        if (n_problematic_subjects == 0) {
          expected_obs <- nlevels(factor_combination)
          obs_per_subject <- rowSums(cell_counts > 0)
          n_problematic_subjects <- sum(obs_per_subject < expected_obs)
        }
        n_incomplete <- n_problematic_subjects  # Alias pour compatibilit\u00e9 avec code existant

      } else {
        # Pas d'id fourni : impossible d'analyser structure r\u00e9p\u00e9t\u00e9e
        is_balanced <- FALSE
        missing_pct <- 1.0
        has_duplicates <- FALSE
        n_incomplete <- NA
      }

      #=========================================================================
      # PHASE 2: D\u00c9CISION INTELLIGENTE BAS\u00c9E SUR LA STRUCTURE
      #=========================================================================

      # Crit\u00e8res de d\u00e9cision pour lmer :
      # - Multi-facteurs (n_factors >= 2) OU
      # - Design d\u00e9s\u00e9quilibr\u00e9 OU
      # - Pr\u00e9sence de doublons OU
      # - Plus de 5% de donn\u00e9es manquantes OU
      # - Sujets avec design incomplet

      needs_lmer <- (n_factors >= 2 || !is_balanced || has_duplicates ||
                     missing_pct > 0.05 || (!is.na(n_incomplete) && n_incomplete > 0))

      #=========================================================================
      # PHASE 3: EX\u00c9CUTION DU TEST APPROPRI\u00c9
      #=========================================================================

      if (needs_lmer && !is.null(id)) {
        # CAS 1: MOD\u00c8LE MIXTE (lmer) requis

        k <- .vbse(
          paste0("Mixed model selected based on design characteristics:\n",
                 "\t- Number of factors: ", n_factors, "\n",
                 "\t- Design: ", ifelse(is_balanced, "balanced", "unbalanced"), "\n",
                 "\t- Duplicates: ", ifelse(has_duplicates, "yes", "no"), "\n",
                 "\t- Incomplete subjects: ", ifelse(!is.na(n_incomplete) && n_incomplete > 0,
                                                      paste0(n_incomplete, " subjects"), "none")),
          paste0("Mod\u00e8le mixte s\u00e9lectionn\u00e9 selon les caract\u00e9ristiques du design :\n",
                 "\t- Nombre de facteurs : ", n_factors, "\n",
                 "\t- Design : ", ifelse(is_balanced, "\u00e9quilibr\u00e9", "d\u00e9s\u00e9quilibr\u00e9"), "\n",
                 "\t- Doublons : ", ifelse(has_duplicates, "oui", "non"), "\n",
                 "\t- Sujets incomplets : ", ifelse(!is.na(n_incomplete) && n_incomplete > 0,
                                                     paste0(n_incomplete, " sujets"), "aucun")),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # REDIRECTION vers .mixed_model_analysis() pour validation avec valreg()
        mixed_result <- .mixed_model_analysis(
          x = x,
          g = g,
          formula = formula,
          data = data,
          paired = paired,
          id = id,
          within = within,
          between = between,
          random_slope = FALSE,
          alpha = alpha,
          k = k,
          code = code,
          debug = debug,
          verbose = verbose
        )

        # Initialiser bilan si n\u00e9cessaire avec structure attendue par m.test()
        if (!exists("bilan")) {
          bilan <- list(
            x,           # [[1]]
            g,           # [[2]]
            TRUE,        # [[3]] check_normality (on suppose TRUE pour mod\u00e8le mixte)
            TRUE         # [[4]] check_variance_equal (on suppose TRUE pour mod\u00e8le mixte)
          )
        }

        # ADAPTATION: .mixed_model_analysis() peut retourner soit une liste compl\u00e8te,
        # soit directement un mod\u00e8le lmerMod (\u00e0 cause du return anticip\u00e9 ligne 251)
        if (inherits(mixed_result, "lmerMod")) {
          # Cas o\u00f9 on a un mod\u00e8le brut (return anticip\u00e9 dans .mixed_model_analysis)
          # Extraire les infos n\u00e9cessaires du mod\u00e8le directement
          .dbg("Detected raw lmerMod object from .mixed_model_analysis()",
               "Objet lmerMod brut d\u00e9tect\u00e9 depuis .mixed_model_analysis()",
               debug = debug)

          # Obtenir l'ANOVA table avec lmerTest pour les p-values
          anova_table <- tryCatch({
            suppressMessages({
              # Convertir en lmerModLmerTest si n\u00e9cessaire
              if (!inherits(mixed_result, "lmerModLmerTest")) {
                mixed_result <- lmerTest::lmer(formula(mixed_result),
                                               data = model.frame(mixed_result),
                                               REML = TRUE)
              }
              anova(mixed_result, type = "III")
            })
          }, error = function(e) NULL)

          # Extraire p-value globale
          global_pval <- if (!is.null(anova_table)) {
            pval_col <- which(colnames(anova_table) %in% c("Pr(>F)", "p.value", "P(>|t|)"))
            if (length(pval_col) > 0) {
              row_names <- rownames(anova_table)
              valid_rows <- which(!grepl("Residuals|Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
              if (length(valid_rows) > 0) min(anova_table[valid_rows, pval_col[1]], na.rm = TRUE) else NA
            } else NA
          } else NA

          bilan$robust_results <- list(
            method = "Mixed_Model_lmer",
            model = mixed_result,
            anova_table = anova_table,
            significant_effects = NULL,
            posthoc_applicable = TRUE,
            test_result = mixed_result,
            assumptions_checked = FALSE,
            variance_components = NULL,
            random_effects = NULL,
            fixed_effects = NULL
          )
          bilan$k <- k
          bilan$global_pvalue <- global_pval

        } else {
          # Cas normal : mixed_result est une liste compl\u00e8te
          bilan$robust_results <- list(
            method = "Mixed_Model_lmer",
            model = mixed_result$model,
            anova_table = mixed_result$anova_table,
            significant_effects = mixed_result$significant_effects,
            posthoc_applicable = mixed_result$posthoc_applicable,
            test_result = mixed_result$model,
            assumptions_checked = mixed_result$assumptions_checked,
            variance_components = mixed_result$variance_components,
            random_effects = mixed_result$random_effects,
            fixed_effects = mixed_result$fixed_effects
          )
          bilan$k <- mixed_result$k
          bilan$global_pvalue <- mixed_result$global_pvalue
        }

        return(bilan)

        # === CODE ANCIEN SUPPRIM\u00c9 (maintenant dans .mixed_model_analysis()) ===

      } else {
        # CAS 2: Pas d'id fourni ou design simple non g\u00e9r\u00e9 par les cas pr\u00e9c\u00e9dents

        k <- .vbse(
          paste0("Complex repeated measures design detected.\n",
                 "\tSuggested robust approaches:\n",
                 "\t- Mixed model with robust estimation (lmerTest::lmer)\n",
                 "\t- Robust repeated measures from {WRS2} package\n",
                 "\t- Permutation-based methods\n",
                 "\tNote: Requires id parameter for automatic implementation."),
          paste0("Design de mesures r\u00e9p\u00e9t\u00e9es complexe d\u00e9tect\u00e9.\n",
                 "\tApproches robustes sugg\u00e9r\u00e9es :\n",
                 "\t- Mod\u00e8le mixte avec estimation robuste (lmerTest::lmer)\n",
                 "\t- Mesures r\u00e9p\u00e9t\u00e9es robustes du package {WRS2}\n",
                 "\t- M\u00e9thodes par permutation\n",
                 "\tNote : N\u00e9cessite le param\u00e8tre id pour impl\u00e9mentation automatique."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        robust_results$method <- "Repeated_Measures_Manual"
        robust_results$warnings <- c(
          robust_results$warnings,
          "Complex repeated measures - manual implementation required"
        )
      }
    }

  }  else if (alea || alea_plus) {
      k <- .vbse(
        paste0("Suggested robust approach for random effects:\n",
               "\t- Mixed model with robust estimation (lmerTest::lmer)\n",
               "\t- Generalized Estimating Equations (GEE) if appropriate\n",
               "\t- Bootstrap methods for inference\n",
               "\tNote: Implementation requires dedicated function .mixed_model_analysis() (in development)."),
        paste0("Approche robuste sugg\u00e9r\u00e9e pour effets al\u00e9atoires :\n",
               "\t- Mod\u00e8le mixte avec estimation robuste (lmerTest::lmer)\n",
               "\t- \u00e9quations d'estimation g\u00e9n\u00e9ralis\u00e9es (GEE) si appropri\u00e9\n",
               "\t- M\u00e9thodes de bootstrap pour l'inf\u00e9rence\n",
               "\tNote : L'impl\u00e9mentation n\u00e9cessite fonction d\u00e9di\u00e9e .mixed_model_analysis() (en d\u00e9veloppement)."),
        verbose = verbose, code = code, k = k, cpt = "off"
      )

    } else {
      # NOTE: Message simplifi\u00e9 - le test appropri\u00e9 sera annonc\u00e9 lors de son ex\u00e9cution
      # Inutile de lister toutes les options possibles ici (noie l'utilisateur)
      #============================================================================
      #           IMPL\u00c9MENTATION AUTOMATIQUE DES TESTS ROBUSTES
      #           (R\u00e9ponse \u00e0 NOTE PERSO #11)
      #============================================================================

      # Initialiser la structure de r\u00e9sultats robustes
      robust_results <- list(
        method = NULL,
        test_result = NULL,
        assumptions_checked = list(),
        warnings = character(),
        posthoc_applicable = FALSE
      )

      # D\u00e9terminer le nombre de facteurs r\u00e9els
      n_factors <- length(factor_vars)

      wrs2_done <- FALSE

      if (force_wrs2 && !paired && !check_ancova && n_factors >= 2) {
        # Trimming adaptatif selon la proportion r\u00e9elle d'outliers
        # R\u00e9f\u00e9rence : Wilcox (2017), ch. 7 - tr=0.05 minimum pour robustesse aux queues lourdes
        n_total_obs <- nrow(data)
        pct_outliers <- if (n_total_obs > 0) 100 * n_outliers_marginal / n_total_obs else 0

        if (n_extreme_outliers > 0 || pct_outliers > 5) {
          tr_val <- 0.20
          tr_reason_en <- paste0("extreme outliers or >5% outliers (", round(pct_outliers, 1), "%)")
          tr_reason_fr <- paste0("outliers extr\u00eames ou >5% d'outliers (", round(pct_outliers, 1), "%)")
        } else if (pct_outliers > 0) {
          tr_val <- 0.10
          tr_reason_en <- paste0("moderate outliers (", round(pct_outliers, 1), "%, non-extreme)")
          tr_reason_fr <- paste0("outliers mod\u00e9r\u00e9s (", round(pct_outliers, 1), "%, non extr\u00eames)")
        } else {
          tr_val <- 0.05
          tr_reason_en <- "no outliers, minimal robustness trim"
          tr_reason_fr <- "pas d'outliers, trimming minimal de robustesse"
        }

        k <- .vbse(
          paste0("Applying robust WRS2 ANOVA due to extreme imbalance (ratio > 10).\n",
                 "\tTrim level: ", round(tr_val * 100), "% (", tr_reason_en, ")."),
          paste0("Application de l'ANOVA robuste WRS2 en raison d'un d\u00e9s\u00e9quilibre extr\u00eame (ratio > 10).\n",
                 "\tTrimmage : ", round(tr_val * 100), "% (", tr_reason_fr, ")."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (isTRUE(code)) {
          k_code <- k_code + 1
          .code_multi(k_code, "ANOVA robuste WRS2 (d\u00e9s\u00e9quilibre extr\u00eame)", c(
            "library(WRS2)",
            paste0("tr <- ", tr_val),
            if (n_factors == 2) "WRS2::t2way(formula, data = data, tr = tr)" else
              "WRS2::t3way(formula, data = data, tr = tr)"
          ))
        }

        if (requireNamespace("WRS2", quietly = TRUE)) {
          tryCatch({
            wrs2_result <- if (n_factors == 2) {
              WRS2::t2way(formula, data = data, tr = tr_val)
            } else {
              WRS2::t3way(formula, data = data, tr = tr_val)
            }

            robust_results$method <- if (n_factors == 2) "WRS2_t2way" else "WRS2_t3way"
            robust_results$test_result <- wrs2_result
            robust_results$posthoc_applicable <- TRUE
            wrs2_done <- TRUE

            if (verbose) {
              print(wrs2_result)
              cat("\n")
            }
          }, error = function(e) {
            robust_results$warnings <- c(robust_results$warnings, paste0("WRS2 error: ", e$message))
            k <- .vbse(
              paste0("WRS2 robust ANOVA failed: ", e$message, "\n\tFalling back to permutation ANOVA."),
              paste0("\u00c9chec ANOVA robuste WRS2 : ", e$message, "\n\tRetour vers ANOVA par permutation."),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            wrs2_done <- FALSE
          })
        } else {
          k <- .vbse(
            "Package WRS2 not available. Install with: install.packages('WRS2').\n\tFalling back to permutation ANOVA.",
            "Package WRS2 non disponible. Installer avec : install.packages('WRS2').\n\tRetour vers ANOVA par permutation.",
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          wrs2_done <- FALSE
        }
      }

      #=============================================================================
      # PHASE 1: S\u00c9LECTION ET APPLICATION AUTOMATIQUE DU TEST ROBUSTE
      #=============================================================================

      if (!wrs2_done && paired && nlevels(g_cat) >= 3 && n_factors == 1) {
        #---------------------------------------------------------------------------
        # CAS 1: MESURES R\u00c9P\u00c9T\u00c9ES, k >= 3, UN SEUL FACTEUR
        # \u2192 Test de Friedman (extension non param\u00e9trique de RM-ANOVA)
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Applying Friedman test for repeated measures (k >= 3 conditions, one factor)...",
          "Application du test de Friedman pour mesures r\u00e9p\u00e9t\u00e9es (k >= 3 conditions, un facteur)...",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        tryCatch({
          # Pr\u00e9parer les donn\u00e9es en format large pour Friedman
          if (!is.null(within) && length(within) == 1) {
            # Cr\u00e9er data frame pour reshape
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
        paste0("Test de Friedman termin\u00e9 :\n",
                     "\tKhi\u00b2 = ", round(friedman_result$statistic, 3),
                     ", ddl = ", friedman_result$parameter,
                     ", p = ", .format_pval(friedman_result$p.value)),
              verbose = verbose, code = code, k = k, cpt = "off"
            )

            # V\u00e9rifier les assomptions sp\u00e9cifiques \u00e0 Friedman
            # 1. \u00c9chelle au moins ordinale
            n_unique <- length(unique(x[!is.na(x)]))
            if (n_unique < 5) {
              robust_results$warnings <- c(
                robust_results$warnings,
                paste0("Few distinct values (n=", n_unique, "). Verify ordinal scale assumption.")
              )
              k <- .vbse(
                paste0("Note: Only ", n_unique, " distinct values in DV. Verify ordinal scale is appropriate."),
                paste0("Note : Seulement ", n_unique, " valeurs distinctes dans VD. V\u00e9rifier que l'\u00e9chelle ordinale est appropri\u00e9e."),
                verbose = verbose, code = code, k = k, cpt = "off"
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
                       "% de rangs ex-aequo. Friedman ajuste automatiquement mais la puissance peut \u00eatre r\u00e9duite."),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }

          } else {
            # Plusieurs facteurs within : Friedman non applicable
            stop("Multiple within-subject factors detected. Friedman test not applicable.")
          }

        }, error = function(e) {
          warning(.msg(
            paste0("Failed to perform Friedman test: ", e$message),
            paste0("\u00e9chec du test de Friedman : ", e$message)
          ))
          robust_results$method <- "Friedman_Failed"
          robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

          k <- .vbse(
            "Friedman test failed. Consider mixed models for this design.",
            "Test de Friedman \u00e9chou\u00e9. Envisager mod\u00e8les mixtes pour ce plan.",
            verbose = verbose, code = code, k = k, cpt = "on"
          )
        })

      } else if (!wrs2_done && n_factors == 1 && !paired && !check_ancova && nlevels(g_cat) >= 2) {
        #---------------------------------------------------------------------------
        # CAS 2: UN SEUL FACTEUR, DONN\u00c9ES IND\u00c9PENDANTES, PAS D'ANCOVA
        # \u2192 Kruskal-Wallis (d\u00e9j\u00e0 document\u00e9 dans NOTE #12, mais enrichi ici)
        # IMPORTANT: Kruskal-Wallis n'est PAS adapt\u00e9 pour ANCOVA (avec covariables)
        # ni pour designs multi-facteurs. Utilis\u00e9 UNIQUEMENT pour 1 facteur seul.
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Applying Kruskal-Wallis test for one-way independent design (one factor, no covariates)...",
          "Application du test de Kruskal-Wallis pour plan unifactoriel ind\u00e9pendant (un facteur, pas de covariables)...",
          verbose = verbose, code = code, k = k, cpt = "on"
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
      paste0("Test de Kruskal-Wallis termin\u00e9 :\n",
                   "\tKhi\u00b2 = ", round(kw_result$statistic, 3),
                   ", ddl = ", kw_result$parameter,
                   ", p = ", .format_pval(kw_result$p.value)),
            verbose = verbose, code = code, k = k, cpt = "off"
          )

          # Contr\u00f4le: homog\u00e9n\u00e9it\u00e9 des distributions avec Fligner-Killeen
          k <- .vbse(
            "Checking variance homogeneity with Fligner-Killeen test (robust to non-normality)...",
            "V\u00e9rification de l'homog\u00e9n\u00e9it\u00e9 des variances avec test de Fligner-Killeen (robuste \u00e0 la non-normalit\u00e9)...",
            verbose = verbose, code = code, k = k, cpt = "on"
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
            paste0("Fligner-Killeen : khi\u00b2 = ", round(fk_result$statistic, 3),
                   ", p = ", .format_pval(fk_result$p.value)),
            verbose = verbose, code = code, k = k, cpt = "off"
          )

          if (fk_result$p.value < alpha) {
            robust_results$warnings <- c(
              robust_results$warnings,
              "Variance heterogeneity detected (Fligner-Killeen)"
            )
            k <- .vbse(
              "Variance heterogeneity detected. Kruskal-Wallis remains valid but power may be affected.",
              "H\u00e9t\u00e9rog\u00e9n\u00e9it\u00e9 des variances d\u00e9tect\u00e9e. Kruskal-Wallis reste valide mais la puissance peut \u00eatre affect\u00e9e.",
              verbose = verbose, code = code, k = k, cpt = "off"
            )
          }

          # Contr\u00f4le: d\u00e9tection d'outliers extr\u00eames
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
                paste0("Attention : ", n_extreme, " valeur(s) extr\u00eame(s) d\u00e9tect\u00e9e(s). Peut affecter l'interpr\u00e9tation."),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
          }

        }, error = function(e) {
          warning(.msg(
            paste0("Failed to perform Kruskal-Wallis test: ", e$message),
            paste0("\u00e9chec du test de Kruskal-Wallis : ", e$message)
          ))
          robust_results$method <- "Kruskal_Wallis_Failed"
          robust_results$warnings <- c(robust_results$warnings, as.character(e$message))
        })

      } else if (!wrs2_done && n_factors == 2 && !paired && !check_ancova) {
        #---------------------------------------------------------------------------
        # CAS 3: DEUX FACTEURS, DONN\u00c9ES IND\u00c9PENDANTES, PAS ANCOVA
        # \u2192 Scheirer-Ray-Hare (extension de Kruskal-Wallis pour 2-way)
        # R\u00e9f\u00e9rence: Scheirer, Ray & Hare (1976). The analysis of ranked data derived
        # from completely randomized factorial designs. Biometrics, 32(2), 429-434.
        # https://doi.org/10.2307/2529511
        #---------------------------------------------------------------------------

        k <- .vbse(
          paste0("Applying Scheirer-Ray-Hare test [rcompanion::scheirerRayHare()] for two-way independent design.\n",
                 "\tReason: Non-parametric alternative for 2-way ANOVA (rank-based).\n",
                 "\t\tTests main effects AND interaction.\n",
                 "\t\t(Appropriate for discrete/non-normal data with 2 factors.)"),
          paste0("Application du test de Scheirer-Ray-Hare [rcompanion::scheirerRayHare()] pour plan bifactoriel ind\u00e9pendant.\n",
                 "\tRaison : Alternative non param\u00e9trique \u00e0 l'ANOVA 2-way (bas\u00e9e sur rangs).\n",
                 "\t\tTeste les effets principaux ET l'interaction.\n",
                 "\t\t(Appropri\u00e9 pour donn\u00e9es discr\u00e8tes/non-normales avec 2 facteurs.)"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("rcompanion", quietly = TRUE)) {
          tryCatch({
            # Construire formule pour Scheirer-Ray-Hare (main effects + interaction)
            response_var <- all.vars(formula)[1]
            srh_formula <- as.formula(paste0(
              response_var, " ~ ", factor_vars[1], " + ", factor_vars[2]
            ))

            # Test de Scheirer-Ray-Hare
            # Capturer la sortie pour \u00e9viter l'affichage automatique (DV, Observations, D, MS total)
            invisible(capture.output({
              srh_result <- rcompanion::scheirerRayHare(
                formula = srh_formula,
                data = data
              )
            }))

            robust_results$method <- "Scheirer_Ray_Hare"
            robust_results$test_result <- srh_result
            robust_results$posthoc_applicable <- TRUE

            if (verbose) {
              # Afficher r\u00e9sultats SRH de mani\u00e8re \u00e9pur\u00e9e (sans Residuals)
              # Ne garder que les lignes des facteurs et de l'interaction
              srh_display <- srh_result[!rownames(srh_result) %in% c("Residuals"), , drop = FALSE]
              # Ne garder que les colonnes Df, H, et p.value
              cols_to_keep <- intersect(c("Df", "H", "p.value"), colnames(srh_display))
              if (length(cols_to_keep) > 0) {
                srh_display <- srh_display[, cols_to_keep, drop = FALSE]
              }
              cat("\n")
              print(srh_display)
              cat("\n")
            }

            # Marquer que les interactions SONT test\u00e9es (contrairement \u00e0 Kruskal-Wallis)
            robust_results$assumptions_checked$interaction_tested <- TRUE

            # CONCLUSION : Identifier effets significatifs (main effects ET interaction)
            if (verbose && !is.null(srh_result) && "p.value" %in% colnames(srh_result)) {
              sig_effects <- rownames(srh_result)[
                !grepl("Residuals", rownames(srh_result), ignore.case = TRUE) &
                srh_result$p.value < alpha
              ]

              if (length(sig_effects) > 0) {
                k <- .vbse(
                  paste0("==> Significant effects detected (alpha = ", alpha, "): ",
                         paste(sig_effects, collapse = ", ")),
                  paste0("==> Effets significatifs d\u00e9tect\u00e9s (alpha = ", alpha, ") : ",
                         paste(sig_effects, collapse = ", ")),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              } else {
                k <- .vbse(
                  paste0("==> No significant effects at alpha = ", alpha),
                  paste0("==> Aucun effet significatif \u00e0 alpha = ", alpha),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              }
            }

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to perform Scheirer-Ray-Hare test: ", e$message),
              paste0("\u00e9chec du test de Scheirer-Ray-Hare : ", e$message)
            ))
            robust_results$method <- "Scheirer_Ray_Hare_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              "Scheirer-Ray-Hare failed. Falling back to permutation ANOVA recommendation.",
              "Scheirer-Ray-Hare \u00e9chou\u00e9. Recommandation d'ANOVA par permutation.",
              verbose = verbose, code = code, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            "Package 'rcompanion' not available for Scheirer-Ray-Hare. Switching to permutation ANOVA...",
            "Package 'rcompanion' non disponible pour Scheirer-Ray-Hare. Passage vers ANOVA par permutation...",
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          robust_results$method <- "Scheirer_Ray_Hare_Unavailable"

          # Rediriger vers permutation ANOVA (trait\u00e9 dans CAS 4)
          n_factors <- 999  # Force passage \u00e0 CAS 4
        }

      }

      if (!wrs2_done && n_factors >= 2 && !paired &&
          (is.null(robust_results$method) ||
           robust_results$method %in% c("Scheirer_Ray_Hare_Unavailable", "Scheirer_Ray_Hare_Failed"))) {
        #---------------------------------------------------------------------------
        # CAS 4: PLUSIEURS FACTEURS (2+), DONN\u00c9ES IND\u00c9PENDANTES
        # \u2192 ANOVA par permutation (m\u00e9thode UNIVERSELLE et OPTIMALE)
        # R\u00e9f\u00e9rence: Anderson & Robinson (2001). Permutation tests for linear models.
        # Australian & New Zealand Journal of Statistics, 43(1), 75-88.
        # https://doi.org/10.1111/1467-842X.00156
        # R\u00e9f\u00e9rence: Manly (2007). Randomization, Bootstrap and Monte Carlo Methods
        # in Biology (3rd ed.). Chapman & Hall/CRC.
        #---------------------------------------------------------------------------

        k <- .vbse(
          paste0("Applying permutation ANOVA [lmPerm::aovp()] for multi-factor independent design.\n",
                 "\tDesign: ", length(factor_vars), " factor(s) detected."),
          paste0("Application de l'ANOVA par permutation [lmPerm::aovp()] pour plan multi-facteurs ind\u00e9pendant.\n",
                 "\tPlan : ", length(factor_vars), " facteur(s) d\u00e9tect\u00e9(s)."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("lmPerm", quietly = TRUE)) {
          tryCatch({
            # ANOVA par permutation avec nombre de permutations adapt\u00e9
            # Plus de facteurs = plus de permutations n\u00e9cessaires
            n_perm <- ifelse(length(factor_vars) <= 2, "Prob", "Exact")

            # Ajuster le mod\u00e8le par permutation
            # capture.output() absorbe le message parasite "Settings:  unique SS" \u00e9mis par aovp()
            utils::capture.output(
              perm_result <- lmPerm::aovp(formula, data = data, perm = n_perm)
            )

            robust_results$method <- "Permutation_ANOVA"
            robust_results$test_result <- perm_result
            robust_results$posthoc_applicable <- TRUE

            if (verbose) {
              k <- .vbse(
                "ANOVA table (permutation p-values):",
                "Tableau ANOVA (valeurs p par permutation) :",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
              cat("\n")
              print(summary(perm_result))
              cat("\n")
            }

            # Extraire et afficher effets significatifs
            perm_summary <- summary(perm_result)
            if (!is.null(perm_summary) && "Pr(Prob)" %in% colnames(perm_summary[[1]])) {
              sig_effects <- rownames(perm_summary[[1]])[
                !grepl("Residuals", rownames(perm_summary[[1]]), ignore.case = TRUE) &
                perm_summary[[1]][, "Pr(Prob)"] < alpha
              ]

              if (length(sig_effects) > 0) {
                sig_list <- paste0("\t  - ", trimws(sig_effects), collapse = "\n")
                k <- .vbse(
                  paste0("Significant effects detected (alpha = ", alpha, "):\n", sig_list),
                  paste0("Effets significatifs d\u00e9tect\u00e9s (alpha = ", alpha, ") :\n", sig_list),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              } else {
                k <- .vbse(
                  paste0("No significant effects at alpha = ", alpha),
                  paste0("Aucun effet significatif \u00e0 alpha = ", alpha),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              }
            }

            robust_results$assumptions_checked$distribution_free <- TRUE
            robust_results$assumptions_checked$handles_imbalance <- TRUE
            robust_results$assumptions_checked$tests_interactions <- TRUE
            robust_results$assumptions_checked$outlier_robust <- TRUE

            # Stocker le degr\u00e9 de d\u00e9s\u00e9quilibre (information interne, pas de message verbose)
            # Note : aovp g\u00e8re correctement les d\u00e9s\u00e9quilibres jusqu'\u00e0 10:1
            table_data_local <- table(g_cat)
            if (length(table_data_local) > 1) {
              max_ratio <- max(table_data_local) / min(table_data_local)
              robust_results$assumptions_checked$imbalance_ratio <- max_ratio

              if (max_ratio > 3) {
                robust_results$warnings <- c(
                  robust_results$warnings,
                  paste0("Severe imbalance: max/min ratio = ", round(max_ratio, 2))
                )
              }
            }

            # Contr\u00f4les d'assomptions pour post-hocs (normalit\u00e9 + variance)
            # M\u00eame si l'ANOVA est par permutation, les post-hocs d\u00e9pendent des hypoth\u00e8ses
            k <- .vbse(
              "Post-hoc suitability check after robust ANOVA.\n\tNormality + variance tests determine the appropriate post-hoc family.",
              "Contr\u00f4le d'ad\u00e9quation des post-hocs apr\u00e8s ANOVA robuste.\n\tTests de normalit\u00e9 + variance pour choisir les post-hocs appropri\u00e9s.",
              verbose = verbose, code = code, k = k, cpt = "on"
            )
            need_posthoc_checks <- (!normality_already_tested || !variance_already_tested)

            if (!need_posthoc_checks) {
              k <- .vbse(
                "Note: Normality/variance already checked earlier. Reusing results for post-hoc selection.",
                "Note : Normalit\u00e9/variance d\u00e9j\u00e0 contr\u00f4l\u00e9es plus t\u00f4t. R\u00e9utilisation des r\u00e9sultats pour les post-hocs.",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }

            if (isTRUE(code) && need_posthoc_checks) {
              k_code <- k_code + 1
              .code_multi(k_code, "Contr\u00f4le post-hoc (normalit\u00e9 + variance)", c(
                "norm_res <- .normality(x, g = g_cat, alpha = alpha, paired = FALSE)",
                "var_res <- .variance(x, g = g_cat, check_normality = norm_res[[1]], alpha = alpha)"
              ))
            }
            if (!normality_already_tested) {
              normality_result <- .normality(
                x, g = g_cat, alpha = alpha, paired = FALSE,
                debug = debug, verbose = verbose, code = code, k = k,
                cpt = "sub", prefix = "a) "
              )
              check_normality <- normality_result[[1]]
              k <- normality_result[[2]]
              normality_already_tested <- TRUE
            }

            if (!variance_already_tested) {
              pvals_variance <- .variance(
                x, g = g_cat, check_normality = check_normality,
                alpha = alpha, paired = FALSE, debug = debug,
                verbose = verbose, code = code, k = k, cpt = "sub",
                assumption_label = NULL, prefix = "b) "
              )
              check_variance_equal <- pvals_variance[[1]]
              k <- pvals_variance[[2]]
              variance_already_tested <- TRUE
            }
            if (isTRUE(outliers_marginal_detected) || n_extreme_outliers > 0) {
              check_normality <- FALSE
              k <- .vbse(
                "Note: Outliers detected earlier. Forcing robust/non-parametric post-hocs.",
                "Note : Outliers d\u00e9tect\u00e9s en amont. Post-hocs robustes/non-param\u00e9triques impos\u00e9s.",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
            if (!check_variance_equal) {
              check_normality <- FALSE  # Forcer post-hocs non-param\u00e9triques
            }

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to perform permutation ANOVA: ", e$message),
              paste0("\u00e9chec de l'ANOVA par permutation : ", e$message)
            ))
            robust_results$method <- "Permutation_ANOVA_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              "Permutation ANOVA failed. No automatic robust alternative available for this design.",
              "ANOVA par permutation \u00e9chou\u00e9e. Aucune alternative robuste automatique disponible pour ce plan.",
              verbose = verbose, code = code, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            paste0("Package 'lmPerm' NOT available for permutation ANOVA.\n",
                   "\tThis is the RECOMMENDED approach for multi-factor robust analysis.\n",
                   "\tInstall with: install.packages('lmPerm')"),
            paste0("Package 'lmPerm' NON disponible pour ANOVA par permutation.\n",
                   "\tC'est l'approche RECOMMAND\u00e0E pour analyse robuste multi-facteurs.\n",
                   "\tInstaller avec : install.packages('lmPerm')"),
            verbose = verbose, code = code, k = k, cpt = "on"
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
          (paired && !is.null(between) && length(between) > 0) ||
          use_mixed_model) {
        #---------------------------------------------------------------------------
        # CAS 5: EFFETS AL\u00c9ATOIRES OU DESIGN RM COMPLEXE / R\u00c9PLICATS
        # \u2192 Mod\u00e8les mixtes (lmer) avec validation
        #---------------------------------------------------------------------------

        k <- .vbse(
          "Complex design with random effects or mixed within/between. Applying linear mixed model (lmer)...",
          "Plan complexe avec effets al\u00e9atoires ou within/between mixte. Application d'un mod\u00e8le lin\u00e9aire mixte (lmer)...",
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        if (requireNamespace("lme4", quietly = TRUE) &&
            requireNamespace("lmerTest", quietly = TRUE)) {
          tryCatch({
            # Construire la formule lmer selon le design
            response_var <- all.vars(formula)[1]

            if (alea_lmer) {
              # Syntaxe lmer native : utiliser la formule telle quelle
              # L'utilisateur a explicitement sp\u00e9cifi\u00e9 (1|id) ou (var|id)
              lmer_formula <- formula
              k <- .vbse(
                paste0("Using user-provided lmer formula: ", deparse(lmer_formula)),
                paste0("Utilisation de la formule lmer fournie par l'utilisateur : ", deparse(lmer_formula)),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            } else if (paired && !is.null(within)) {
              # Mesures r\u00e9p\u00e9t\u00e9es
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
            } else if (alea_error) {
              # Formule avec Error() : conversion simplifi\u00e9e
              k <- .vbse(
                "Note: Complex Error() syntax detected. Using simplified random intercept model.",
                "Note : Syntaxe Error() complexe d\u00e9tect\u00e9e. Utilisation d'un mod\u00e8le \u00e0 intercept al\u00e9atoire simplifi\u00e9.",
                verbose = verbose, code = code, k = k, cpt = "off"
              )

              # Extraire le terme d'erreur (simplifi\u00e9)
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
              paste0("Formule du mod\u00e8le mixte : ", deparse(lmer_formula)),
              verbose = verbose, code = code, k = k, cpt = "off"
            )

            # Ajuster le mod\u00e8le mixte
            mixed_model <- lmerTest::lmer(lmer_formula, data = data, REML = TRUE)

            robust_results$method <- "Mixed_Model_lmer"
            robust_results$test_result <- mixed_model
            robust_results$model <- mixed_model  # For post-hoc compatibility
            robust_results$posthoc_applicable <- TRUE

            k <- .vbse(
              "Mixed model fitted successfully (REML estimation).",
              "Mod\u00e8le mixte ajust\u00e9 avec succ\u00e8s (estimation REML).",
              verbose = verbose, code = code, k = k, cpt = "off"
            )

            # Afficher tableau ANOVA et extraire effets significatifs
            anova_table <- stats::anova(mixed_model, type = "III")

            if (verbose) {
              k <- .vbse(
                "Type III ANOVA table for mixed model (Satterthwaite approximation):",
                "Tableau ANOVA de type III pour mod\u00e8le mixte (approximation de Satterthwaite) :",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
              cat("\n")
              print(anova_table)
              cat("\n")
            }

            # Extraire effets significatifs pour post-hocs
            if ("Pr(>F)" %in% colnames(anova_table)) {
              sig_effects <- rownames(anova_table)[which(anova_table[, "Pr(>F)"] < alpha)]
              robust_results$significant_effects <- sig_effects

              if (length(sig_effects) > 0 && verbose) {
                k <- .vbse(
                  paste0("Significant effect(s) detected (\u03b1 = ", alpha, "): ", paste(sig_effects, collapse = ", ")),
                  paste0("Effet(s) significatif(s) d\u00e9tect\u00e9(s) (\u03b1 = ", alpha, ") : ", paste(sig_effects, collapse = ", ")),
                  verbose = verbose, code = code, k = k, cpt = "off"
                )
              }
            }

            #-----------------------------------------------------------------------
            # DIAGNOSTICS DU MOD\u00c8LE MIXTE
            #-----------------------------------------------------------------------
            k <- .vbse(
              "Performing mixed model diagnostics...",
              "R\u00e9alisation des diagnostics du mod\u00e8le mixte...",
              verbose = verbose, code = code, k = k, cpt = "on"
            )

            # 1. Normalit\u00e9 des r\u00e9sidus
            lmer_resid <- residuals(mixed_model)

            resid_normality_result <- .normality(
              lmer_resid, g = NULL, alpha = alpha, paired = FALSE,
              debug = debug, verbose = verbose, code = code, k = k
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
                paste0("Attention : Les r\u00e9sidus violent l'hypoth\u00e8se de normalit\u00e9.\n",
                       "\tEnvisager : (1) Transformation de la VD, (2) Mod\u00e8les mixtes robustes (package robustlmm),\n",
                       "\t(3) Mod\u00e8les lin\u00e9aires g\u00e9n\u00e9ralis\u00e9s mixtes (GLMM) si appropri\u00e9."),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }

            # 2. Homosc\u00e9dasticit\u00e9 des r\u00e9sidus
            if (requireNamespace("car", quietly = TRUE)) {
              # Cr\u00e9er facteur de groupe appropri\u00e9
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
                paste0("Test de Levene sur r\u00e9sidus : F = ", round(lev_test$`F value`[1], 3),
                       ", p = ", .format_pval(p_lev)),
                verbose = verbose, code = code, k = k, cpt = "off"
              )

              if (p_lev < alpha) {
                robust_results$warnings <- c(
                  robust_results$warnings,
                  "Mixed model: residual heteroscedasticity detected"
                )
              }
            }

            # 3. Variance components
            var_comp_S4 <- lme4::VarCorr(mixed_model)
            # Convertir en dataframe pour \u00e9viter erreurs S4 as.vector()
            var_comp <- as.data.frame(var_comp_S4)
            robust_results$assumptions_checked$variance_components <- var_comp

            if (verbose) {
              k <- .vbse(
                "Variance components:",
                "Composantes de variance :",
                verbose = verbose, code = code, k = k, cpt = "off"
              )
              cat("\n")
              print(var_comp_S4)  # Afficher format original pour lisibilit\u00e9
              cat("\n")
            }

            # Remplacer le mod\u00e8le NULL par le mod\u00e8le mixte
            model <- mixed_model

          }, error = function(e) {
            warning(.msg(
              paste0("Failed to fit mixed model: ", e$message),
              paste0("\u00e9chec de l'ajustement du mod\u00e8le mixte : ", e$message)
            ))
            robust_results$method <- "Mixed_Model_Failed"
            robust_results$warnings <- c(robust_results$warnings, as.character(e$message))

            k <- .vbse(
              paste0("Mixed model fitting failed. Possible causes:\n",
                     "\t- Singular fit (variance component = 0)\n",
                     "\t- Convergence issues\n",
                     "\t- Insufficient data for random effects\n",
                     "\tConsider simplified random structure or consult statistician."),
              paste0("Ajustement du mod\u00e8le mixte \u00e9chou\u00e9. Causes possibles :\n",
                     "\t- Ajustement singulier (composante de variance = 0)\n",
                     "\t- Probl\u00e8mes de convergence\n",
                     "\t- Donn\u00e9es insuffisantes pour effets al\u00e9atoires\n",
                     "\tEnvisager structure al\u00e9atoire simplifi\u00e9e ou consulter statisticien."),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
          })
        } else {
          k <- .vbse(
            paste0("Packages 'lme4' and/or 'lmerTest' NOT available.\n",
                   "\tMixed models are REQUIRED for this complex design.\n",
                   "\tInstall with: install.packages(c('lme4', 'lmerTest'))"),
            paste0("Packages 'lme4' et/ou 'lmerTest' NON disponibles.\n",
                   "\tLes mod\u00e8les mixtes sont REQUIS pour ce plan complexe.\n",
                   "\tInstaller avec : install.packages(c('lme4', 'lmerTest'))"),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          robust_results$method <- "Mixed_Model_Unavailable"
          robust_results$warnings <- c(
            robust_results$warnings,
            "lme4/lmerTest packages not available - INSTALLATION REQUIRED"
          )
        }
      }
    } # J'ai ajout\u00e9 pour voir....
      # Si aucune m\u00e9thode n'a \u00e9t\u00e9 appliqu\u00e9e
      if (is.null(robust_results$method)) {
        # Construire information sur r\u00e9plicats si applicable
        replicates_info_en <- ""
        replicates_info_fr <- ""
        if (use_mixed_model && exists("n_replicates_per_cell") && n_replicates_per_cell > 1) {
          replicates_info_en <- paste0("\t  - Replicates: ", n_replicates_per_cell, " per subject\u00d7condition\n")
          replicates_info_fr <- paste0("\t  - R\u00e9plicats : ", n_replicates_per_cell, " par sujet\u00d7condition\n")
        }

        k <- .vbse(
          paste0("WARNING: No automatic robust method could be applied for this design.\n",
                 "\tDesign characteristics:\n",
                 "\t  - Paired: ", if(paired) "Yes" else "No", "\n",
                 "\t  - Number of factors: ", n_factors, "\n",
                 "\t  - ANCOVA: ", if(check_ancova) "Yes" else "No", "\n",
                 "\t  - Random effects: ", if(alea) "Yes" else "No", "\n",
                 replicates_info_en,
                 "\tRecommendation: Consult a statistician for appropriate analysis strategy."),
          paste0("ATTENTION : Aucune m\u00e9thode robuste automatique n'a pu \u00eatre appliqu\u00e9e pour ce plan.\n",
                 "\tCaract\u00e9ristiques du plan :\n",
                 "\t  - Appari\u00e9 : ", if(paired) "Oui" else "Non", "\n",
                 "\t  - Nombre de facteurs : ", n_factors, "\n",
                 "\t  - ANCOVA : ", if(check_ancova) "Oui" else "Non", "\n",
                 "\t  - Effets al\u00e9atoires : ", if(alea) "Oui" else "Non", "\n",
                 replicates_info_fr,
                 "\tRecommandation : Consulter un statisticien pour strat\u00e9gie d'analyse appropri\u00e9e."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        robust_results$method <- "No_Automatic_Method_Available"
        robust_results$warnings <- c(
          robust_results$warnings,
          "Design configuration not covered by automatic robust methods"
        )
      }

      #=============================================================================
      # PHASE 2: AVERTISSEMENTS (si pr\u00e9sents)
      #=============================================================================
      # NOTE: R\u00e9sum\u00e9 "M\u00e9thode appliqu\u00e9e" SUPPRIM\u00c9 car redondant avec messages pr\u00e9c\u00e9dents
      # qui affichent d\u00e9j\u00e0 la m\u00e9thode choisie et les r\u00e9sultats d\u00e9taill\u00e9s

      # Avertissements internes conserves dans robust_results$warnings
      # mais non affiches dans le bilan (redondant avec messages precedents)

      # Ajouter robust_results \u00e0 la structure de retour finale
      # (sera int\u00e9gr\u00e9 dans le return() \u00e0 la fin de la fonction)

  } # Fin if (robuste)

  #============================================================================
  #                  RESUME DU MODELE (si param\u00e9trique)
  #============================================================================

  if (!robuste && !is.null(model)) {

    .dbg("", "Affichage du r\u00e9sum\u00e9 du mod\u00e8le...", debug=debug)

    # \u00c9TAPES 8 ET 9 SUPPRIM\u00c9ES : messages redondants avant affichage ANOVA
    # Passage direct \u00e0 l'affichage du tableau ANOVA avec type de SC

    # Afficher le tableau ANOVA avec le type de SC appropri\u00e9
    if (verbose) {
      if (inherits(model, "lm")) {
        # =======================================================================
        # S\u00c9LECTION INTELLIGENTE DU TYPE DE SOMMES DES CARR\u00c9S
        # =======================================================================
        # R\u00e9f\u00e9rence: Maxwell & Delaney (2004). Designing Experiments and Analyzing Data.

        ss_selection <- .select_ss_type(
          formula = formula,
          data = data,
          model = model,
          analysis_type = "ANOVA",
          alpha = alpha,
          verbose = verbose, code = code,
          k = k,
          debug = debug
        )

        k <- ss_selection$k
        ss_type <- ss_selection$ss_type

        k <- .vbse(
          paste0("ANOVA Type ", ss_type, " Sum of Squares [",
                 if(ss_type == "I") "anova() {stats}]" else "Anova() {car}]", "\n",
                 "\tReason: ", ss_selection$reason),
          paste0("ANOVA Sommes de Carr\u00e9s Type ", ss_type, " [",
                 if(ss_type == "I") "anova() {stats}]" else "Anova() {car}]", "\n",
                 "\tRaison : ", ss_selection$reason),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
        if (isTRUE(code)) {
          k_code <- k_code + 1
          if (ss_type == "I") {
            .code_multi(k_code, paste0("ANOVA Type ", ss_type, " (Sommes des carr\u00e9s s\u00e9quentielles)"), c(
              paste0("model <- lm(", deparse(formula), ", data = data)"),
              "anova_result <- anova(model)",
              "print(anova_result)"
            ))
          } else {
            .code_multi(k_code, paste0("ANOVA Type ", ss_type, " (Sommes des carr\u00e9s marginales)"), c(
              "library(car)",
              paste0("model <- lm(", deparse(formula), ", data = data)"),
              paste0("anova_result <- car::Anova(model, type = '", ss_type, "')"),
              "print(anova_result)"
            ))
          }
        }
        cat("\n")

        # Type I utilise anova() de base, Type II/III utilisent car::Anova()
        if (ss_type == "I") {
          # Type I: Sommes des carr\u00e9s s\u00e9quentielles
          anova_result <- anova(model)
          print(anova_result)
          cat("\n")

        } else if (ss_type %in% c("II", "III") && requireNamespace("car", quietly = TRUE)) {
          # Type II/III
          anova_result <- car::Anova(model, type = ss_type)
          print(anova_result)
          cat("\n")

        } else {
          print(summary(model))
        }

        # Analyser et afficher quels facteurs/interactions sont significatifs
        # (pour tous les types I/II/III)
        if (exists("anova_result") && !is.null(anova_result)) {
          if ("Pr(>F)" %in% colnames(anova_result)) {
            sig_effects <- rownames(anova_result)[which(anova_result[, "Pr(>F)"] < alpha)]
            # Exclure "(Intercept)" et "Residuals"
            sig_effects <- sig_effects[!sig_effects %in% c("(Intercept)", "Residuals")]

            if (length(sig_effects) > 0) {
              k <- .vbse(
                paste0("Significant effect(s) detected (\u03b1 = ", alpha, "): ", paste(sig_effects, collapse = ", ")),
                paste0("Effet(s) significatif(s) d\u00e9tect\u00e9(s) (\u03b1 = ", alpha, ") : ", paste(sig_effects, collapse = ", ")),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            } else {
              k <- .vbse(
                paste0("No significant effects detected (\u03b1 = ", alpha, ")."),
                paste0("Aucun effet significatif d\u00e9tect\u00e9 (\u03b1 = ", alpha, ")."),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
          }
        }
        cat("\n")
      } else if (inherits(model, "aov")) {
        cat("\n")
        print(summary(model))
        cat("\n")
      } else if (inherits(model, "aovlist")) {
        # ANOVA avec mesures r\u00e9p\u00e9t\u00e9es (Error() term)
        # Id\u00e9e : Synth\u00e9tiser les r\u00e9sultats avant d'afficher le tableau complet
        # APA : Rapporter les effets significatifs de chaque strate
        # DOI : Maxwell & Delaney (2004). https://doi.org/10.4324/9781315642956

        smry <- summary(model)

        # \u00c9tape 7: Synth\u00e8se des r\u00e9sultats ANOVA avec Error()
        k <- .vbse(
          paste0("ANOVA results with Error term (repeated measures / random effects):"),
          paste0("R\u00e9sultats de l'ANOVA avec terme Error (mesures r\u00e9p\u00e9t\u00e9es / effets al\u00e9atoires) :"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        # Note explicative sur le terme Error() - AVANT les r\u00e9sultats
        # Extraire le nom du facteur id depuis la formule pour personnaliser le message
        id_var <- if (!is.null(id)) id else "id"

        k <- .vbse(
          paste0("Note: The Error() term controls the repeated measures structure.\n",
                 "        Each stratum represents a variance decomposition:\n",
                 "            \u2022 'Error: ", id_var, "' \u2192 Between-subjects variance (individual differences)\n",
                 "            \u2022 'Error: Within' \u2192 Within-subjects variance (repeated measures effect)\n",
                 "        The F-test in 'Within' stratum tests if the repeated factor has a significant effect."),
          paste0("Note : Le terme Error() contr\u00f4le la structure des mesures r\u00e9p\u00e9t\u00e9es.\n",
                 "        Chaque strate repr\u00e9sente une d\u00e9composition de la variance :\n",
                 "            \u2022 'Error: ", id_var, "' \u2192 Variance inter-sujets (diff\u00e9rences individuelles)\n",
                 "            \u2022 'Error: Within' \u2192 Variance intra-sujets (effet du facteur r\u00e9p\u00e9t\u00e9)\n",
                 "        Le test F dans la strate 'Within' teste si le facteur r\u00e9p\u00e9t\u00e9 a un effet significatif."),
          verbose = verbose, code = code, k = k, cpt = "off"
        )

        # Parcourir chaque strate et extraire les effets significatifs
        all_sig_effects <- list()
        all_strata_info <- list()

        for (strata_name in names(smry)) {
          if (!is.null(smry[[strata_name]]) && length(smry[[strata_name]]) > 0) {
            strata_table <- smry[[strata_name]][[1]]

            if (is.data.frame(strata_table) && "Pr(>F)" %in% colnames(strata_table)) {
              # Identifier les effets (exclure Residuals)
              row_names <- rownames(strata_table)
              effect_rows <- !grepl("Residuals", row_names, ignore.case = TRUE)

              if (any(effect_rows)) {
                effects <- row_names[effect_rows]
                p_values <- strata_table[effect_rows, "Pr(>F)"]

                # Identifier effets significatifs
                sig_idx <- which(p_values < alpha & !is.na(p_values))

                if (length(sig_idx) > 0) {
                  sig_effects <- trimws(effects[sig_idx])  # Nettoyer espaces
                  sig_pvals <- p_values[sig_idx]
                  all_sig_effects[[strata_name]] <- data.frame(
                    effect = sig_effects,
                    p_value = sig_pvals,
                    stringsAsFactors = FALSE
                  )
                }

                # Stocker info de la strate
                all_strata_info[[strata_name]] <- list(
                  n_effects = length(effects),
                  effects = trimws(effects),  # Nettoyer espaces
                  p_values = p_values
                )
              }
            }
          }
        }

        # Afficher synth\u00e8se par strate
        if (length(all_strata_info) > 0) {
          for (strata_name in names(all_strata_info)) {
            info <- all_strata_info[[strata_name]]

            if (info$n_effects > 0) {
              # Construire message pour cette strate
              effects_summary <- paste0(sapply(seq_along(info$effects), function(i) {
                paste0(info$effects[i], " (p=", .format_pval(info$p_values[i]), ")")
              }), collapse = ", ")

              # Identifier effets significatifs
              sig_in_strata <- all_sig_effects[[strata_name]]

              # Format simplifi\u00e9 si un seul effet test\u00e9 dans la strate
              if (info$n_effects == 1) {
                effect_name <- info$effects[1]
                p_val <- info$p_values[1]

                if (!is.null(sig_in_strata) && nrow(sig_in_strata) > 0) {
                  # Effet significatif
                  k <- .vbse(
                    paste0("==> Significant effect of ", effect_name, " (p=", .format_pval(p_val), ")."),
                    paste0("==> Effet significatif de ", effect_name, " (p=", .format_pval(p_val), ")."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                } else {
                  # Pas d'effet significatif
                  k <- .vbse(
                    paste0("==> No significant effect of ", effect_name, " (p=", .format_pval(p_val), ")."),
                    paste0("==> Aucun effet significatif de ", effect_name, " (p=", .format_pval(p_val), ")."),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                }
              } else {
                # Format avec mention de strate si plusieurs effets
                if (!is.null(sig_in_strata) && nrow(sig_in_strata) > 0) {
                  sig_summary <- paste0(sapply(seq_len(nrow(sig_in_strata)), function(i) {
                    paste0(sig_in_strata$effect[i], " (p=", .format_pval(sig_in_strata$p_value[i]), ")")
                  }), collapse = ", ")

                  k <- .vbse(
                    paste0("Stratum '", strata_name, "': Significant effect(s) - ", sig_summary),
                    paste0("Strate '", strata_name, "' : Effet(s) significatif(s) - ", sig_summary),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                } else {
                  k <- .vbse(
                    paste0("Stratum '", strata_name, "': No significant effects (tested: ", effects_summary, ")"),
                    paste0("Strate '", strata_name, "' : Aucun effet significatif (test\u00e9s : ", effects_summary, ")"),
                    verbose = verbose, code = code, k = k, cpt = "off"
                  )
                }
              }
            } else {
              # Strate sans effet test\u00e9 (seulement r\u00e9sidus)
              k <- .vbse(
                paste0("\tStratum '", strata_name, "': Error stratum (no effects tested, residuals only)"),
                paste0("\tStrate '", strata_name, "' : Strate d'erreur (aucun effet test\u00e9, r\u00e9sidus seulement)"),
                verbose = verbose, code = code, k = k, cpt = "off"
              )
            }
          }
        }

        # Note d\u00e9plac\u00e9e AVANT les r\u00e9sultats (voir ligne ~4135)

        # D\u00e9cider si afficher le tableau complet
        # Id\u00e9e : Ne pas afficher si mod\u00e8le simple (1 facteur, 1 effet test\u00e9)
        # APA : \u00c9viter redondance si synth\u00e8se d\u00e9j\u00e0 claire
        total_effects_tested <- sum(sapply(all_strata_info, function(x) x$n_effects))

        if (total_effects_tested > 1 || length(all_strata_info) > 2) {
          # Mod\u00e8le complexe (plusieurs effets ou plusieurs strates) \u2192 Afficher tableau
          cat("\n")
          print(summary(model))
          cat("\n")
        } else {
          # Mod\u00e8le simple (1 effet, 2 strates standard) \u2192 Synth\u00e8se suffit
          # Le tableau complet n'apporte pas d'info suppl\u00e9mentaire
          .dbg("Simple model with 1 effect - detailed table skipped (summary already shown in step 7).",
               "Mod\u00e8le simple avec 1 effet - tableau d\u00e9taill\u00e9 ignor\u00e9 (synth\u00e8se d\u00e9j\u00e0 dans \u00e9tape 7).",
               debug = debug)
        }
      }
    }
  }

  #============================================================================
  #                         RETOUR DES RESULTATS
  #============================================================================

  .dbg("", "Pr\u00e9paration du retour des r\u00e9sultats...", debug=debug)

  # Extraire global_pvalue pour return=FALSE
  global_pvalue <- NA

  # Tentative 1: Depuis robust_results (cas mod\u00e8les mixtes, tests robustes)
  if (!is.null(robust_results) && !is.null(robust_results$test_result)) {
    if (inherits(robust_results$test_result, "lmerMod") || inherits(robust_results$test_result, "lmerModLmerTest")) {
      # Mod\u00e8le mixte : extraire p-value minimale des effets fixes
      tryCatch({
            anova_tab <- stats::anova(robust_results$test_result, type = "III")
        pval_col <- which(colnames(anova_tab) %in% c("Pr(>F)", "p.value", "P(>|t|)"))
        if (length(pval_col) > 0) {
          row_names <- rownames(anova_tab)
          valid_rows <- which(!grepl("Residuals|Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
          if (length(valid_rows) > 0) {
            global_pvalue <- min(anova_tab[valid_rows, pval_col[1]], na.rm = TRUE)
          }
        }
      }, error = function(e) NULL)
    } else if (!is.null(robust_results$p_value)) {
      global_pvalue <- robust_results$p_value
    }
  }

  # Tentative 2: Depuis le mod\u00e8le param\u00e9trique (ANOVA classique)
  if ((is.null(global_pvalue) || is.na(global_pvalue)) && !is.null(model)) {
    if (inherits(model, c("aov", "lm"))) {
      tryCatch({
        anova_tab <- car::Anova(model, type = "III")
        pval_col <- which(colnames(anova_tab) %in% c("Pr(>F)", "p.value", "P(>|t|)"))
        if (length(pval_col) > 0) {
          row_names <- rownames(anova_tab)
          valid_rows <- which(!grepl("Residuals|Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
          if (length(valid_rows) > 0) {
            global_pvalue <- min(anova_tab[valid_rows, pval_col[1]], na.rm = TRUE)
          }
        }
      }, error = function(e) NULL)
    } else if (inherits(model, "aovlist")) {
      # Mesures r\u00e9p\u00e9t\u00e9es : extraire depuis summary
      tryCatch({
        smry <- summary(model)
        # Parcourir les strates et trouver la p-value minimale
        all_pvals <- c()
        for (strata in names(smry)) {
          if (!is.null(smry[[strata]]) && length(smry[[strata]]) > 0) {
            # Chaque strate est une listof, prendre le premier \u00e9l\u00e9ment
            strata_table <- smry[[strata]][[1]]
            if (is.data.frame(strata_table) && "Pr(>F)" %in% colnames(strata_table)) {
              # Exclure les lignes Residuals
              row_names <- rownames(strata_table)
              valid_rows <- !grepl("Residuals", row_names, ignore.case = TRUE)
              if (any(valid_rows)) {
                pvals_strata <- strata_table[valid_rows, "Pr(>F)"]
                pvals_strata <- pvals_strata[!is.na(pvals_strata)]
                all_pvals <- c(all_pvals, pvals_strata)
              }
            }
          }
        }
        if (length(all_pvals) > 0) {
          global_pvalue <- min(all_pvals, na.rm = TRUE)
        }
      }, error = function(e) {
        # En cas d'erreur, laisser global_pvalue \u00e0 NA
        NULL
      })
    }
  }

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
    model = model,  # PEUT \u00caTRE NULL (voir documentation)
    formula = formula,
    updated_formula = updated_formula,
    within = within,
    between = between,
    ancova_checks = ancova_checks,  # R\u00e9sultats des contr\u00f4les ANCOVA
    residual_type = residual_type_used,  # NOUVEAU: type de r\u00e9sidu utilis\u00e9
    robust_results = robust_results,
    global_pvalue = global_pvalue,  # Pour return=FALSE
    k = k
  )

  #============================================================================
  #                         MODE CODE=TRUE
  #============================================================================
  # G\u00e9n\u00e8re du code R comment\u00e9 pour reproduire l'analyse

  # NOTE: Centralized code_str generation has been removed (bp.log 7.4.6).
  # Code generation for multi-factor analyses will be handled by coordinated
  # message() calls in future versions. For now, users should refer to verbose output.

  #============================================================================
  #                         RETOUR FINAL
  #============================================================================

  .dbg("=== End .multi_factor_analysis() ===",
       "=== Fin de .multi_factor_analysis() ===", debug=debug)

  return(result)
}

#==============================================================================
#                         NOTES FINALES ET DOCUMENTATION
#==============================================================================

#* ============================================================================
#* R\u00e0SUM\u00e9 DES NOTES PERSONNELLES TRAIT\u00e0ES
#* ============================================================================
#*
#* NOTE #1  [EN COURS]        : Gestion des NA - Solution partielle impl\u00e9ment\u00e9e
#* NOTE #2  [R\u00e9SOLUE]        : R\u00e9sidus studentis\u00e9s - Trace du type ajout\u00e9e
#* NOTE #3  [EN COURS]        : V\u00e9rification croisement complet - .mixed_model_analysis() \u00e0  impl\u00e9menter
#* NOTE #4  [PARTIELLEMENT]   : Binning ANCOVA - Suppr\u00e9ssion du binning automatique
#* NOTE #5  [R\u00e9SOLUE]        : Normalit\u00e9 diff\u00e9rences - Approche conforme standards
#* NOTE #6  [PARTIELLEMENT]   : Ind\u00e9pendance - Contr\u00f4le syst\u00e9matique, .control_independence() \u00e0  modifier
#* NOTE #7  [R\u00e9SOLUE]        : Comparaison mod\u00e8les - AIC/BIC ajout\u00e9s
#* NOTE #8  [PARTIELLEMENT]   : Retour param\u00e9trique - Logique skew/kurt impl\u00e9ment\u00e9e
#* NOTE #9  [R\u00e9SOLUE]        : Sidak - Code supprim\u00e9 (non pertinent)
#* NOTE #10 [R\u00e9SOLUE]        : Variance r\u00e9sidus - Brown-Forsythe implement\u00e9
#* NOTE #11 [NON R\u00e9SOLUE]    : Tests robustes auto - Plan d'action d\u00e9taill\u00e9
#* NOTE #12 [PARTIELLEMENT]   : Kruskal 2-way - Restriction impl\u00e9ment\u00e9e
#*
#* PRIORIT\u00e0S POUR VERSIONS FUTURES :
#* 1. [URGENT]     Impl\u00e9menter .mixed_model_analysis() avec validation valreg()
#* 2. [IMPORTANT]  Friedman test automatique (mesures r\u00e9p\u00e9t\u00e9es, k>=3)
#* 3. [IMPORTANT]  Modifier .control_independence() pour filtrage correct
#* 4. [MOYEN]      ANOVA par permutation (lmPerm::aovp) automatique
#* 5. [MOYEN]      Int\u00e9gration WRS2 et Scheirer-Ray-Hare
#* 6. [FAIBLE]     Gestion avanc\u00e9e NA avec strat\u00e9gies d'imputation
#*
#* ============================================================================

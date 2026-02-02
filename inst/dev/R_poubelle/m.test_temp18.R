################################################################################
#
#                           m.test() - VERSION 18 CORRECTED
#
################################################################################

################################################################################
# FONCTION AUXILIAIRE : .n_levels_g()
################################################################################

#' Calculate the number of levels in grouping variable(s)
#' @param g A factor, vector, or data.frame of factors
#' @return Integer representing the number of levels (or combinations for data.frame)
#' @keywords internal
.n_levels_g <- function(g) {
  if (is.data.frame(g)) {
    # Pour data.frame : produit du nombre de niveaux de chaque facteur
    prod(sapply(g, function(col) {
      if (is.factor(col)) {
        nlevels(droplevels(col))
      } else {
        length(unique(na.omit(col)))
      }
    }))
  } else {
    # Pour vecteur simple
    if (is.factor(g)) {
      nlevels(droplevels(g))
    } else {
      length(unique(na.omit(g)))
    }
  }
}

################################################################################
# FONCTION AUXILIAIRE : .normalize_from_formula()
################################################################################

#' Extract and normalize x, g, and model.frame from formula
#' @param formula Model formula
#' @param data Data frame
#' @param id Subject identifier (optional)
#' @param wt_names Within-subject factor names (optional)
#' @param na_strategy How to handle NA: "preserve" or "omit"
#' @return List with x, g, and mf (model.frame)
#' @keywords internal
.normalize_from_formula <- function(formula, data, id = NULL, wt_names = NULL,
                                    na_strategy = "preserve") {

  # Supprimer Error() pour model.frame
  f_tmp <- .drop_error_term(formula)

  # Créer model.frame
  mf <- model.frame(
    f_tmp,
    data = data,
    na.action = if(na_strategy == "preserve") na.pass else na.omit,
    drop.unused.levels = TRUE
  )

  # Extraire x (réponse)
  x <- model.response(mf)

  # Extraire g (prédicteurs, sans id)
  pred_cols <- colnames(mf)[-1]

  if (!is.null(id) && id %in% pred_cols) {
    pred_cols <- setdiff(pred_cols, id)
  }

  if (length(pred_cols) == 1) {
    g <- mf[[pred_cols[1]]]
    if (!is.factor(g)) g <- factor(g)
  } else if (length(pred_cols) > 1) {
    g <- mf[, pred_cols, drop = FALSE]
  } else {
    g <- NULL
  }

  return(list(x = x, g = g, mf = mf))
}

################################################################################
# FONCTION PRINCIPALE : m.test()
################################################################################

#' Automatic statistical comparison for means, medians, and variances
#'
#' @description
#' This function performs automatic statistical analysis by selecting the most
#' appropriate tests based on data structure and statistical assumptions. It handles:
#' \itemize{
#'   \item One-factor and multi-factor designs (ANOVA, ANCOVA)
#'   \item Univariate and multivariate responses (MANOVA)
#'   \item Paired and unpaired data (repeated measures)
#'   \item Covariates and random effects (mixed models)
#'   \item Parametric and non-parametric tests
#'   \item Post-hoc comparisons and bootstrap validation
#' }
#'
#' @param x Quantitative response. Can be:
#'   \itemize{
#'     \item a numeric vector;
#'     \item a matrix/data frame (multivariate response for MANOVA);
#'     \item a character or integer vector of column names/indices if \code{data} is supplied.
#'   }
#'   Ignored when \code{formula} is used.
#' @param g Grouping factor (or vector). May also be a data frame of factors
#'   (several explanatory variables) when \code{formula} is not used. Ignored when
#'   \code{formula} is used.
#' @param data A data frame that contains the variables referenced by \code{x}, \code{g},
#'   and/or \code{formula}. If omitted, symbols are looked up in the calling
#'   environment. Formulas like \code{data$Y ~ data$G} are also accepted provided
#'   the object \code{data} is visible.
#' @param formula A model formula such as \code{Y ~ G1 + G2} (or
#'   \code{cbind(Y1, Y2, ...) ~ ...} for MANOVA). Notation \code{Y ~ G | id} is internally
#'   converted to \code{Y ~ G + Error(id/...)} for repeated-measures designs.
#'   Covariates can be included as numeric predictors.
#' @param paired Logical. Forces paired/repeated analysis. It is also enabled
#'   automatically when \code{id} is supplied, or in specific repeated-measures patterns.
#' @param id Subject identifier (a column name or index in \code{data}, or a vector
#'   the same length as \code{x}). Supplying \code{id} implies \code{paired = TRUE}
#'   and enables mixed models with random effects.
#' @param wt Within-subject factor(s) for repeated measures. Can be a name,
#'   index, vector, or a data frame of factors present (or injected) in \code{data}.
#' @param within Alias of \code{wt}.
#' @param between Between-subject factor(s). Reserved for explicit specification
#'   in mixed designs. Transmitted to multi-factor analysis functions.
#' @param alpha Global p-value threshold (default \code{0.05}).
#' @param control Name of the category used as control in directed post-hoc
#'   procedures (e.g., Dunnett-type comparisons), when applicable.
#' @param verbose Logical. Print the step-by-step reasoning.
#' @param plot Logical. Draw distribution plots (boxplots, violins, etc.).
#' @param return Logical. If \code{TRUE}, returns a summary object (p-values, compact
#'   letter displays, etc.). If \code{FALSE}, returns raw model output.
#' @param boot Logical. Enable bootstrap for means/medians where relevant.
#' @param iter Number of bootstrap iterations when \code{boot = TRUE}. If \code{0},
#'   a default is chosen as \code{iter <- 1/alpha * 5}.
#' @param conf Confidence level for bootstrap intervals.
#' @param maxcat Maximum number of allowed groups; beyond this, some procedures
#'   may fail.
#' @param silent Logical. Suppress secondary warnings.
#' @param code Logical. Print a simplified R "recipe" that reproduces the analysis
#'   (disables \code{verbose}).
#' @param debug Logical. Emit detailed diagnostic messages (forces \code{code = FALSE}).
#'
#' @return
#' A list containing test results, p-values, group comparisons, and optionally
#' bootstrap confidence intervals. The exact structure depends on the analysis type.
#'
#' @details
#' \strong{NA Handling:}
#' \itemize{
#'   \item For one-factor analyses: Listwise deletion (complete cases only)
#'   \item For multi-factor analyses: NA preserved (may be handled by mixed models)
#' }
#'
#' @importFrom fda.usc fanova.hetero
#' @importFrom agricolae kurtosis skewness SNK.test
#' @importFrom lawstat levene.test
#' @importFrom WRS2 med1way medpb t1way lincon
#' @importFrom FSA dunnTest
#' @importFrom onewaytests bf.test
#' @importFrom vioplot vioplot
#' @importFrom DescTools DunnettTest
#' @import methods
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: One-factor unpaired
#' data1 <- .simul(n = 12, k_H = 3, seed = 123)
#' m.test(A ~ F, data = data1, verbose = TRUE)
#'
#' # Example 2: Repeated measures
#' data2 <- .simul(n = 12, seed = 456)
#' m.test(A ~ G, data = data2, id = "idG", verbose = TRUE)
#'
#' # Example 3: Mixed design
#' data3 <- .simul(n = 20, seed = 789)
#' m.test(A ~ F * Temps + Error(id/Temps), data = data3, verbose = TRUE)
#' }
#'
m.test <- function(x = NULL, g = NULL, data = NULL, formula = NULL,
                   paired = FALSE, id = NULL, wt = NULL, within = NULL, between = NULL,
                   alpha = 0.05, control = NULL, verbose = TRUE, plot = TRUE,
                   return = TRUE, boot = TRUE, iter = 0, conf = 0.95,
                   maxcat = 50, silent = TRUE,
                   code = FALSE, debug = FALSE) {

  ################################################################################
  # BLOC 1 : INITIALISATION
  ################################################################################

  call_expr <- match.call()
  k <- 0

  if (code == TRUE) verbose <- FALSE
  if (debug == TRUE) {
    verbose <- FALSE
    code <- FALSE
  }

  if (iter == 0) iter <- 1/alpha * 5

  ################################################################################
  # BLOC 2 : VALIDATION DES PARAMÈTRES PRINCIPAUX
  ################################################################################

  if (!is.null(data) && !is.data.frame(data)) {
    .exit("'data' must be a data.frame.",
          "'data' doit être un data.frame.",
          verbose = verbose, return = return)
  }

  ################################################################################
  # BLOC 3 : DÉTECTION ET NORMALISATION DE LA FORMULE
  ################################################################################

  check_manova <- FALSE
  check_anova <- FALSE

  if (is.null(formula) && inherits(x, "formula")) {
    .dbg("x appears to be a formula. Transferring to formula.",
         "`x` semble être une formule. Transfert vers `formula`.",
         debug = debug)
    formula <- x
    x <- NULL
  }

  if ("formula" %in% names(call_expr)) {
    formula <- call_expr$formula
  }

  ################################################################################
  # BLOC 4 : ANALYSE DE LA VARIABLE RÉPONSE 'x'
  ################################################################################

  if (!is.null(x)) {

    if (is.data.frame(x) || is.matrix(x)) {

      num_cols <- sapply(x, is.numeric)
      if (!all(num_cols)) {
        .exit("'x' contains non-numeric columns.",
              "'x' présente des colonnes non numériques.",
              verbose = verbose, return = return)
      }

      if (is.data.frame(g) && nrow(g) != nrow(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas les mêmes dimensions.",
              verbose = verbose, return = return)
      }

      if (ncol(x) == 1) {
        x <- as.vector(x[, 1])
      } else {
        check_manova <- TRUE
      }

    } else if (is.vector(x) && is.null(data)) {

      if (!is.numeric(x)) {
        .exit("'x' non-numeric.",
              "'x' non numérique.",
              verbose = verbose, return = return)
      }

      if (is.vector(g) && length(g) != length(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas la même taille.",
              verbose = verbose, return = return)
      }

      if (is.data.frame(g) && nrow(g) != length(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas les mêmes dimensions.",
              verbose = verbose, return = return)
      }

    } else if (!is.null(data)) {

      if (is.numeric(x)) {
        if (!all(x %in% seq_len(ncol(data)))) {
          .exit("x contains invalid indices for data.",
                "x contient des indices invalides pour data.",
                verbose = verbose, return = return)
        }
        x <- data[, x, drop = FALSE]

      } else if (is.character(x)) {
        if (!all(x %in% colnames(data))) {
          .exit("x contains column names that do not exist in data.",
                "x contient des noms de colonnes qui n'existent pas dans data.",
                verbose = verbose, return = return)
        }
        x <- data[, x, drop = FALSE]

      } else {
        .exit("When 'data' is provided, 'x' must be numeric indices or column names.",
              "Quand 'data' est fourni, 'x' doit être des indices ou des noms de colonnes.",
              verbose = verbose, return = return)
      }

      if (ncol(x) == 1) {
        x <- as.vector(x[, 1])
      } else {
        check_manova <- TRUE
      }
    }
  }

  ################################################################################
  # BLOC 5 : ANALYSE DE LA VARIABLE GROUPANTE 'g'
  ################################################################################

  if (!is.null(g) && is.null(formula)) {

    if (is.data.frame(g)) {
      check_anova <- TRUE

      for (i in seq_along(g)) {
        if (!is.factor(g[[i]])) {
          g[[i]] <- factor(g[[i]])
        }
      }

    } else if (!is.null(data)) {

      if (is.numeric(g) && length(g) == 1) {
        if (!g %in% seq_len(ncol(data))) {
          .exit("g contains an invalid index for data.",
                "g contient un indice invalide pour data.",
                verbose = verbose, return = return)
        }
        g <- data[, g]

      } else if (is.character(g) && length(g) == 1) {
        if (!g %in% colnames(data)) {
          .exit("g does not exist in data.",
                "g n'existe pas dans data.",
                verbose = verbose, return = return)
        }
        g <- data[, g]

      } else if ((is.numeric(g) || is.character(g)) && length(g) > 1) {

        if (is.numeric(g)) {
          if (!all(g %in% seq_len(ncol(data)))) {
            .exit("g contains invalid indices for data.",
                  "g contient des indices invalides pour data.",
                  verbose = verbose, return = return)
          }
          g <- data[, g, drop = FALSE]
        } else {
          if (!all(g %in% colnames(data))) {
            .exit("g contains column names that do not exist in data.",
                  "g contient des noms de colonnes qui n'existent pas dans data.",
                  verbose = verbose, return = return)
          }
          g <- data[, g, drop = FALSE]
        }

        check_anova <- TRUE
      }

      if (!is.data.frame(g) && !is.factor(g)) {
        g <- factor(g)
      }
    }

    if (!is.data.frame(g) && !is.factor(g)) {
      g <- factor(g)
    }
  }

  ################################################################################
  # BLOC 6 : GESTION DES ALIAS DE PARAMÈTRES
  ################################################################################

  if (!is.null(within)) {
    if (!is.null(wt)) {
      .exit("Cannot specify both 'wt' and 'within'.",
            "Impossible de spécifier à la fois 'wt' et 'within'.",
            verbose = verbose, return = return)
    }
    wt <- within
  }

  ################################################################################
  # BLOC 7 : NORMALISATION DU PARAMÈTRE 'wt'
  ################################################################################

  wt_names <- NULL

  if (!is.null(wt)) {

    if (is.data.frame(wt)) {
      wt_names <- colnames(wt)
      if (!is.null(data)) {
        data <- cbind(data, wt)
      } else {
        data <- wt
      }

    } else if (!is.null(data)) {

      if (is.numeric(wt)) {
        if (!all(wt %in% seq_len(ncol(data)))) {
          .exit("wt contains invalid indices for data.",
                "wt contient des indices invalides pour data.",
                verbose = verbose, return = return)
        }
        wt_names <- colnames(data)[wt]

      } else if (is.character(wt)) {
        if (!all(wt %in% colnames(data))) {
          .exit("wt contains column names that do not exist in data.",
                "wt contient des noms de colonnes qui n'existent pas dans data.",
                verbose = verbose, return = return)
        }
        wt_names <- wt
      }

    } else {
      wt_temp_name <- paste0("wt_", sample(1000:9999, 1))
      if (is.null(data)) {
        data <- data.frame(temp_wt = wt)
      } else {
        data[[wt_temp_name]] <- wt
      }
      wt_names <- wt_temp_name
    }

    if (!is.null(wt_names)) paired <- TRUE
  }

  ################################################################################
  # BLOC 8 : NORMALISATION DU PARAMÈTRE 'id'
  ################################################################################

  pair_vector <- NULL

  if (!is.null(id)) {

    paired <- TRUE

    if (is.null(data)) {
      pair_temp_name <- paste0("id_", sample(1000:9999, 1))
      data <- data.frame(temp_id = id)
      colnames(data)[1] <- pair_temp_name
      id <- pair_temp_name

    } else {

      if (is.numeric(id) && length(id) == 1) {
        if (!id %in% seq_len(ncol(data))) {
          .exit("id contains an invalid index for data.",
                "id contient un indice invalide pour data.",
                verbose = verbose, return = return)
        }
        id <- colnames(data)[id]

      } else if (is.character(id) && length(id) == 1) {
        if (!id %in% colnames(data)) {
          .exit("id does not exist in data.",
                "id n'existe pas dans data.",
                verbose = verbose, return = return)
        }

      } else if (length(id) > 1) {
        if (is.null(data)) {
          data <- data.frame(pair_id = id)
        } else if (nrow(data) != length(id)) {
          .exit("Length of 'id' does not match data rows.",
                "La longueur de 'id' ne correspond pas au nombre de lignes de data.",
                verbose = verbose, return = return)
        } else {
          pair_temp_name <- paste0("id_", sample(1000:9999, 1))
          data[[pair_temp_name]] <- id
          id <- pair_temp_name
        }
      }
    }
  }

  ################################################################################
  # BLOC 9 : TRAITEMENT DU CAS FORMULE
  ################################################################################

  if (!is.null(formula)) {

    if (!inherits(formula, "formula")) {
      .exit("The 'formula' parameter must be a valid formula.",
            "Le paramètre 'formula' doit être une formule valide.",
            verbose = verbose, return = return)
    }

    # --- 9.1 : Détection MANOVA ---
    lhs <- deparse(formula[[2]])
    if (grepl("cbind\\(", lhs)) {
      check_manova <- TRUE
      .dbg("Multivariate response detected (cbind).",
           "Réponse multivariée détectée (cbind).",
           debug = debug)
    }

    # --- 9.2 : Détection multi-facteurs ---
    rhs <- attr(terms(formula), "term.labels")
    rhs_clean <- rhs[!grepl("^Error\\(", rhs)]

    if (length(rhs_clean) > 1 || any(grepl("[*:]", rhs_clean))) {
      check_anova <- TRUE
      .dbg("Multi-factor design detected.",
           "Design multi-facteurs détecté.",
           debug = debug)
    }

    # --- 9.3 : Gestion notation | ---
    if (grepl("\\|", deparse(formula))) {
      split_formula <- strsplit(deparse(formula), "\\|")[[1]]

      if (length(split_formula) == 2) {
        left_side <- trimws(split_formula[1])
        right_side <- trimws(split_formula[2])

        formula <- as.formula(paste(left_side, "+ Error(", right_side, ")"))

        .dbg(paste0("Formula adjusted with 'Error()': ", deparse(formula)),
             paste0("Formule ajustée avec 'Error()' : ", deparse(formula)),
             debug = debug)
      } else {
        .exit("The formula contains a malformed '|'.",
              "La formule contient un '|' mal formé.",
              return = return, verbose = verbose)
      }
    }

    # --- 9.4 : Environnement ---
    if (!is.null(data) && is.data.frame(data)) {
      env <- new.env(parent = parent.frame())
      list2env(data, env)

      vars <- all.vars(formula)
      for (variable in vars) {
        if (!exists(variable, envir = env, inherits = FALSE)) {
          assign(variable, get(variable, envir = parent.frame(), inherits = TRUE),
                 envir = env)
        }
      }
      environment(formula) <- env

    } else if (!is.null(data) && !is.data.frame(data)) {
      .exit("The 'data' parameter must be a valid data frame.",
            "Le paramètre 'data' doit être un data.frame valide.",
            verbose = verbose, return = return)
    }

    # --- 9.5 : EXTRACTION x et g (AVANT modification Error) ---
    # Extraire avec formule SANS modifications id/Error
    formula_for_extraction <- formula

    # Supprimer Error() s'il existe pour l'extraction
    f_tmp <- if (!is.null(data)) .strip_data_dollar_safe(formula_for_extraction, data) else formula_for_extraction
    f_mf  <- .drop_error_term(f_tmp)

    mf <- tryCatch(
      model.frame(
        f_mf,
        data = if (!is.null(data)) data else parent.frame(),
        na.action = na.pass,
        drop.unused.levels = TRUE
      ),
      error = function(e) {
        .exit("Error evaluating the formula. Check your data and formula.",
              "Erreur lors de l'évaluation de la formule. Vérifiez vos données et votre formule.",
              return = return, verbose = verbose)
      }
    )

    # NE PAS TRIER - garder l'ordre original de data comme v17
    # Le tri peut perturber la détection des patterns dans .multi_factor_analysis()

    # Extraire x
    x <- model.response(mf)

    if (is.matrix(x) && ncol(x) > 1) {
      check_manova <- TRUE
    }

    # Extraire g (SANS id)
    pred_cols <- colnames(mf)[-1]

    if (!is.null(id) && id %in% pred_cols) {
      pred_cols <- setdiff(pred_cols, id)
    }

    if (length(pred_cols) == 1) {
      g <- mf[[pred_cols[1]]]
      if (!is.factor(g)) g <- factor(g)
    } else if (length(pred_cols) > 1) {
      g <- mf[, pred_cols, drop = FALSE]
      check_anova <- TRUE
    } else {
      .exit("No predictors found in formula after excluding id.",
            "Aucun prédicteur trouvé dans la formule après exclusion de id.",
            verbose = verbose, return = return)
    }

    # --- 9.6 : CONSTRUCTION FORMULE POUR ANALYSE (AVEC Error si nécessaire) ---
    formula_for_analysis <- formula
    has_error_already <- grepl("Error\\(", deparse(formula))

    if (!is.null(id) && !has_error_already) {

      if (!is.null(wt_names)) {
        # Cas avec wt
        for (wt_var in wt_names) {
          if (!wt_var %in% colnames(data)) {
            .exit(paste0("The specified variable `", wt_var, "` does not exist in data."),
                  paste0("La variable spécifiée `", wt_var, "` n'existe pas dans `data`."),
                  verbose = verbose, return = return)
          }
        }

        for (wt_var in wt_names) {
          if (!wt_var %in% all.vars(formula_for_analysis)) {
            formula_for_analysis <- update(formula_for_analysis, paste0(". ~ . + ", wt_var))
          }
        }

        if (length(wt_names) > 1) {
          interaction_terms <- paste(wt_names, collapse = " * ")
          formula_for_analysis <- update(formula_for_analysis,
                                         paste0(". ~ . + Error(", id, "/(", interaction_terms, "))"))
        } else {
          formula_for_analysis <- update(formula_for_analysis,
                                         paste0(". ~ . + Error(", id, "/", wt_names[1], ")"))
        }

      } else {
        # Cas sans wt : utiliser la logique v17
        # Extraire les prédicteurs de la formule ORIGINALE
        pred_vars <- all.vars(formula[[3]])

        if (length(pred_vars) == 1) {
          # Un seul prédicteur : Error(id/predictor)
          # Ex: A ~ G avec id="idG" → A ~ G + Error(idG/G)
          formula_for_analysis <- update(formula_for_analysis,
                                         paste0(". ~ . + Error(", id, "/", pred_vars[1], ")"))
        } else {
          # Plusieurs prédicteurs : Error(id) seulement
          formula_for_analysis <- update(formula_for_analysis,
                                         paste0(". ~ . + Error(", id, ")"))
        }
      }

      .dbg(paste0("Formula for analysis: ", deparse(formula_for_analysis)),
           paste0("Formule pour l'analyse : ", deparse(formula_for_analysis)),
           debug = debug)
    }

    # Mettre à jour formula vers formula_for_analysis
    formula <- formula_for_analysis
  }

  ################################################################################
  # BLOC 10 : CONTRÔLES FINAUX
  ################################################################################

  if (is.null(x) || (is.null(g) && is.null(formula))) {
    .exit("Both 'x' and 'g' (or 'formula') must be specified.",
          "'x' et 'g' (ou 'formula') doivent être spécifiés.",
          verbose = verbose, return = return)
  }

  lev_g <- .n_levels_g(g)

  if (lev_g > maxcat) {
    .exit(paste0("Too many groups (", lev_g, " > ", maxcat, ")."),
          paste0("Trop de groupes (", lev_g, " > ", maxcat, ")."),
          verbose = verbose, return = return)
  }

  if (lev_g < 2) {
    .exit("At least 2 groups are required for comparison.",
          "Au moins 2 groupes sont nécessaires pour une comparaison.",
          verbose = verbose, return = return)
  }

  # Note: Si des doublons (id, group) existent, ils seront gérés par les sous-fonctions
  # d'analyse qui peuvent agréger ou prendre la première observation automatiquement.

  ################################################################################
  # BLOC 11 : GESTION DES NA
  ################################################################################

  if (!check_anova && !check_manova) {
    if (is.vector(x) && is.vector(g)) {
      complete_cases <- complete.cases(x, g)
      n_na <- sum(!complete_cases)

      if (n_na > 0) {
        if (verbose && !silent) {
          .vbse(
            paste0("Warning: ", n_na, " observations with missing values removed (listwise deletion)."),
            paste0("Attention : ", n_na, " observations avec valeurs manquantes supprimées (suppression listwise)."),
            verbose = verbose, k = k, cpt = "off"
          ) -> k
        }

        x <- x[complete_cases]
        g <- g[complete_cases]

        if (!is.null(data)) {
          data <- data[complete_cases, , drop = FALSE]
        }
      }
    }
  }

  if (check_anova || check_manova) {
    .dbg("NA values preserved for multi-factor/MANOVA analysis.",
         "Valeurs NA conservées pour analyse multi-facteurs/MANOVA.",
         debug = debug)
  }

  ################################################################################
  # BLOC 12 : VÉRIFICATION DES CROISEMENTS
  ################################################################################

  if (check_anova == TRUE && check_manova == FALSE) {

    if (debug) {
      cat("\n")
      cat("===========================\n")
      cat("BLOC 12 : Vérification des croisements\n")
      cat("===========================\n")
    }

    data_model <- model.frame(.formulator_safe(formula, mode = "linear", debug = debug),
                              data = data)

    f_lin2 <- .formulator_safe(formula, mode = "linear", debug = debug)
    f_lin2_noErr <- .drop_error_term(f_lin2)
    variables_explicatives <- attr(terms(f_lin2_noErr), "term.labels")

    non_numeriques <- variables_explicatives[
      !sapply(variables_explicatives, function(var) is.numeric(data_model[[var]]))
    ]

    if (length(non_numeriques) > 0) {
      variables_categoriques <- non_numeriques
      combinaisons <- table(data_model[variables_categoriques])

      if (min(combinaisons) == 0) {
        .dbg("Some crossed categories are not specified. Check your data.",
             "Certaines catégories croisées ne sont pas renseignées. Vérifiez vos données.",
             debug = debug)
      } else {
        .dbg("The crossed categories are correctly specified.",
             "Les catégories croisées sont bien renseignées.",
             debug = debug)
      }
    } else {
      .dbg("All explanatory variables are numeric.",
           "Toutes les variables explicatives sont numériques.",
           debug = debug)
    }
  }

  ################################################################################
  # BLOC 13 : DÉTECTION TYPE D'ANALYSE
  ################################################################################

  lev_g <- .n_levels_g(g)
  has_error <- (!is.null(formula) && grepl("Error\\(", deparse(formula)))

  if (isTRUE(paired) && (has_error || (!is.null(id) && lev_g > 2))) {
    check_anova <- TRUE
  }

  if (debug) {
    cat("\n")
    cat("===========================\n")
    cat("BLOC 13 : Flags de routage\n")
    cat("===========================\n")
    cat("check_manova  = ", check_manova, "\n")
    cat("check_anova   = ", check_anova, "\n")
    cat("paired        = ", paired, "\n")
    cat("lev_g         = ", lev_g, "\n")
    cat("has_error     = ", has_error, "\n")
    cat("\n")
  }

  ################################################################################
  # BLOC 14 : ROUTAGE
  ################################################################################

  .dbg(NULL, "========================", debug = debug)
  .dbg(NULL, "   ROUTAGE ANALYSE", debug = debug)
  .dbg(NULL, "========================", debug = debug)

  route <- NULL

  if (isTRUE(paired)) {

    if (isTRUE(check_manova) || isTRUE(check_anova)) {
      route <- if (isTRUE(check_manova)) "MANOVA(repeated)" else "ANOVA multi-facteurs (repeated)"

      .vbse(paste0("Routing: ", route),
            paste0("Routage : ", route),
            verbose = verbose, k = k, cpt = "off") -> k

      if (isTRUE(check_manova)) {
        bilan <- .manova_analysis(x = x, g = g, formula = formula, data = data,
                                  alpha = alpha, paired = TRUE, id = id, between = between,
                                  k = k, code = code, debug = debug, verbose = verbose)
      } else {
        bilan <- .multi_factor_analysis(x = x, g = g, formula = formula, data = data,
                                        alpha = alpha, paired = TRUE, id = id, between = between,
                                        k = k, code = code, debug = debug, verbose = verbose)
      }

    } else {
      lev <- nlevels(droplevels(factor(g)))

      if (lev > 2 && is.null(id)) {
        .exit("For paired designs with >2 levels, 'id' is required.",
              "Pour un plan apparié avec >2 modalités, 'id' est requis.",
              verbose = verbose, return = return)
      }

      route <- if (lev == 2) "Paired (2 levels)" else "Paired (k>2 levels)"

      .vbse(paste0("Routing: ", route),
            paste0("Routage : ", route),
            verbose = verbose, k = k, cpt = "off") -> k

      bilan <- .one_factor_analysis(x = x, g = g, formula = formula, data = data,
                                    alpha = alpha, paired = TRUE, id = id,
                                    k = k, code = code, debug = debug, verbose = verbose)
    }

  } else {

    if (isTRUE(check_manova)) {
      .vbse("Routing: MANOVA",
            "Routage : MANOVA",
            verbose = verbose, k = k, cpt = "off") -> k

      bilan <- .manova_analysis(x = x, g = g, formula = formula, data = data,
                                alpha = alpha, paired = FALSE, id = id, between = between,
                                k = k, code = code, debug = debug, verbose = verbose)

    } else if (isTRUE(check_anova)) {
      .vbse("Routing: multi-factor ANOVA",
            "Routage : ANOVA multi-facteurs",
            verbose = verbose, k = k, cpt = "off") -> k

      bilan <- .multi_factor_analysis(x = x, g = g, formula = formula, data = data,
                                      alpha = alpha, paired = FALSE, id = id, between = between,
                                      k = k, code = code, debug = debug, verbose = verbose)

    } else {
      .vbse("Routing: one-factor (unpaired)",
            "Routage : 1 facteur (non apparié)",
            verbose = verbose, k = k, cpt = "off") -> k

      bilan <- .one_factor_analysis(x = x, g = g, formula = formula, data = data,
                                    alpha = alpha, paired = FALSE, id = id,
                                    k = k, code = code, debug = debug, verbose = verbose)
    }
  }

  ################################################################################
  # BLOC 15 : RÉCUPÉRATION RÉSULTATS
  ################################################################################

  x <- bilan[[1]]
  g <- bilan[[2]]
  check_normality <- bilan[[3]]
  check_variance_equal <- bilan[[4]]
  k <- bilan$k

  ################################################################################
  # BLOC 16 : GRAPHIQUES
  ################################################################################

  if ((plot == TRUE) && (length(unique(g)) < maxcat)) {
    .dbg(NULL, "Représentations graphiques.", debug = debug)

    boxplot(x ~ g, col = "cyan")
    vioplot(x ~ g, col = "#00AA0077", add = TRUE)
    stripchart(x ~ g, col = "#FF000088", pch = 16, vertical = TRUE, add = TRUE,
               method = "jitter", jitter = 1/(length(unique(g)) + 2))
  }

  ################################################################################
  # BLOC 17 : POST-HOC
  ################################################################################

  .dbg(NULL, "Tests posts-hocs.", debug = debug)

  if (return == TRUE) {
    synth <- .posthoc(x, g, alpha = alpha, normal = check_normality,
                      var.equal = check_variance_equal,
                      control = control, code = code, debug = debug,
                      verbose = verbose, paired = paired,
                      boot = boot, iter = iter, conf = conf, k = k)
    return(synth)
  } else {
    return(bilan)
  }
}

################################################################################
#                            FIN DE m.test() v18
################################################################################

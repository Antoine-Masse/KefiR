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
#' @note For data.frames, only categorical/factor variables are counted.
#'       Numeric variables (covariates) are excluded to avoid false positives
#'       in "too many groups" checks when ANCOVA models are used.
.n_levels_g <- function(g) {
  if (is.data.frame(g)) {
    # Pour data.frame : produit du nombre de niveaux UNIQUEMENT pour les variables catégorielles
    # Les variables numériques (covariables ANCOVA) sont EXCLUES du comptage
    categorical_cols <- sapply(g, function(col) {
      is.factor(col) || is.character(col)
    })

    if (!any(categorical_cols)) {
      # Aucune variable catégorielle : retourner 1 (pas de groupes)
      return(1)
    }

    # Calculer le produit uniquement pour les colonnes catégorielles
    prod(sapply(g[, categorical_cols, drop = FALSE], function(col) {
      if (is.factor(col)) {
        nlevels(droplevels(col))
      } else {
        length(unique(na.omit(col)))
      }
    }))
  } else {
    # Pour vecteur simple : comportement inchangé
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
#' @param boot_type Character. Force the bootstrap type: \code{"mean"}, \code{"median"}, or \code{"meanbp"}. 
#'   If \code{NULL} (default), automatically determined based on normality tests.
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
#' @importFrom WRS2 med1way medpb2 t1way lincon
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
                   return = TRUE, boot = TRUE, boot_type = NULL, iter = 0, conf = 0.95,
                   maxcat = 50, silent = TRUE,
                   code = FALSE, debug = FALSE) {

  ################################################################################
  # BLOC 1 : INITIALISATION
  ################################################################################

  call_expr <- match.call()
  k <- 0
  pvals <- NULL

  # NOTE PERSO: Ajustement des paramètres par défaut selon return
  # Si return=FALSE et verbose/plot non spécifiés → mettre à FALSE par défaut
  # (cf. Cahier des charges priorité 2)
  if (return == FALSE) {
    if (!("verbose" %in% names(call_expr))) {
      verbose <- FALSE
    }
    if (!("plot" %in% names(call_expr))) {
      plot <- FALSE
    }
  }

  if (code == TRUE) verbose <- FALSE
  if (debug == TRUE) {
    verbose <- FALSE
    code <- FALSE
  }

  if (iter == 0) iter <- 1/alpha * 5

  ################################################################################
  # BLOC 2 : VALIDATION DES PARAMÈTRES PRINCIPAUX
  ################################################################################

  # NOTE: Validation du data.frame reportée après détection du format d'appel
  # pour permettre les appels sans data (format x,g)

  ################################################################################
  # BLOC 3 : DÉTECTION ET NORMALISATION DE LA FORMULE
  ################################################################################

  check_manova <- FALSE
  check_anova <- FALSE

  # NOTE PERSO: Robustesse pour syntaxes non-standard (Cahier des charges priorité 2)
  # Cas 1: m.test(A ~ F, data) → formule capturée dans x, data dans g
  # Cas 2: m.test(A ~ data$F, data) → formule capturée dans x
  # Cas 3: m.test(formula, data) → géré automatiquement

  if (is.null(formula) && inherits(x, "formula")) {
    .dbg("x appears to be a formula. Transferring to formula.",
         "`x` semble être une formule. Transfert vers `formula`.",
         debug = debug)
    formula <- x
    x <- NULL

    # Si g contient un data.frame, c'est probablement l'argument data mal positionné
    if (!is.null(g) && is.data.frame(g) && is.null(data)) {
      .dbg("Detected data argument in 'g' position. Transferring to 'data'.",
           "Argument 'data' détecté en position 'g'. Transfert vers 'data'.",
           debug = debug)
      data <- g
      g <- NULL
    }
  }

  # Cas supplémentaire : formule a été passée via argument nommé
  if ("formula" %in% names(call_expr)) {
    formula <- call_expr$formula
  }

  # Cas formule avec data$ : essayer d'extraire le data.frame depuis l'environnement parent
  if (!is.null(formula) && is.null(data)) {
    .dbg("Formula without data argument detected. Attempting to extract data from formula environment.",
         "Formule sans argument data détectée. Tentative d'extraction depuis l'environnement de la formule.",
         debug = debug)

    # Vérifier si la formule contient des références du type data$var
    formula_str <- deparse(formula)
    if (grepl("\\$", formula_str)) {
      # Extraire le nom du data.frame (avant le $)
      data_name <- gsub("^.*?([a-zA-Z_][a-zA-Z0-9_.]*)\\$.*$", "\\1", formula_str)

      # Essayer de récupérer le data.frame depuis l'environnement parent
      parent_env <- parent.frame()
      if (exists(data_name, envir = parent_env)) {
        data_candidate <- get(data_name, envir = parent_env)
        if (is.data.frame(data_candidate)) {
          .dbg(paste0("Found data.frame '", data_name, "' in parent environment. Using it as data."),
               paste0("Data.frame '", data_name, "' trouvé dans environnement parent. Utilisation comme data."),
               debug = debug)
          data <- data_candidate

          # Nettoyer la formule pour enlever les data$
          formula_clean <- as.formula(gsub(paste0(data_name, "\\$"), "", formula_str))
          formula <- formula_clean
        }
      }
    }
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
              verbose = verbose, code = code, return = return)
      }

      if (is.data.frame(g) && nrow(g) != nrow(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas les mêmes dimensions.",
              verbose = verbose, code = code, return = return)
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
              verbose = verbose, code = code, return = return)
      }

      if (is.vector(g) && length(g) != length(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas la même taille.",
              verbose = verbose, code = code, return = return)
      }

      if (is.data.frame(g) && nrow(g) != length(x)) {
        .exit("'x' and 'g' length differ.",
              "'x' et 'g' n'ont pas les mêmes dimensions.",
              verbose = verbose, code = code, return = return)
      }

    } else if (!is.null(data)) {

      if (is.numeric(x)) {
        if (!all(x %in% seq_len(ncol(data)))) {
          .exit("x contains invalid indices for data.",
                "x contient des indices invalides pour data.",
                verbose = verbose, code = code, return = return)
        }
        x <- data[, x, drop = FALSE]

      } else if (is.character(x)) {
        if (!all(x %in% colnames(data))) {
          .exit("x contains column names that do not exist in data.",
                "x contient des noms de colonnes qui n'existent pas dans data.",
                verbose = verbose, code = code, return = return)
        }
        x <- data[, x, drop = FALSE]

      } else {
        .exit("When 'data' is provided, 'x' must be numeric indices or column names.",
              "Quand 'data' est fourni, 'x' doit être des indices ou des noms de colonnes.",
              verbose = verbose, code = code, return = return)
      }

      if (ncol(x) == 1) {
        x <- as.vector(x[, 1])
      } else {
        check_manova <- TRUE
      }
    }
  }

  # Validation différée du paramètre data (après détection du format d'appel)
  if (!is.null(data) && !is.data.frame(data)) {
    .exit("'data' must be a data.frame.",
          "'data' doit être un data.frame.",
          verbose = verbose, code = code, return = return)
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
                verbose = verbose, code = code, return = return)
        }
        g <- data[, g]

      } else if (is.character(g) && length(g) == 1) {
        if (!g %in% colnames(data)) {
          .exit("g does not exist in data.",
                "g n'existe pas dans data.",
                verbose = verbose, code = code, return = return)
        }
        g <- data[, g]

      } else if ((is.numeric(g) || is.character(g)) && length(g) > 1) {

        if (is.numeric(g)) {
          if (!all(g %in% seq_len(ncol(data)))) {
            .exit("g contains invalid indices for data.",
                  "g contient des indices invalides pour data.",
                  verbose = verbose, code = code, return = return)
          }
          g <- data[, g, drop = FALSE]
        } else {
          if (!all(g %in% colnames(data))) {
            .exit("g contains column names that do not exist in data.",
                  "g contient des noms de colonnes qui n'existent pas dans data.",
                  verbose = verbose, code = code, return = return)
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
            verbose = verbose, code = code, return = return)
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
                verbose = verbose, code = code, return = return)
        }
        wt_names <- colnames(data)[wt]

      } else if (is.character(wt)) {
        if (!all(wt %in% colnames(data))) {
          .exit("wt contains column names that do not exist in data.",
                "wt contient des noms de colonnes qui n'existent pas dans data.",
                verbose = verbose, code = code, return = return)
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
                verbose = verbose, code = code, return = return)
        }
        id <- colnames(data)[id]

      } else if (is.character(id) && length(id) == 1) {
        if (!id %in% colnames(data)) {
          .exit("id does not exist in data.",
                "id n'existe pas dans data.",
                verbose = verbose, code = code, return = return)
        }

      } else if (length(id) > 1) {
        if (is.null(data)) {
          data <- data.frame(pair_id = id)
        } else if (nrow(data) != length(id)) {
          .exit("Length of 'id' does not match data rows.",
                "La longueur de 'id' ne correspond pas au nombre de lignes de data.",
                verbose = verbose, code = code, return = return)
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
            verbose = verbose, code = code, return = return)
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
    # DISTINGUER 2 syntaxes avec | :
    #   1) Syntaxe lmer : Y ~ F + (1|Subject) ou Y ~ F + (var|Subject) → PRÉSERVER
    #   2) Syntaxe simple : Y ~ F | Bloc → TRANSFORMER en Y ~ F + Error(Bloc)
    formula_str_9_3 <- deparse(formula)
    is_lmer_syntax <- grepl("[(][^)]+[|][^)]+[)]", formula_str_9_3)  # (1|id), (var|id)

    if (grepl("\\|", formula_str_9_3) && !is_lmer_syntax) {
      # Syntaxe simple Y ~ F | Bloc : transformer en Error()
      split_formula <- strsplit(formula_str_9_3, "\\|")[[1]]

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
    } else if (is_lmer_syntax) {
      # Syntaxe lmer native (1|Subject) : extraire Subject comme id et simplifier la formule
      # Cette syntaxe indique des mesures répétées/appariées, pas nécessairement un modèle mixte
      # Équivalent à : m.test(Y ~ Group, id = "Subject")

      # Extraire l'id depuis (1|id) ou (var|id)
      lmer_match <- regmatches(formula_str_9_3, regexpr("[|][^)]+[)]", formula_str_9_3))
      if (length(lmer_match) > 0) {
        extracted_id <- trimws(gsub("[|)]", "", lmer_match[1]))

        if (!is.null(data) && extracted_id %in% names(data)) {
          # Assigner id si non déjà défini
          if (is.null(id)) {
            id <- extracted_id
            k <- .vbse(
              paste0("lmer syntax (1|", id, ") detected: treating as repeated measures with id='", id, "'"),
              paste0("Syntaxe lmer (1|", id, ") détectée : traitement comme mesures répétées avec id='", id, "'"),
              verbose = verbose, code = code, k = k, cpt = "on"
            )
          }

          # Simplifier la formule en retirant le terme (1|id)
          formula_clean <- gsub("\\s*\\+\\s*[(][^)]+[|][^)]+[)]", "", formula_str_9_3)
          formula_clean <- gsub("[(][^)]+[|][^)]+[)]\\s*\\+\\s*", "", formula_clean)
          formula_clean <- gsub("[(][^)]+[|][^)]+[)]", "", formula_clean)
          formula <- as.formula(formula_clean)

          .dbg(paste0("Formula simplified from '", formula_str_9_3, "' to '", deparse(formula), "'"),
               paste0("Formule simplifiée de '", formula_str_9_3, "' à '", deparse(formula), "'"),
               debug = debug)
        } else {
          k <- .vbse(
            paste0("Warning: Variable '", extracted_id, "' from lmer syntax not found in data. Keeping formula as-is."),
            paste0("Attention : Variable '", extracted_id, "' de la syntaxe lmer non trouvée dans data. Conservation de la formule."),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
        }
      }
    }

    # --- 9.4 : EXTRACTION id depuis Error() si présent dans formule ---
    # Idée : Si formule contient Error(var), extraire var comme id
    # APA : Maxwell et al. (2018). Designing Experiments, 3rd ed. Chapter 11.
    # DOI : https://doi.org/10.4324/9781315642956
    # Le terme Error() spécifie la variable de groupement (sujet) pour plans appariés

    has_error_already <- grepl("Error\\(", deparse(formula))

    if (has_error_already && is.null(id)) {
      # Extraire la variable dans Error()
      error_match <- regmatches(deparse(formula), regexpr("Error\\(([^/)]+)", deparse(formula)))
      if (length(error_match) > 0) {
        id_from_error <- gsub("Error\\(", "", error_match)
        id_from_error <- trimws(id_from_error)

        .dbg(paste0("Extracted id='", id_from_error, "' from Error() term in formula."),
             paste0("Extraction de id='", id_from_error, "' depuis le terme Error() dans la formule."),
             debug = debug)

        # --- DISTINCTION Error(ID_sujet) vs Error(facteur) ---
        # Idée : Error(var) peut signifier 2 choses:
        #   1) var = ID sujet (beaucoup de modalités, peu d'obs par modalité) => mesures répétées
        #   2) var = facteur catégoriel (peu de modalités, beaucoup d'obs) => strate d'erreur aov()
        # APA : Maxwell & Delaney (2004). Designing Experiments, 2nd ed. Ch 11-12.
        # DOI : https://doi.org/10.4324/9781315642956

        is_subject_id <- TRUE  # Par défaut, supposer ID sujet

        if (!is.null(data) && id_from_error %in% colnames(data)) {
          n_levels <- nlevels(as.factor(data[[id_from_error]]))
          n_obs <- nrow(data)
          ratio <- n_levels / n_obs

          # Heuristique : ratio élevé (>5%) => ID sujet, ratio faible (<5%) => facteur
          if (ratio < 0.05) {
            # Cas facteur catégoriel : GARDER Error() dans formule, NE PAS extraire comme id
            is_subject_id <- FALSE

            .dbg(paste0("Error(", id_from_error, ") has ", n_levels, " levels (",
                       round(ratio*100, 2), "% of observations).\n",
                       "\tThis appears to be a CATEGORICAL FACTOR, not a subject ID.\n",
                       "\t==> Keeping Error() term in formula for standard RM-ANOVA with aov().\n",
                       "\t==> NOT extracting as id parameter (paired will remain FALSE unless explicitly set)."),
                 paste0("Error(", id_from_error, ") a ", n_levels, " niveaux (",
                       round(ratio*100, 2), "% des observations).\n",
                       "\tCeci semble être un FACTEUR CATÉGORIEL, pas un ID sujet.\n",
                       "\t==> Conservation du terme Error() dans la formule pour RM-ANOVA standard avec aov().\n",
                       "\t==> PAS d'extraction comme paramètre id (paired restera FALSE sauf si défini explicitement)."),
                 debug = debug)
          } else {
            # Cas ID sujet : extraire comme id
            .dbg(paste0("Error(", id_from_error, ") has ", n_levels, " levels (",
                       round(ratio*100, 2), "% of observations).\n",
                       "\tThis appears to be a SUBJECT ID variable.\n",
                       "\t==> Extracting as id parameter for paired analysis."),
                 paste0("Error(", id_from_error, ") a ", n_levels, " niveaux (",
                       round(ratio*100, 2), "% des observations).\n",
                       "\tCeci semble être une variable ID SUJET.\n",
                       "\t==> Extraction comme paramètre id pour analyse appariée."),
                 debug = debug)
          }
        }

        # Assigner à id SEULEMENT si c'est un vrai ID sujet
        if (is_subject_id) {
          id <- id_from_error
          paired <- TRUE
        } else {
          # Facteur catégoriel : ne PAS extraire, laisser dans Error()
          id_from_error <- NULL
          # Ne PAS définir paired=TRUE ici
          # (paired peut être TRUE si défini explicitement par l'utilisateur)

          # IMPORTANT: Marquer pour router vers RM-ANOVA avec aov()
          # Idée : Error(facteur) doit être traité par aov(), pas par t-test/Wilcoxon
          # APA : Maxwell & Delaney (2004). Designing Experiments, 2nd ed. Ch 11.
          # DOI : https://doi.org/10.4324/9781315642956
          # Le terme Error(facteur) spécifie strates d'erreur = RM-ANOVA obligatoire
          has_error_term_factor <- TRUE
        }
      }
    }

    # --- 9.4b : VALIDATION Error() TERM (AVANT model.frame) ---
    # Idée : Error() doit porter sur variable de groupement (facteur ou ID discret)
    #        REJETER les variables continues (beaucoup de valeurs uniques)
    # APA : Maxwell, Delaney & Kelley (2018). Designing Experiments and Analyzing Data, 3rd ed.
    #       Chapitre 11-12 : "Error term must specify grouping factor (subjects, blocks)"
    # DOI : https://doi.org/10.4324/9781315642956
    # Référence : Barr et al. (2013). J Mem Lang 68(3):255-278. DOI: 10.1016/j.jml.2012.11.001

    # Fonction helper pour vérifier si variable est groupement valide
    .is_valid_grouping_var <- function(var_data, var_name) {
      # Si déjà facteur ou character => OK
      if (is.factor(var_data) || is.character(var_data)) {
        return(list(valid = TRUE, convert = FALSE))
      }

      # Si numérique, vérifier si c'est un ID discret ou une variable continue
      if (is.numeric(var_data) || is.integer(var_data)) {
        n_total <- length(var_data)
        n_unique <- length(unique(var_data))
        pct_unique <- n_unique / n_total

        # Critère : si >50% valeurs uniques OU >50 valeurs uniques => variable continue
        if (pct_unique > 0.5 || n_unique > 50) {
          return(list(
            valid = FALSE,
            reason_en = paste0("Variable '", var_name, "' appears to be a CONTINUOUS variable (",
                              n_unique, " unique values out of ", n_total, " observations = ",
                              round(pct_unique*100, 1), "%).\n",
                              "\tError() term requires a GROUPING variable (subject ID, block, etc.), not a continuous measure.\n",
                              "\tSolution: Use a categorical grouping variable (subject ID, block) for Error() term."),
            reason_fr = paste0("La variable '", var_name, "' semble être une variable CONTINUE (",
                              n_unique, " valeurs uniques sur ", n_total, " observations = ",
                              round(pct_unique*100, 1), "%).\n",
                              "\tLe terme Error() nécessite une variable de GROUPEMENT (ID sujet, bloc, etc.), pas une mesure continue.\n",
                              "\tSolution : Utiliser une variable de groupement catégorielle (ID sujet, bloc) pour Error().")
          ))
        }

        # Sinon, c'est probablement un ID discret => OK, convertir en facteur
        return(list(valid = TRUE, convert = TRUE, n_groups = n_unique))
      }

      # Autre type => invalide
      return(list(
        valid = FALSE,
        reason_en = paste0("Variable '", var_name, "' has unsupported type for Error() term: ", class(var_data)[1]),
        reason_fr = paste0("La variable '", var_name, "' a un type non supporté pour le terme Error() : ", class(var_data)[1])
      ))
    }

    # Validation et conversion de id (si déjà extrait)
    if (!is.null(id) && !is.null(data) && id %in% colnames(data)) {
      check_result <- .is_valid_grouping_var(data[[id]], id)

      if (!check_result$valid) {
        .exit(
          check_result$reason_en,
          check_result$reason_fr,
          return = return, verbose = verbose
        )
      }

      if (check_result$convert) {
        data[[id]] <- factor(data[[id]])
        .dbg(paste0("Converted '", id, "' to factor for Error() term (", check_result$n_groups, " groups)."),
             paste0("Conversion de '", id, "' en facteur pour Error() (", check_result$n_groups, " groupes)."),
             debug = debug)
      }
    }

    # --- 9.5 : Environnement ---
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
            verbose = verbose, code = code, return = return)
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
        # Extraire le message d'erreur original
        error_msg <- conditionMessage(e)

        # Cas spécial: erreur liée à Error()
        if (grepl("Error", error_msg, ignore.case = FALSE)) {
          .exit(
            paste0("The Error() term in your formula could not be processed.\n",
                   "\tOriginal error: ", error_msg, "\n",
                   "\tPlease ensure:\n",
                   "\t- The variable inside Error() exists in your data\n",
                   "\t- The variable is a grouping factor (subject ID, block, etc.)\n",
                   "\t- Example: Error(subject_id) where subject_id is a column in your data"),
            paste0("Le terme Error() dans votre formule n'a pas pu être traité.\n",
                   "\tErreur originale : ", error_msg, "\n",
                   "\tAssurez-vous que :\n",
                   "\t- La variable dans Error() existe dans vos données\n",
                   "\t- La variable est un facteur de groupement (ID sujet, bloc, etc.)\n",
                   "\t- Exemple : Error(id_sujet) où id_sujet est une colonne de vos données"),
            return = return, verbose = verbose
          )
        }

        # Créer un message plus informatif pour autres erreurs
        .exit(
          paste0("Error evaluating the formula.\n",
                 "\tOriginal error: ", error_msg, "\n",
                 "\tPlease check:\n",
                 "\t- All variables exist in your data\n",
                 "\t- Variable names are spelled correctly\n",
                 "\t- Data types are appropriate for the formula"),
          paste0("Erreur lors de l'évaluation de la formule.\n",
                 "\tErreur originale : ", error_msg, "\n",
                 "\tVérifiez s'il vous plaît :\n",
                 "\t- Toutes les variables existent dans vos données\n",
                 "\t- Les noms de variables sont correctement orthographiés\n",
                 "\t- Les types de données sont appropriés pour la formule"),
          return = return, verbose = verbose
        )
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
            verbose = verbose, code = code, return = return)
    }

    # --- 9.6 : CONSTRUCTION FORMULE POUR ANALYSE (AVEC Error si nécessaire) ---
    formula_for_analysis <- formula

    if (!is.null(id) && !has_error_already) {

      if (!is.null(wt_names)) {
        # Cas avec wt
        for (wt_var in wt_names) {
          if (!wt_var %in% colnames(data)) {
            .exit(paste0("The specified variable `", wt_var, "` does not exist in data."),
                  paste0("La variable spécifiée `", wt_var, "` n'existe pas dans `data`."),
                  verbose = verbose, code = code, return = return)
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
        # Cas sans wt : plan apparié simple
        # Référence académique : Maxwell, Delaney & Kelley (2018), Chapter 11-12
        # Pour un plan apparié (within-subject), chaque id doit apparaître dans TOUTES
        # les modalités du/des facteur(s) intra-sujet.
        #
        # SYNTAXE CORRECTE :
        # - Plan apparié (within-subject) : Error(id)
        #   Chaque sujet traverse tous les niveaux => pas d'emboîtement
        #
        # SYNTAXE INCORRECTE :
        # - Error(id/F) décrit un effet IMBRIQUÉ (nested) : id dans F
        #   Cela impliquerait un plan ENTRE-sujets (between-subject)
        #   où chaque id n'apparaît que dans UN niveau de F
        #   => INCOMPATIBLE avec un plan apparié
        #
        # DONC : Toujours utiliser Error(id) pour un plan apparié,
        # indépendamment du nombre de facteurs intra-sujet

        formula_for_analysis <- update(formula_for_analysis,
                                       paste0(". ~ . + Error(", id, ")"))
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
          verbose = verbose, code = code, return = return)
  }

  lev_g <- .n_levels_g(g)

  if (lev_g > maxcat) {
    .exit(paste0("Too many groups (", lev_g, " > ", maxcat, ")."),
          paste0("Trop de groupes (", lev_g, " > ", maxcat, ")."),
          verbose = verbose, code = code, return = return)
  }

  if (lev_g < 2) {
    .exit("At least 2 groups are required for comparison.",
          "Au moins 2 groupes sont nécessaires pour une comparaison.",
          verbose = verbose, code = code, return = return)
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
            verbose = verbose, code = code, k = k, cpt = "off"
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
  # BLOC 11B : FILTRAGE DES GROUPES AVEC EFFECTIFS INSUFFISANTS
  ################################################################################
  # Référence : Pour des tests statistiques fiables (bartlett.test, bootstrap, etc.),
  # chaque groupe doit avoir au moins min_group_size observations.
  # Les groupes trop petits sont écartés AVANT l'analyse pour éviter les erreurs.

  min_group_size <- 3  # Minimum requis par groupe

  # Identifier les variables factorielles pour le filtrage
  filter_factor_vars <- character(0)

  if (is.data.frame(g)) {
    # Cas multi-facteurs : identifier les colonnes factorielles
    for (col in names(g)) {
      if (is.factor(g[[col]]) || is.character(g[[col]])) {
        filter_factor_vars <- c(filter_factor_vars, col)
      }
    }
  } else if (is.factor(g) || is.character(g)) {
    # Cas un seul facteur : créer un nom temporaire
    filter_factor_vars <- ".group_var"
  }

  # Appliquer le filtrage si nous avons des facteurs
  if (length(filter_factor_vars) > 0) {

    # Préparer data pour le filtrage
    if (!is.null(data)) {
      data_for_filter <- data
    } else {
      # Créer un data.frame temporaire si data est NULL
      data_for_filter <- data.frame(.x_temp = x)
    }

    # Si g est un vecteur simple, l'ajouter temporairement à data
    if (length(filter_factor_vars) == 1 && filter_factor_vars[1] == ".group_var") {
      data_for_filter$.group_var <- g
    }

    # Appeler la fonction de filtrage
    filter_result <- .filter_small_groups(
      data = data_for_filter,
      factor_vars = filter_factor_vars,
      min_n = min_group_size,
      verbose = verbose,
      k = k,
      code = code
    )

    k <- filter_result$k

    # Si des groupes ont été filtrés, mettre à jour les données
    if (filter_result$filtered) {

      if (!is.null(data)) {
        data <- filter_result$data

        # Supprimer la colonne temporaire si créée
        if (".group_var" %in% names(data)) {
          data$.group_var <- NULL
        }

        # Recalculer x et g à partir des données filtrées si formule présente
        if (!is.null(formula)) {
          # Re-extraire x et g depuis les données filtrées
          f_tmp <- .strip_data_dollar_safe(formula, data)
          f_mf  <- .drop_error_term(f_tmp)

          mf <- model.frame(
            f_mf,
            data = data,
            na.action = na.pass,
            drop.unused.levels = TRUE
          )

          x <- model.response(mf)
          pred_cols <- colnames(mf)[-1]

          if (!is.null(id) && id %in% pred_cols) {
            pred_cols <- setdiff(pred_cols, id)
          }

          if (length(pred_cols) == 1) {
            g <- mf[[pred_cols[1]]]
            if (!is.factor(g)) g <- factor(g)
          } else if (length(pred_cols) > 1) {
            g <- mf[, pred_cols, drop = FALSE]
          }
        }
      } else {
        # Sans data, mettre à jour x et g directement depuis filter_result
        filtered_data <- filter_result$data
        x <- filtered_data$.x_temp
        g <- droplevels(filtered_data$.group_var)
      }

      # Recalculer le nombre de niveaux après filtrage
      lev_g <- .n_levels_g(g)

      # Vérifier qu'il reste au moins 2 groupes
      if (lev_g < 2) {
        .exit("After removing groups with insufficient observations, fewer than 2 groups remain.",
              "Après suppression des groupes avec effectifs insuffisants, il reste moins de 2 groupes.",
              verbose = verbose, code = code, return = return)
      }
    }
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

  # --- ROUTAGE SPÉCIAL POUR Error(facteur) ---
  # Idée : Si Error(facteur_catégoriel) détecté (ratio < 5%), forcer vers RM-ANOVA
  # APA : Maxwell & Delaney (2004). Error() spécifie strates → aov() obligatoire
  # DOI : https://doi.org/10.4324/9781315642956
  if (exists("has_error_term_factor") && has_error_term_factor && has_error) {
    .dbg(paste0("Error(categorical_factor) detected with formula: ", deparse(formula), "\n",
                "\tForcing routing to RM-ANOVA analysis (aov with Error term).\n",
                "\tThis is NOT a simple comparison - Error() defines strata structure."),
         paste0("Error(facteur_catégoriel) détecté avec formule : ", deparse(formula), "\n",
                "\tForçage du routage vers analyse RM-ANOVA (aov avec terme Error).\n",
                "\tCe n'est PAS une simple comparaison - Error() définit la structure en strates."),
         debug = debug)

    # Forcer vers .multi_factor_analysis() qui sait gérer aov() avec Error()
    check_anova <- TRUE
    # paired reste FALSE (aov() gère Error() directement)
  }

  # Détection préliminaire des réplicats pour forcer routage vers .multi_factor_analysis()
  # Idée : Si plan apparié avec réplicats, forcer check_anova=TRUE pour routage correct
  # APA : Barr et al. (2013). Random effects structure for confirmatory hypothesis testing.
  #       J Mem Lang 68(3):255-278. DOI: 10.1016/j.jml.2012.11.001
  # Note : RM-ANOVA suppose UNE obs par sujet×condition; avec réplicats => modèle mixte requis
  if (isTRUE(paired) && !is.null(id) && !is.null(data) && id %in% colnames(data)) {
    # Compter observations par sujet
    id_counts <- table(data[[id]])
    n_subjects <- length(id_counts)

    if (is.vector(g) || is.factor(g)) {
      # Cas simple : un facteur
      n_conditions <- length(unique(g))
      expected_obs_per_subject <- n_conditions

      # Si obs par sujet > nombre de conditions => réplicats
      if (all(id_counts > expected_obs_per_subject)) {
        .dbg(paste0("Replicates detected: ", unique(id_counts)[1], " obs per subject, ",
                   n_conditions, " conditions => ", unique(id_counts)[1] / n_conditions, " replicates."),
             paste0("Réplicats détectés: ", unique(id_counts)[1], " obs par sujet, ",
                   n_conditions, " conditions => ", unique(id_counts)[1] / n_conditions, " réplicats."),
             debug = debug)

        # Forcer vers .multi_factor_analysis() pour gestion modèle mixte
        check_anova <- TRUE

        .dbg("Forcing check_anova=TRUE to route to .multi_factor_analysis() for mixed model handling.",
             "Forçage check_anova=TRUE pour router vers .multi_factor_analysis() pour modèle mixte.",
             debug = debug)
      }
    }
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

      .dbg(paste0("Routing: ", route),
           paste0("Routage : ", route),
           debug = debug)

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
              verbose = verbose, code = code, return = return)
      }

      route <- if (lev == 2) "Paired (2 levels)" else "Paired (k>2 levels)"

      .dbg(paste0("Routing: ", route),
           paste0("Routage : ", route),
           debug = debug)

      bilan <- .one_factor_analysis(x = x, g = g, formula = formula, data = data,
                                    alpha = alpha, paired = TRUE, id = id,
                                    k = k, code = code, debug = debug, verbose = verbose)
    }

  } else {

    if (isTRUE(check_manova)) {
      .dbg("Routing: MANOVA",
           "Routage : MANOVA",
           debug = debug)

      bilan <- .manova_analysis(x = x, g = g, formula = formula, data = data,
                                alpha = alpha, paired = FALSE, id = id, between = between,
                                k = k, code = code, debug = debug, verbose = verbose)

    } else if (isTRUE(check_anova)) {
      .dbg("Routing: multi-factor ANOVA",
           "Routage : ANOVA multi-facteurs",
           debug = debug)

      bilan <- .multi_factor_analysis(x = x, g = g, formula = formula, data = data,
                                      alpha = alpha, paired = FALSE, id = id, between = between,
                                      k = k, code = code, debug = debug, verbose = verbose)

    } else {
      .dbg("Routing: one-factor (unpaired)",
           "Routage : 1 facteur (non apparié)",
           debug = debug)

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
  chosen_test <- bilan$chosen_test  # Quel test non-paramétrique a été utilisé (med1way, t1way, kruskal, ou NULL)

  ################################################################################
  # BLOC 16 : GRAPHIQUES (DÉPLACÉ APRÈS POST-HOCS - voir BLOC 18)
  ################################################################################

  # NOTE PERSO: Priorité 5 - Graphiques avec lettres de significativité
  # L'ancien code graphique (boxplot simple) a été déplacé APRÈS les post-hocs
  # pour permettre la superposition des lettres de significativité.
  # Voir BLOC 18 ci-dessous pour la nouvelle implémentation.

  ################################################################################
  # BLOC 17 : POST-HOC
  ################################################################################

  .dbg(NULL, "Tests posts-hocs.", debug = debug)

  # NOTE PERSO: Gestion return=FALSE (Cahier des charges priorité 2)
  # Si return=FALSE → renvoyer UNIQUEMENT p-value globale, PAS de post-hocs
  if (return == FALSE) {
    # Extraire la p-value globale du bilan
    global_pvalue <- bilan$global_pvalue

    # Si global_pvalue manquant, essayer de l'extraire du modèle
    if (is.null(global_pvalue) || is.na(global_pvalue)) {
      .dbg("Warning: global_pvalue not found in bilan, trying to extract from model.",
           "Attention : global_pvalue introuvable dans bilan, extraction depuis le modèle.",
           debug = debug)

      # Tentative d'extraction depuis le modèle
      if (!is.null(bilan$model) && inherits(bilan$model, c("aov", "lm"))) {
        tryCatch({
          anova_table <- car::Anova(bilan$model, type = "III")
          # Prendre la p-value de la première ligne (facteur principal)
          if (nrow(anova_table) >= 1) {
            # Chercher la colonne Pr(>F) ou p.value
            pval_col <- which(colnames(anova_table) %in% c("Pr(>F)", "p.value", "P(>|t|)"))
            if (length(pval_col) > 0) {
              # Exclure les lignes Residuals/Intercept
              row_names <- rownames(anova_table)
              valid_rows <- which(!grepl("Residuals|Intercept|^\\(Intercept\\)", row_names, ignore.case = TRUE))
              if (length(valid_rows) > 0) {
                global_pvalue <- anova_table[valid_rows[1], pval_col[1]]
              }
            }
          }
        }, error = function(e) {
          .dbg(paste0("Could not extract p-value from model: ", e$message),
               paste0("Impossible d'extraire p-value du modèle : ", e$message),
               debug = debug)
        })
      }

      # Si toujours NA, essayer depuis robust_results
      if ((is.null(global_pvalue) || is.na(global_pvalue)) && !is.null(bilan$robust_results)) {
        if (!is.null(bilan$robust_results$p_value)) {
          global_pvalue <- bilan$robust_results$p_value
        }
      }

      # Dernier recours : NA
      if (is.null(global_pvalue) || is.na(global_pvalue)) {
        global_pvalue <- NA
      }
    }
    return(as.numeric(global_pvalue))
  }

  # SINON (return == TRUE): Faire les post-hocs et retourner bilan complet

  # Post-hocs: For MANOVA, use specialized multivariate post-hocs
  if (is.matrix(x)) {
    .dbg("Running MANOVA post-hocs (discriminant analysis + protected ANOVAs).",
         "Exécution des post-hocs MANOVA (analyse discriminante + ANOVAs protégées).",
         debug = debug)

    # Appeler .posthoc_MANOVA() pour analyses post-hoc multivariées appropriées
    posthoc_manova <- tryCatch({
      .posthoc_MANOVA(
        x = x,
        g = g,
        manova_result = bilan,
        alpha = alpha,
        method = "both",  # Discriminant + protected ANOVAs
        verbose = verbose, code = code,
        debug = debug,
        k = k
      )
    }, error = function(e) {
      k <<- .vbse(
        paste0("Note: MANOVA post-hocs skipped due to error: ", e$message),
        paste0("Note: Post-hocs MANOVA ignorés en raison d'une erreur: ", e$message),
        verbose = verbose, code = code, k = k, cpt = "on"
      )
      return(NULL)
    })

    # Ajouter les résultats post-hocs au bilan
    if (!is.null(posthoc_manova)) {
      bilan$posthoc_manova <- posthoc_manova
    }

    # NOTE PERSO: Priorité 5 - Graphics MANOVA
    # Pour MANOVA, créer graphiques séparés pour chaque variable dépendante
    if (plot == TRUE) {
      # Gérer le cas où g est un data.frame (multi-facteurs)
      g_for_plot <- if (is.data.frame(g)) {
        # Créer facteur combiné pour graphique
        interaction(g, drop = TRUE, sep = ":")
      } else {
        g
      }

      n_groups_plot <- length(unique(g_for_plot))

      if (n_groups_plot < maxcat) {
        .dbg("Creating MANOVA plots (one per dependent variable).",
             "Création graphiques MANOVA (un par variable dépendante).",
             debug = debug)

        # Créer un graphique par colonne de x (variable dépendante)
        n_vars <- ncol(x)
        var_names <- colnames(x)
        if (is.null(var_names)) {
          var_names <- paste0("V", seq_len(n_vars))
        }

        # Layout pour plusieurs graphiques (2 colonnes)
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)

        n_rows <- ceiling(n_vars / 2)
        par(mfrow = c(n_rows, 2), mar = c(4, 4, 3, 1))

        for (i in seq_len(n_vars)) {
          # Extraire post-hocs pour cette VD si disponibles
          posthoc_i <- NULL
          if (!is.null(posthoc_manova) && !is.null(posthoc_manova$protected_anovas)) {
            if (!is.null(posthoc_manova$protected_anovas[[var_names[i]]])) {
              posthoc_i <- posthoc_manova$protected_anovas[[var_names[i]]]
            }
          }

          # Appeler .plot_with_letters() pour cette VD
          .plot_with_letters(
            x = x[, i],
            g = g_for_plot,
            posthoc_result = posthoc_i,
            main = paste0("MANOVA - ", var_names[i]),
            ylab = var_names[i],
            xlab = "Groupes",
            verbose = verbose, code = code,
            debug = debug
          )
        }
      } else {
        .dbg(paste0("MANOVA plots skipped: too many groups (", n_groups_plot, " >= ", maxcat, ")"),
             paste0("Graphiques MANOVA ignorés : trop de groupes (", n_groups_plot, " >= ", maxcat, ")"),
             debug = debug)
      }
    }

    # Message de version final (MANOVA)
    if (verbose == TRUE) {
      cat("\n")
      cat(.msg(
        "m.test() - version 02 beta 2025 - report any issues to antoine.masse@u-bordeaux.fr.\n",
        "m.test() - version 02 bêta 2025 - envoyer ce bilan à antoine.masse@u-bordeaux.fr en cas d'anomalie.\n"
      ))
    } else if (code == TRUE && !is.null(bilan)) {
      # En mode code, ajouter le message en commentaire
      cat("# m.test() - version 02 bêta 2025 - envoyer ce bilan à antoine.masse@u-bordeaux.fr en cas d'anomalie.\n")
    }

    # Nettoyage du bilan MANOVA : conserver uniquement les éléments utiles à l'utilisateur
    # Supprimer les éléments internes et restructurer pour cohérence avec autres retours
    if (!is.null(bilan)) {
      # Éléments à conserver pour MANOVA
      manova_result <- list(
        groups = if (!is.null(bilan$posthoc_manova$groups)) bilan$posthoc_manova$groups else NULL,
        p.value = if (!is.null(bilan$posthoc_manova$p.value)) bilan$posthoc_manova$p.value else NULL,
        test_statistics = bilan$test_statistics,
        global_pvalue = bilan$global_pvalue,
        discriminant_analysis = if (!is.null(bilan$posthoc_manova$discriminant_analysis))
                                  bilan$posthoc_manova$discriminant_analysis else NULL,
        protected_anovas = if (!is.null(bilan$posthoc_manova$protected_anovas))
                             bilan$posthoc_manova$protected_anovas else NULL,
        method = "MANOVA (Wilks' Lambda)"
      )
      class(manova_result) <- "posthoc"
      return(manova_result)
    }

    return(bilan)

  } else{
    # Univariate case: check if mixed model, ANCOVA, or regular ANOVA/t-test

    ############################################################################
    # BLOC 17.5 : ANNONCE DES POST-HOCS (Étape 8)
    ############################################################################

    # NOTE PERSO: Routage vers .posthoc_mixed_model() si modèle mixte utilisé
    is_mixed_model <- !is.null(bilan$robust_results) &&
                      !is.null(bilan$robust_results$method) &&
                      bilan$robust_results$method == "Mixed_Model_lmer"

    # NOTE PERSO: Routage vers .posthoc_ANCOVA() si ANCOVA détectée
    # (cf. Cahier des charges priorité 3)
    is_ancova <- !is.null(bilan$check_ancova) && isTRUE(bilan$check_ancova)

    # NOTE BP-008: Afficher le message uniquement pour les ANOVA/t-tests standards avec >2 groupes
    # Pour 2 groupes, la comparaison a déjà été présentée dans l'analyse principale
    # Les fonctions spécialisées (.posthoc_ANCOVA, .posthoc_mixed_model) créent leur propre étape
    # NOTE IMPORTANTE: Dans KefiR, les post-hocs sont TOUJOURS effectués, même si test global non significatif
    # Cela donne une vue complète à l'utilisateur qui peut juger par lui-même
    # Note: Le titre "Posthoc - Tests post-hoc..." est affiché par .posthoc() lui-même

    if (is_mixed_model) {
      .dbg("Mixed model detected, routing to .posthoc_mixed_model()",
           "Modèle mixte détecté, routage vers .posthoc_mixed_model()",
           debug = debug)

      # Extraire le modèle et les effets significatifs
      mixed_model <- bilan$robust_results$model
      sig_effects <- bilan$robust_results$significant_effects

      synth <- .posthoc_mixed_model(
        mixed_model = mixed_model,
        significant_effects = sig_effects,
        alpha = alpha,
        method = "tukey",
        conf.level = conf,
        verbose = verbose, code = code,
        debug = debug,
        k = k
      )

    } else if (is_ancova) {
      .dbg("ANCOVA detected, routing to .posthoc_ANCOVA()",
           "ANCOVA détectée, routage vers .posthoc_ANCOVA()",
           debug = debug)

      # Pour ANCOVA, on a besoin du modèle aov
      # Extraire ou reconstruire le modèle depuis le bilan
      ancova_model <- NULL

      if (!is.null(bilan$model)) {
        ancova_model <- bilan$model
      } else if (!is.null(formula) && !is.null(data)) {
        # Reconstruire le modèle
        ancova_model <- tryCatch({
          aov(formula, data = data)
        }, error = function(e) {
          k <<- .vbse(
            paste0("Warning: Could not reconstruct ANCOVA model for post-hocs: ", e$message),
            paste0("Attention : Impossible de reconstruire le modèle ANCOVA pour post-hocs : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          NULL
        })
      }

      if (!is.null(ancova_model)) {
        synth <- .posthoc_ANCOVA(
          ancova_model = ancova_model,
          factor_name = NULL,  # Auto-détection du premier facteur
          alpha = alpha,
          method = "tukey",  # Défaut académique
          conf.level = conf,
          verbose = verbose,
          debug = debug,
          code = code,
          k = k
        )
      } else {
        # Fallback si modèle non disponible
        k <- .vbse(
          paste0("Warning: ANCOVA model not available for post-hocs.\n",
                 "Falling back to standard post-hoc tests (may be inappropriate for ANCOVA).\n",
                 "Consider using emmeans package directly for proper ANCOVA comparisons."),
          paste0("Attention : Modèle ANCOVA non disponible pour post-hocs.\n",
                 "Utilisation des tests post-hoc standards (peut être inapproprié pour ANCOVA).\n",
                 "Considérez l'usage direct du package emmeans pour comparaisons ANCOVA appropriées."),
          verbose = verbose, code = code, k = k, cpt = "on"
        )

        synth <- tryCatch({
          .posthoc(x, g, alpha = alpha, normal = check_normality,
                   var.equal = check_variance_equal,
                   control = control, debug = debug,
                   verbose = verbose, code = code, paired = paired,
                   boot = boot, boot_type = boot_type, iter = iter, conf = conf, k = k,
                   chosen_test = chosen_test)
        }, error = function(e) {
          k <<- .vbse(
            paste0("Note: Post-hoc comparisons skipped: ", e$message),
            paste0("Note : Comparaisons post-hoc ignorées : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          return(NULL)
        })
      }
    } else {
      # Standard ANOVA/t-test post-hocs
      .dbg("Using standard post-hoc tests (.posthoc)",
           "Utilisation des tests post-hoc standards (.posthoc)",
           debug = debug)

      # Appeler .posthoc() même pour 2 groupes (en mode silencieux) pour obtenir
      # l'objet complet avec bootstrap, p-values, etc.
      # Les messages ont déjà été affichés par .one_factor_analysis()
      number_of_groups <- length(unique(g))

      if (number_of_groups == 2) {
        # Pour 2 groupes, les comparaisons sont déjà affichées dans .one_factor_analysis()
        # MAIS on appelle quand même .posthoc() en mode silencieux (verbose=FALSE)
        # pour obtenir l'objet complet avec bootstrap, p-values, etc.
        .dbg("Calling .posthoc() for 2 groups in silent mode (messages already displayed in .one_factor_analysis())",
             "Appel de .posthoc() pour 2 groupes en mode silencieux (messages déjà affichés dans .one_factor_analysis())",
             debug = debug)

        synth <- tryCatch({
          .posthoc(x, g, alpha = alpha, normal = check_normality,
                   var.equal = check_variance_equal,
                   control = control, code = FALSE, debug = debug,
                   verbose = FALSE,  # MODE SILENCIEUX pour éviter réaffichage (et code=FALSE pour éviter duplication)
                   paired = paired,
                   boot = boot, boot_type = boot_type, iter = iter, conf = conf, k = k,
                   chosen_test = chosen_test)
        }, error = function(e) {
          k <<- .vbse(
            paste0("Note: Post-hoc analysis failed for 2 groups: ", e$message),
            paste0("Note : Analyse post-hoc échouée pour 2 groupes : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          return(NULL)
        })
      } else {
        # Plus de 2 groupes : appeler .posthoc() TOUJOURS (même si test global non significatif)
        synth <- tryCatch({
          .posthoc(x, g, alpha = alpha, normal = check_normality,
                   var.equal = check_variance_equal,
                   control = control, debug = debug,
                   verbose = verbose, code = code, paired = paired,
                   boot = boot, boot_type = boot_type, iter = iter, conf = conf, k = k,
                   chosen_test = chosen_test)
        }, error = function(e) {
          k <<- .vbse(
            paste0("Note: Post-hoc comparisons skipped: ", e$message),
            paste0("Note : Comparaisons post-hoc ignorées : ", e$message),
            verbose = verbose, code = code, k = k, cpt = "on"
          )
          return(NULL)
        })
      }

      # Synchroniser global_pvalue pour le retour
      # .posthoc() définit synth$p.value mais pas synth$global_pvalue
      # IMPORTANT: global_pvalue doit être un SCALAIRE, pas une matrice
      if (!is.null(synth) && is.null(synth$global_pvalue)) {
        if (!is.null(synth$p.value)) {
          # Si p.value est une matrice, extraire le minimum (hors diagonale)
          if (is.matrix(synth$p.value)) {
            pmat <- synth$p.value
            diag(pmat) <- NA
            synth$global_pvalue <- min(pmat, na.rm = TRUE)
          } else if (length(synth$p.value) == 1) {
            synth$global_pvalue <- as.numeric(synth$p.value)
          } else {
            # Vecteur de p-values : prendre le minimum
            synth$global_pvalue <- min(synth$p.value, na.rm = TRUE)
          }
        } else if (!is.null(pvals) && !is.na(pvals)) {
          synth$global_pvalue <- as.numeric(pvals)
        }
      }
    }

    ############################################################################
    # BLOC 18 : GRAPHIQUES AVEC LETTRES (Priorité 5)
    ############################################################################

    # NOTE PERSO: Graphics placés APRÈS post-hocs pour superposer lettres
    # Univariate case: un seul graphique avec lettres de significativité

    if (plot == TRUE) {
      # Gérer le cas où g est un data.frame (multi-facteurs)
      g_for_plot <- if (is.data.frame(g)) {
        # Créer facteur combiné pour graphique
        interaction(g, drop = TRUE, sep = ":")
      } else {
        g
      }

      n_groups_plot <- length(unique(g_for_plot))

      if (n_groups_plot < maxcat) {
        .dbg("Creating plot with significance letters.",
             "Création graphique avec lettres de significativité.",
             debug = debug)

        # Déterminer titre selon type d'analyse
        plot_title <- ""
        if (is_ancova) {
          plot_title <- "ANCOVA (adjusted means)"
        } else if (paired) {
          plot_title <- "Paired comparisons"
        } else if (is_mixed_model) {
          plot_title <- "Mixed Model (EMMs)"
        } else {
          plot_title <- "ANOVA"
        }

        # Gérer format multi-effets de .posthoc_mixed_model()
        # Si synth contient $note et plusieurs sous-listes : extraire le premier effet
        posthoc_for_plot <- synth
        if (!is.null(synth$note) && !inherits(synth, "posthoc")) {
          # Multi-effets : extraire le premier effet significatif
          effect_names <- setdiff(names(synth), c("k", "note"))
          if (length(effect_names) > 0) {
            posthoc_for_plot <- synth[[effect_names[1]]]
            plot_title <- paste0(plot_title, " - ", effect_names[1])
          }
        }

        # Appeler .plot_with_letters() avec résultats post-hocs
        .plot_with_letters(
          x = x,
          g = g_for_plot,
          posthoc_result = posthoc_for_plot,
          main = plot_title,
          ylab = "Value",
          xlab = "Groups",
          verbose = verbose, code = code,
          debug = debug
        )
      } else {
        k <- .vbse(
          paste0("Plot skipped: too many groups (", n_groups_plot, " >= ", maxcat, ")"),
          paste0("Graphique ignoré : trop de groupes (", n_groups_plot, " >= ", maxcat, ")"),
          verbose = verbose, code = code, k = k, cpt = "on"
        )
      }
    }

    # Ajouter chosen_test à l'objet synth pour le retour
    if (!is.null(synth)) {
      synth$chosen_test <- chosen_test  # NULL=parametric, "med1way", "t1way", or "kruskal"

      # ========================================================================
      # STANDARDISATION DU RETOUR : Ajouter descriptive_stats si absent
      # ========================================================================
      # Structure standard : groups, p.value, global_pvalue, descriptive_stats, method

      # 1. Ajouter global_pvalue si absent
      if (is.null(synth$global_pvalue)) {
        if (!is.null(bilan$global_pvalue)) {
          synth$global_pvalue <- bilan$global_pvalue
        } else if (!is.null(synth$p.value) && is.matrix(synth$p.value)) {
          # Prendre la plus petite p-value non-diagonale
          pmat <- synth$p.value
          diag(pmat) <- NA
          synth$global_pvalue <- min(pmat, na.rm = TRUE)
        }
      }

      # 2. Créer descriptive_stats si absent (harmonisation)
      if (is.null(synth$descriptive_stats)) {
        # Déterminer le type de statistique selon le test utilisé
        stat_type <- if (!is.null(chosen_test) && chosen_test %in% c("med1way", "kruskal")) {
          "median"  # Tests non-paramétriques
        } else {
          "mean"    # Tests paramétriques
        }

        # Si ANCOVA avec adjusted_means, utiliser celles-ci
        if (!is.null(synth$adjusted_means)) {
          synth$descriptive_stats <- data.frame(
            group = synth$adjusted_means$group,
            adjusted_mean = synth$adjusted_means$adjusted_mean,
            SE = synth$adjusted_means$SE,
            type = "adjusted_mean"
          )
        } else if (!is.data.frame(g)) {
          # Calculer statistiques descriptives par groupe
          g_factor <- droplevels(factor(g))
          levels_g <- levels(g_factor)

          if (stat_type == "median") {
            stats_vals <- tapply(x, g_factor, median, na.rm = TRUE)
            stats_iqr <- tapply(x, g_factor, IQR, na.rm = TRUE)
            synth$descriptive_stats <- data.frame(
              group = levels_g,
              median = as.numeric(stats_vals[levels_g]),
              IQR = as.numeric(stats_iqr[levels_g]),
              type = "median"
            )
          } else {
            stats_mean <- tapply(x, g_factor, mean, na.rm = TRUE)
            stats_sd <- tapply(x, g_factor, sd, na.rm = TRUE)
            stats_n <- tapply(x, g_factor, function(z) sum(!is.na(z)))
            stats_se <- stats_sd / sqrt(stats_n)
            synth$descriptive_stats <- data.frame(
              group = levels_g,
              mean = as.numeric(stats_mean[levels_g]),
              SD = as.numeric(stats_sd[levels_g]),
              SE = as.numeric(stats_se[levels_g]),
              n = as.numeric(stats_n[levels_g]),
              type = "mean"
            )
          }
        }
      }

      # 3. Ajouter method si absent
      if (is.null(synth$method)) {
        synth$method <- if (!is.null(chosen_test)) {
          switch(chosen_test,
            "med1way" = "Robust ANOVA (med1way)",
            "t1way" = "Robust ANOVA (t1way)",
            "kruskal" = "Kruskal-Wallis",
            "Parametric ANOVA/t-test"
          )
        } else {
          "Parametric ANOVA/t-test"
        }
      }

      # Nettoyage des éléments internes (pas utiles à l'utilisateur)
      synth$k <- NULL           # Compteur interne .vbse()
      synth$verbose <- NULL     # Paramètre interne
      synth$debug <- NULL       # Paramètre interne
      synth$code <- NULL        # Paramètre interne

      # Supprimer adjusted_means du niveau principal (déjà dans descriptive_stats)
      # Garder dans $details pour utilisateurs avancés si besoin
      if (!is.null(synth$adjusted_means) && !is.null(synth$descriptive_stats)) {
        if (is.null(synth$details)) synth$details <- list()
        synth$details$adjusted_means <- synth$adjusted_means
        synth$adjusted_means <- NULL
      }

      # Supprimer les éléments ANCOVA spécifiques du niveau principal
      # (pairwise_comparisons, emm_result, pairs_result, note → déjà dans $details)
      synth$pairwise_comparisons <- NULL
      synth$emm_result <- NULL
      synth$pairs_result <- NULL
      synth$note <- NULL

      # ========================================================================
      # RÉORDONNER les éléments pour structure cohérente
      # ========================================================================
      # Ordre standard : groups, p.value, global_pvalue, descriptive_stats, method, [bootstrap], [details], [chosen_test]

      synth_ordered <- list(
        groups = synth$groups,
        p.value = synth$p.value,
        global_pvalue = synth$global_pvalue,
        descriptive_stats = synth$descriptive_stats,
        method = synth$method
      )

      # Ajouter éléments optionnels s'ils existent
      if (!is.null(synth$bootstrap)) synth_ordered$bootstrap <- synth$bootstrap
      if (!is.null(synth$details)) synth_ordered$details <- synth$details
      if (!is.null(synth$chosen_test)) synth_ordered$chosen_test <- synth$chosen_test

      # Conserver la classe
      class(synth_ordered) <- class(synth)
      synth <- synth_ordered
    }

    # Message de version final
    if (verbose == TRUE) {
      cat("\n")
      cat(.msg(
        "m.test() - version 02 beta 2025 - report any issues to antoine.masse@u-bordeaux.fr.\n",
        "m.test() - version 02 bêta 2025 - envoyer ce bilan à antoine.masse@u-bordeaux.fr en cas d'anomalie.\n"
      ))
    } else if (code == TRUE && !is.null(synth)) {
      # En mode code, ajouter le message en commentaire
      cat("# m.test() - version 02 bêta 2025 - envoyer ce bilan à antoine.masse@u-bordeaux.fr en cas d'anomalie.\n")
    }

    return(synth)
  }
}

################################################################################
#                            FIN DE m.test() v18
################################################################################

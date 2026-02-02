#' Function to validate a regression model
#'
#' @author Antoine MASSE (2025)
#'
#' @param reg A regression model
#' @param verbose To see the detailed balance sheet
#' @param nvar The maximum number of variables allowed
#' @param boot For checking the model by bootstrap with bootreg() (FALSE or TRUE)
#' @param alpha Maximum value accepted for the p-values of the model and its coefficients
#' @param conf.level Confidence interval accepted to validate the regression model by bootstrap
#' @param plot For seeing the graphical analysis of bootreg() (FALSE or TRUE)
#' @param data optional, the data.frame of data if complex model.
#' @param raintest_alpha Minimal value of p-value accepted for Rainbow test
#' @param dwtest_alpha Minimal value of p-value accepted for Durbin-Watson test
#' @param shapiro_alpha Minimal value of p-value accepted for Shapiro-Wilk test
#' @param bptest_alpha Minimal value of p-value accepted for Breush-Pagan test
#' @param k Integer. Message counter for integration with m.test() (default: 0)
#' @param orderDW Character. Name of variable for ordering data before Durbin-Watson test. If NULL, uses current order with warning (default: NULL)
#' @param tolerance Character. Normality tolerance level: "basic" or "extrem" following .normality() logic (default: "basic")
#' @param debug Logical. If TRUE, enables debug mode for troubleshooting.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{valid}: Logical. TRUE if model passes all validation checks, FALSE otherwise
#'   \item \code{k}: Integer. Updated message counter for continued numbering in m.test() pipeline
#' }
#' @details This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).
#' valreg will therefore validate the regression model, control the p-values and, possibly (boot argument), control the reliability by bootstrap with the bootreg function.
#'
#' @section Normality Testing:
#' Now uses .normality() function for consistency with m.test() pipeline. Supports two tolerance levels:
#' \itemize{
#'   \item \code{tolerance="basic"}: Standard criteria (|Skewness| ≤ 1, |Kurtosis| ≤ 1.5)
#'   \item \code{tolerance="extrem"}: Relaxed criteria (|Skewness| ≤ 2, |Kurtosis| ≤ 7) for automatic use in mixed models
#' }
#'
#' @section Durbin-Watson Test:
#' The \code{orderDW} parameter controls how autocorrelation is assessed:
#' \itemize{
#'   \item If \code{NULL} (default): Test is performed but result is informative only (does not invalidate model)
#'   \item If specified: Test uses the named variable for ordering and result affects validation
#' }
#' @import lmtest
#' @import tseries
#' @import lmerTest
#' @importFrom stats cooks.distance
#' @importFrom lme4 ranef
#' @importFrom car vif
#' @importFrom graphics plot
#' @importFrom MHTdiscrete Sidak.p.adjust
#' @importFrom MASS fitdistr
#' @importFrom stats as.formula model.frame anova AIC resid fitted lm pf shapiro.test ks.test sd predict
#' @export
#'
#' @examples
#' # Example 1: Linear model (basic usage)
#' data(iris)
#' reg <- lm(Sepal.Length~.,data=iris[,1:4])
#' result <- valreg(reg, verbose=TRUE)
#' print(result$valid)  # TRUE if model is valid
#'
#' # Example 2: Mixed model with extreme tolerance
#' library(lme4)
#' data(sleepstudy)
#' reg_mixed <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' result <- valreg(reg_mixed, verbose=TRUE, tolerance="extrem")
#' print(result$valid)
#'
#' # Example 3: Integration with m.test() pipeline
#' library(lme4)
#' data(iris)
#' reg_mixed <- lmer(Sepal.Length ~ Petal.Length + (1|Species), data = iris)
#' result <- valreg(reg_mixed, verbose=TRUE, k=5)  # Start counter at 5
#' print(result$k)  # Updated counter value
valreg <- function(reg, verbose = TRUE, nvar = 5, boot = TRUE, alpha = 0.05,
                   conf.level = 0.95, plot = FALSE, data = c(),
                   raintest_alpha = 0.05, dwtest_alpha = 0.03,
                   shapiro_alpha = 0.05, bptest_alpha = 0.05,
                   k = 0, orderDW = NULL, tolerance = "basic", debug = FALSE) {
  error <- TRUE
  if (k != 0) {exit_k <- TRUE
  } else {exit_k <- FALSE}
  counter <- k  # Start counting from k instead of 1

  locale <- Sys.getlocale("LC_MESSAGES")
  lang <- ifelse(grepl("fr", locale), "fr", "en")
  if (is.null(nvar)) {
	nvar <- sqrt(length(resid(reg)))

  }
  lrt <- function(model) {
    formula <- formula(model)
    response <- as.character(formula[[2]])
    terms <- attr(terms(formula), "term.labels")
    random_effects <- grep("\\|", terms, value = TRUE)
    if (length(random_effects) == 0) {
      stop(.msg("The provided model has no random effects.", "Le modèle fourni n'a pas d'effets aléatoires."))
    }
    random_effects_str <- paste("(1 |", gsub(".*\\|", "", random_effects), ")", sep = "")
    formula_nul <- as.formula(paste(response, "~", paste(random_effects_str, collapse = " + ")))
    data_used <- model.frame(model)
    # Supprimer warnings boundary singular fit
    model_complet <- suppressWarnings(suppressMessages(
      lmerTest::lmer(formula, data = data_used, REML = TRUE)
    ))
    model_nul <- suppressWarnings(suppressMessages(
      lmerTest::lmer(formula_nul, data = data_used, REML = TRUE)
    ))
    lrt_result <- suppressMessages(anova(model_nul, model_complet))
    p_value <- lrt_result$`Pr(>Chisq)`[2]
    return(p_value)
  }

  lrt2 <- function(model) {
    formula <- formula(model)
    response <- as.character(formula[[2]])
    terms <- attr(terms(formula), "term.labels")
    random_effects <- grep("\\|", terms, value = TRUE)
    if (length(random_effects) == 0) {
      stop(.msg("The provided model has no random effects.", "Le modèle fourni n'a pas d'effets aléatoires."))
    }
    fixed_effects <- setdiff(terms, random_effects)
    if (length(fixed_effects) == 0) {
      stop(.msg("The provided model has no fixed effects.", "Le modèle fourni n'a pas d'effets fixes."))
    }
    fixed_effects_str <- paste(fixed_effects, collapse = " + ")
    formula_lm <- as.formula(paste(response, "~", fixed_effects_str))
    data_used <- model.frame(model)
    # Supprimer warnings boundary singular fit
    model_complet <- suppressWarnings(suppressMessages(
      lmerTest::lmer(formula, data = data_used, REML = TRUE)
    ))
    model_lm <- lm(formula_lm, data = data_used)
    aic_complet <- AIC(model_complet)
    aic_lm <- AIC(model_lm)
    delta_aic <- aic_complet - aic_lm
    return(delta_aic)
  }

  dwtest_mix <- function(model) {
    residuals <- resid(model)
    data_used <- model.frame(model)
    results <- list()
    fitted_values <- fitted(model)
    sorted_residuals_fitted <- residuals[order(fitted_values)]
    lm_model_fitted <- lm(sorted_residuals_fitted ~ 1)
    results$fitted <- dwtest(lm_model_fitted)$p.value
    random_effects <- ranef(model)
    random_effects_vars <- names(random_effects)
    for (var in random_effects_vars) {
      if (var %in% names(data_used)) {
        sorted_residuals_var <- residuals[order(data_used[[var]])]
        lm_model_var <- lm(sorted_residuals_var ~ 1)
        results[[var]] <- dwtest(lm_model_var)$p.value
      } else {
        warning(.msg(paste("Variable", var, "does not exist in the data."),
                    paste("Variable", var, "n'existe pas dans les données.")))
      }
    }
    if (length(random_effects_vars) > 1) {
      valid_vars <- random_effects_vars[random_effects_vars %in% names(data_used)]
      if (length(valid_vars) > 0) {
        sorted_residuals_crossed <- residuals[do.call(order, data_used[valid_vars])]
        lm_model_crossed <- lm(sorted_residuals_crossed ~ 1)
        results$crossed <- dwtest(lm_model_crossed)$p.value
      } else {
        warning(.msg("No valid random effect variables found for crossed sorting.",
                    "Aucune variable d'effet aléatoire valide trouvée pour le tri croisé."))
      }
    }
    return(unlist(results))
  }

  bptest_lmer <- function(model) {
    if (!inherits(model, "lmerMod")) {
      stop(.msg("The provided model is not an lmer model.", "Le modèle fourni n'est pas un modèle lmer."))
    }

    # Approche alternative: utiliser les résidus et fitted values directement
    # au lieu de recréer un lm() qui peut avoir des problèmes de scope
    residuals_lmer <- residuals(model)
    fitted_lmer <- fitted(model)

    # Créer un data.frame temporaire avec les valeurs nécessaires
    temp_data <- data.frame(
      resid_sq = residuals_lmer^2,
      fitted = fitted_lmer
    )

    # Test de Breusch-Pagan sur régression résidus² ~ fitted values
    # C'est équivalent au BP test classique
    temp_model <- lm(resid_sq ~ fitted, data = temp_data)
    bp_test <- bptest(temp_model)

    return(bp_test)
  }

  cooks_distance_lmer <- function(model) {
    # Tentative d'utilisation directe de cooks.distance()
    tryCatch({
      cooksd <- cooks.distance(model)
      return(cooksd)
    }, error = function(e) {
      # En cas d'erreur, appliquer la solution manuelle
      message("⚠️  cooks.distance() a échoué. Utilisation d'une solution manuelle robuste.")

      # Résidus standardisés (resid type pearson)
      resid_std <- residuals(model, type = "pearson") / sigma(model)

      # Matrice de design des effets fixes
      X <- getME(model, "X")

      # Calcul des leviers (diag(X * (X'X)^-1 * X'))
      XtX_inv <- solve(crossprod(X))   # (X'X)^-1
      leverage <- rowSums((X %*% XtX_inv) * X)

      # Éviter les valeurs trop proches de 1 pour les leviers
      leverage[leverage >= 1] <- 0.9999

      # Calcul des distances de Cook
      cooks_d <- (resid_std^2 * leverage) / (1 - leverage)^2 / 2

      return(cooks_d)
    })
  }

  count_terms_excluding_intercept <- function(model) {
    if (!inherits(model, c("lmerMod", "lm"))) {
      stop(.msg("The model must be of type lmerMod (fitted with lme4) or lm.",
               "Le modèle doit être de type lmerMod (ajusté avec lme4) ou lm."))
    }
    terms_model <- terms(model)
    term_labels <- attr(terms_model, "term.labels")
    if (inherits(model, "lmerMod")) {
      fixed_effects <- names(fixef(model))
      term_labels <- term_labels[term_labels %in% fixed_effects]
    }
    num_terms <- length(term_labels)
    return(num_terms)
  }

  norm_lmer <- function(model) {
    random_effects <- ranef(model)
    p_values <- c()
    for (group in names(random_effects)) {
      for (effect in names(random_effects[[group]])) {
        re_values <- random_effects[[group]][[effect]]
        n <- length(re_values)
        if (n <= 50) {
          test_result <- shapiro.test(re_values)
        } else if (n <= 1000) {
          test_result <- jarque.bera.test(re_values)
        } else {
          test_result <- ks.test(re_values, "pnorm", mean = mean(re_values), sd = sd(re_values))
        }
        p_values <- c(p_values, test_result$p.value)
      }
    }
    if (length(p_values) > 1) {
      corrected_p_values <- MHTdiscrete::Sidak.p.adjust(p_values)
    } else {
      corrected_p_values <- p_values
    }
    return(min(corrected_p_values))
  }

  white <- function(model) {
    # Obtenir les coefficients de manière compatible S3/S4
    model_coefs <- if (inherits(model, "lmerMod")) {
      fixef(model)
    } else {
      model$coefficients
    }

    if (any(is.na(model_coefs))) {
      stop(.msg("NA in coefficients.", "NA dans les coefficients."))
    }

    # Obtenir les données de manière compatible S3/S4
    data_used <- if (inherits(model, "lmerMod")) {
      model.frame(model)
    } else {
      model$model
    }

    terms <- attr(terms(model), "term.labels")
    form <- formula(model)
    if (length(terms) > 1) {
      formula <- paste0(deparse(form), "+", paste(sapply(terms, function(x) paste0("I(", x, "^2)")), collapse = "+"), "+", paste(combn(terms, 2, function(x) paste0("I(", x[1], "*", x[2], ")")), collapse = "+"))
    } else {
      formula <- paste0(deparse(form), "+", paste0("I(", terms, "^2)"))
    }
    model <- lm(as.formula(formula), data = data_used)
    return(bptest(model))
  }
		.dbg(NULL, "valreg() - analyse du nombre de variables.",debug=debug)
  nvar1 <- length(coef(reg))
  if (nvar1 > nvar) {
    .vbse(ang=paste0("More than ", nvar, " variables"),
          fr=paste0("Plus de ", nvar, " variables"),
          k=counter, cpt="off", verbose=verbose)
    error <- FALSE
  } else if (length(fitted(reg)) < (nvar + 1)) {
    .vbse(ang="Not enough values in the subset",
          fr="Pas assez de valeurs dans le sous-ensemble",
          k=counter, cpt="off", verbose=verbose)
    error <- FALSE
  } else if (inherits(reg, "lm") && (length(summary(reg)$coefficients[, "Pr(>|t|)"])) != (length(coef(reg)))) {
    error <- FALSE
  } else {
    counter <- .vbse(ang="Analysis of the p-values of the model and its coefficients.",
                     fr="Analyse des p-values du modèle et de ses coefficients.",
                     k=counter, cpt="on", verbose=verbose)
    if (inherits(reg, "lm")) {
      pval_mdl <- pf(summary(reg)$fstatistic[1], summary(reg)$fstatistic[2], summary(reg)$fstatistic[3], lower.tail = FALSE)
      if (pval_mdl > alpha) {
        .vbse(ang=paste0("Warning!\n\tBad significance of the model. p-value: ", .format_pval(pval_mdl)),
              fr=paste0("Attention !\n\tMauvaise signification du modèle. p-value : ", .format_pval(pval_mdl)),
              k=counter, cpt="off", verbose=verbose)
        error <- FALSE
      } else {
        .vbse(ang=paste0("Good significance of the model. p-value: ", .format_pval(pval_mdl)),
              fr=paste0("Bonne signification du modèle. p-value : ", .format_pval(pval_mdl)),
              k=counter, cpt="off", verbose=verbose)
      }
    } else if (inherits(reg, "lmerMod")) {
      .vbse(ang="Please note, mixed model.\n\t==> Calculation of a global p-value for comparison of the model with random effects alone by Likelihood Ratio Test (LRT).",
            fr="Veuillez noter, modèle mixte.\n\t==> Calcul d'une p-value globale pour la comparaison du modèle avec effets aléatoires seuls par test du rapport de vraisemblance (LRT).",
            k=counter, cpt="off", verbose=verbose)
      pval_mdl <- lrt(reg)
      if (pval_mdl > alpha) {
        .vbse(ang=paste0("Warning!\n\tThe fixed effect adds nothing to the random effect. p-value: ", .format_pval(pval_mdl)),
              fr=paste0("Attention !\n\tL'effet fixe n'ajoute rien à l'effet aléatoire. p-value : ", .format_pval(pval_mdl)),
              k=counter, cpt="off", verbose=verbose)
        error <- FALSE
      } else {
        .vbse(ang=paste0("The fixed effect combined with the random effect is more relevant than the random effect alone. p-value: ", .format_pval(pval_mdl)),
              fr=paste0("L'effet fixe combiné à l'effet aléatoire est plus pertinent que l'effet aléatoire seul. p-value : ", .format_pval(pval_mdl)),
              k=counter, cpt="off", verbose=verbose)
      }
      aic_mdl <- lrt2(reg)
      if (aic_mdl > 0) {
        .vbse(ang="Warning!\n\tThe random effect does not improve the model compared to fixed effects alone.",
              fr="Attention !\n\tL'effet aléatoire n'améliore pas le modèle par rapport aux effets fixes seuls.",
              k=counter, cpt="off", verbose=verbose)
        error <- FALSE
      } else {
        .vbse(ang="The random effect improves the model compared to fixed effects alone.",
              fr="L'effet aléatoire améliore le modèle par rapport aux effets fixes seuls.",
              k=counter, cpt="off", verbose=verbose)

        # Tenter ranova() avec gestion d'erreur robuste (problèmes de scope possibles)
        pval_of_ranova <- tryCatch({
          # Recréer le modèle avec les données du model.frame pour éviter problèmes scope
          formula_orig <- formula(reg)
          data_used <- model.frame(reg)

          # Refitter le modèle avec données explicites
          model_refit <- suppressWarnings(suppressMessages(
            lmerTest::lmer(formula_orig, data = data_used, REML = TRUE)
          ))

          # Appliquer ranova() sur le modèle refitté
          ranova_result <- suppressMessages(suppressWarnings(ranova(model_refit)))
          max(ranova_result$`Pr(>Chisq)`, na.rm = TRUE)
        }, error = function(e) {
          # Si ranova() échoue malgré tout, ne pas bloquer la validation
          # Note: Message supprimé pour ne pas encombrer la sortie
          return(0.5)  # Valeur neutre > alpha - pas d'échec de validation
        })

        if (pval_of_ranova < alpha) {
          .vbse(ang="\tConfirmation by ranova().",
                fr="\tConfirmation par ranova().",
                k=counter, cpt="off", verbose=verbose)
        } else {
          .vbse(ang="\tNon-confirmation by ranova().",
                fr="\tNon-confirmation par ranova().",
                k=counter, cpt="off", verbose=verbose)
          error <- FALSE
        }
      }
    }

	# Fonction pour obtenir les p-values des coefficients
	get_pvals <- function(model) {
		summary_lmerModLmerTest <- lmerTest:::summary.lmerModLmerTest
		if (inherits(reg, "lm")) {
		  pvals <- summary(model)$coefficients[, "Pr(>|t|)"]
		} else {
		    # Supprimer les messages de coercition
			suppressMessages({
				model <- lmerTest:::as_lmerModLmerTest(model)
			})
			pvals <- summary_lmerModLmerTest(model)$coefficients[, "Pr(>|t|)"]
		}
		return(pvals)
	}

    pval_coeff <- get_pvals(reg)
    if (max(pval_coeff) > alpha) {
      .vbse(ang=paste0("\tWarning!\n\tBad significance of the coefficients. max(p.value): ", .format_pval(max(pval_coeff))),
            fr=paste0("\tAttention !\n\tMauvaise signification des coefficients. max(p.value) : ", .format_pval(max(pval_coeff))),
            k=counter, cpt="off", verbose=verbose)
      error <- FALSE
    } else {
      .vbse(ang=paste0("Good significance of the coefficients. max(pval_coeff): ", .format_pval(max(pval_coeff))),
            fr=paste0("Bonne signification des coefficients. max(pval_coeff) : ", .format_pval(max(pval_coeff))),
            k=counter, cpt="off", verbose=verbose)
    }

    # =========================================================================
    # NORMALITY CHECK - Now using .normality() for consistency with m.test()
    # =========================================================================
    normality_result <- .normality(
      x = resid(reg),
      g = NULL,  # Residuals are tested as single group
      alpha = shapiro_alpha,
      tolerance = tolerance,
      paired = FALSE,
      debug = FALSE,
      verbose = verbose, code = code,
      k = counter,
      cpt = "on"
    )

    check_normality <- normality_result[[1]]
    counter <- normality_result[[2]]

    if (!check_normality) {
      error <- FALSE
    }

    # Additional check for random effects in mixed models
    if (inherits(reg, "lmerMod")) {
      pvalt <- norm_lmer(reg)
      if (pvalt < shapiro_alpha) {
        counter <- .vbse(ang=paste0("Warning!\n\tNon-normality of random effects for at least one random variable. p.value: ", .format_pval(pvalt)),
                        fr=paste0("Attention !\n\tNon-normalité des effets aléatoires pour au moins une variable aléatoire. p.value : ", .format_pval(pvalt)),
                        k=counter, cpt="on", verbose=verbose)
        error <- FALSE
      } else {
        counter <- .vbse(ang=paste0("Normality of random effects for at least one random variable. p.value: ", .format_pval(pvalt)),
                        fr=paste0("Normalité des effets aléatoires pour au moins une variable aléatoire. p.value : ", .format_pval(pvalt)),
                        k=counter, cpt="on", verbose=verbose)
      }
    }

    if (inherits(reg, "lm")) {
      counter <- .vbse(ang="Analysis of the adequacy of model (Equivalence between the global model and the model established on the best points).",
                      fr="Analyse de l'adéquation du modèle...\t... (équivalence entre le modèle global et le modèle établi sur les meilleurs points).",
                      k=counter, cpt="on", verbose=verbose)
      if ("(Intercept)" %in% names(coef(reg))) {
        nvar2 <- length(coef(reg)) - 1
      } else {
        nvar2 <- length(coef(reg))
      }
      if (nvar2 > 1) {
        raintest(reg, order.by = "mahalanobis")$p.value -> pvalt
        if (pvalt < raintest_alpha) {
          .vbse(ang=paste0("Warning!\n\tRainbow test ordered by mahalanobis (raintest()) - Bad adequacy. p.value: ", .format_pval(pvalt)),
                fr=paste0("Attention !\n\tTest de Rainbow ordonné par mahalanobis (raintest()  de {lmtest}) - Mauvaise adéquation. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
          error <- FALSE
        } else {
          .vbse(ang=paste0("Rainbow test (raintest()) - Good adequacy. p.value: ", .format_pval(pvalt)),
                fr=paste0("Test de Rainbow (raintest()  de {lmtest}) - Bonne adéquation. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
        }
      } else {
        if (length(data) == 0) {
          data <- eval(getCall(reg)$data,environment(formula(reg)))
        }
        z <- as.formula(paste0("~", attr(reg$terms, "term.labels")))
        raintest(reg, order.by = z, data = data)$p.value -> pvalt
        if (pvalt < raintest_alpha) {
          .vbse(ang=paste0("Warning!\n\tRainbow test ordered by X (raintest()) - Bad adequacy. p.value: ", .format_pval(pvalt)),
                fr=paste0("Attention !\n\tTest de Rainbow ordonné par X (raintest() de {lmtest}) - Mauvaise adéquation. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
          error <- FALSE
        } else {
          .vbse(ang=paste0("Rainbow test (raintest()) - Good adequacy. p.value: ", .format_pval(pvalt)),
                fr=paste0("Test de Rainbow (raintest() de {lmtest}) - Bonne adéquation. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
        }
      }
    }

    # =========================================================================
    # DURBIN-WATSON TEST - orderDW parameter for data ordering
    # =========================================================================
    counter <- .vbse(ang="Analysis of independence of the residuals",
                    fr="Analyse de l'indépendance des résidus",
                    k=counter, cpt="on", verbose=verbose)

    dw_count_in_error <- FALSE  # Track if DW should invalidate model

    if (inherits(reg, "lm")) {
      if (length(data) > 0) {
        ypred <- predict(reg, data)
      } else {
        dt_temp <- eval(getCall(reg)$data,environment(formula(reg)))
        ypred <- predict(reg, dt_temp)
      }
      try(dwtest(reg, ypred)$p.value) -> pvalt

      if (is.null(orderDW)) {
        # No ordering variable specified - DW result shown but doesn't invalidate
        .vbse(ang="Warning: Durbin-Watson test performed without explicit ordering variable.\n\tConsider specifying 'orderDW' parameter for meaningful autocorrelation check.\n\tResult shown but NOT used for model validation.",
              fr="Attention : Test de Durbin-Watson effectué sans variable de tri explicite.\n\tConsidérez spécifier le paramètre 'orderDW' pour un contrôle significatif d'autocorrélation.\n\tRésultat affiché mais NON pris en compte pour validation du modèle.",
              k=counter, cpt="off", verbose=verbose)
        .vbse(ang=paste0("Durbin-Watson test (dwtest()) - p.value (informative only): ", .format_pval(pvalt)),
              fr=paste0("Test de Durbin-Watson (dwtest()) - p.value (informatif seulement) : ", .format_pval(pvalt)),
              k=counter, cpt="off", verbose=verbose)
        dw_count_in_error <- FALSE  # Don't count this in error
      } else {
        # Ordering variable specified - normal interpretation
        .vbse(ang=paste0("Durbin-Watson test with ordering by: ", orderDW),
              fr=paste0("Test de Durbin-Watson avec tri par : ", orderDW),
              k=counter, cpt="off", verbose=verbose)
        if (pvalt < dwtest_alpha) {
          .vbse(ang=paste0("Warning!\n\tDurbin-Watson test (dwtest()) - Bad independence of the residuals. p.value: ", .format_pval(pvalt)),
                fr=paste0("Attention !\n\tTest de Durbin-Watson (dwtest()) - Mauvaise indépendance des résidus. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
          dw_count_in_error <- TRUE
        } else {
          .vbse(ang=paste0("Durbin-Watson test (dwtest()) - Good independence of the residuals. p.value: ", .format_pval(pvalt)),
                fr=paste0("Test de Durbin-Watson (dwtest()) - Bonne indépendance des résidus. p.value : ", .format_pval(pvalt)),
                k=counter, cpt="off", verbose=verbose)
        }
      }
    } else if (inherits(reg, "lmerMod")) {
      pvalt <- min(dwtest_mix(reg))
      .vbse(ang="Warning: only by sorting the values as a function of Y and according to random effect variables.",
            fr="Attention : uniquement en triant les valeurs en fonction de Y et en fonction des variables d'effet aléatoire.",
            k=counter, cpt="off", verbose=verbose)
      if (pvalt < dwtest_alpha) {
        .vbse(ang=paste0("Warning!\n\tDurbin-Watson test (dwtest()) - Bad independence of the residuals. p.value: ", .format_pval(pvalt)),
              fr=paste0("Attention !\n\tTest de Durbin-Watson (dwtest()) - Mauvaise indépendance des résidus. p.value : ", .format_pval(pvalt)),
              k=counter, cpt="off", verbose=verbose)
        dw_count_in_error <- TRUE
      } else {
        .vbse(ang=paste0("Durbin-Watson test (dwtest()) - Good independence of the residuals. p.value: ", .format_pval(pvalt)),
              fr=paste0("Test de Durbin-Watson (dwtest()) - Bonne indépendance des résidus. p.value : ", .format_pval(pvalt)),
              k=counter, cpt="off", verbose=verbose)
      }
    }

    # Only count DW error if orderDW was specified or mixed model
    if (dw_count_in_error) {
      error <- FALSE
    }

    counter <- .vbse(ang="- Analysis of variance of residuals.", fr="- Analyse de la variance des résidus.", k=counter, cpt="on", verbose=verbose)
    if (length(coef(reg)) >= 2 & length(resid(reg)) <= 1000) {
      if (inherits(reg, "lm")) {
        bptest(reg)$p.value -> pvalt
      } else if (inherits(reg, "lmerMod")) {
			bptest_lmer(reg)$p.value -> pvalt
      }
      if (pvalt < bptest_alpha) {
        .vbse(ang=paste0("\tWarning!\n\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value: ", .format_pval(pvalt)), fr=paste0("Attention !\n\tTest de Breush-Pagan (bptest()) - Variance non-constante des résidus. p.value : ", .format_pval(pvalt)), k=counter, cpt="off", verbose=verbose)
		if (inherits(reg, "lmerMod")) {
			.vbse(ang="\tHeteroscedasticity is not always an issue in a mixed model. See {nlme}.", fr="\tL'hétéroscédasticité n'est pas toujours un problème sur un modèle mixte. Voyez {nlme}.", k=counter, cpt="off", verbose=verbose)
		}

        error <- FALSE
      } else {
        .vbse(ang=paste0("\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value: ", .format_pval(pvalt)), fr=paste0("Test de Breush-Pagan (bptest()) - Variance constante des résidus. p.value : ", .format_pval(pvalt)), k=counter, cpt="off", verbose=verbose)
      }
    } else if (length(resid(reg)) > 1000) {
      .vbse(ang="\tWarning!\n\tToo many values to justify checking the constancy of the variance. Breush-Pagan test bptest() of {lmtest} may be hypersensitive.", fr="\tAttention !\n\tTrop de valeurs pour justifier la vérification de la constance de la variance. Le test de Breush-Pagan bptest() de {lmtest} peut être hypersensible.", k=counter, cpt="off", verbose=verbose)
    }
    if (inherits(reg, "lm")) {
      cooks.distance(reg) -> cooksd
    } else if (inherits(reg, "lmerMod")) {
      cooks_distance_lmer(reg) -> cooksd
      unlist(cooksd) -> cooksd
    }
    counter <- .vbse(ang="- Analysis of leverage effect.", fr="- Analyse de l'effet de levier.", k=counter, cpt="on", verbose=verbose)
    if (max(cooksd, na.rm = TRUE) > 1) {
      .vbse(ang=paste0("Warning!\n\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance()): ", round(max(cooksd, na.rm = TRUE), 4)), fr=paste0("Attention !\n\tDistance de Cook (cooks.distance()) - Effet de levier. max(cooks.distance()) : ", round(max(cooksd, na.rm = TRUE), 4)), k=counter, cpt="off", verbose=verbose)
      error <- FALSE
    } else {
      .vbse(ang=paste0("Cook's distance (cooks.distance()) - No leverage effect. max(cooks.distance()): ", round(max(cooksd, na.rm = TRUE), 4)), fr=paste0("Distance de Cook (cooks.distance()) - Pas d'effet de levier. max(cooks.distance()) : ", round(max(cooksd, na.rm = TRUE), 4)), k=counter, cpt="off", verbose=verbose)
    }

    if (count_terms_excluding_intercept(reg) >= 2) {
      vif_reg <- suppressMessages(vif(reg))
      counter <- .vbse(ang="- Multicollinearity test (VIF).", fr="- Test de multicolinéarité (VIF).", k=counter, cpt="on", verbose=verbose)
      if (max(vif_reg, na.rm = TRUE) > 5) {
        .vbse(ang=paste0("Warning!\n\tThe variance inflation factor (VIF) indicates collinear variables with car::vif(). max(vif()): ", round(max(vif_reg, na.rm = TRUE), 2),
                         "\n\tRecommendations: Identify redundant variables (correlations > 0.8); Remove one or combine them (PCA, composite index); or use regularized regression (ridge, lasso, elastic net)."),
              fr=paste0("Attention !\n\tLe facteur d'inflation de la variance (VIF) indique des variables colinéaires avec car::vif(). max(vif()) : ", round(max(vif_reg, na.rm = TRUE), 2),
                        "\n\tRecommandations : Identifier les variables redondantes (corrélations > 0,8) ; Supprimer l'une d'elles ou les combiner (ACP, indice synthétique) ; ou utiliser une régression régularisée (ridge, lasso, elastic net)."),
              k=counter, cpt="off", verbose=verbose)
        error <- FALSE
      } else {
        .vbse(ang=paste0("No significant multicollinearity problem with car::vif(). max(vif()): ", round(max(vif_reg, na.rm = TRUE), 2)), fr=paste0("Pas de problème significatif de multicolinéarité avec car::vif(). max(vif()) : ", round(max(vif_reg, na.rm = TRUE), 2)), k=counter, cpt="off", verbose=verbose)
      }
    }

    if (boot == TRUE & inherits(reg, "lm")) {
      counter <- .vbse(ang="- Analysis of solidity of model by bootstrap.", fr="- Analyse de la solidité du modèle par bootstrap.", k=counter, cpt="on", verbose=verbose)
      if (length(data) > 0) {
        bootreg(reg, verbose = FALSE, alpha = alpha, conf.level = conf.level, plot = plot, data = data) -> bootres
      } else {
        bootreg(reg, verbose = FALSE, alpha = alpha, conf.level = conf.level, plot = plot) -> bootres
      }
      if (bootres == FALSE) {
        .vbse(ang="Warning!\n\tBootstrap (bootreg()) - Fragility of the model in bootstrap. Please, use bootreg()", fr="Attention !\n\tBootstrap (bootreg()) - Fragilité du modèle en bootstrap. Veuillez utiliser bootreg()", k=counter, cpt="off", verbose=verbose)
        error <- FALSE
      } else {
        .vbse(ang="Bootstrap (bootreg()) - Solidity of the model in bootstrap.", fr="Bootstrap (bootreg()) - Solidité du modèle en bootstrap.", k=counter, cpt="off", verbose=verbose)
      }
    }
  }

  # Return both validation result and updated counter for m.test() integration
  if (exit_k==TRUE) {
	  result <- list(
		valid = error,
		k = counter
	  )
	} else {result <- error}		
  return(result)
}

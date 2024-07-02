#' Function to validate a regression model
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
#'
#' @return This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).
#' @return valreg will therefore validate the regression model, control the p-values and, possibly (boot argument), control the reliability by bootstrap with the bootreg function.
#' @import lmtest
#' @import tseries
#' @importFrom stats cooks.distance
#' @importFrom lme4 lmer ranef
#' @importFrom car vif
#' @importFrom graphics plot
#' @importFrom MHTdiscrete Sidak.p.adjust
#' @importFrom MASS fitdistr
#' @export
#'
#' @examples
#' # Example 1: Linear model
#' data(iris)
#' reg <- lm(Sepal.Length~.,data=iris[,1:4])
#' valreg(reg,verbose=TRUE)
#' 
#' # Example 2: Mixed model
#' library(lme4)
#' data(sleepstudy)
#' reg_mixed <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' valreg(reg_mixed, verbose=TRUE)
valreg <- function(reg, verbose=TRUE, nvar=5, boot=TRUE, alpha=0.05, conf.level=0.95,
                   plot=FALSE, data=c(), raintest_alpha=0.05, dwtest_alpha=0.03,
                   shapiro_alpha=0.05, bptest_alpha=0.05) {
  error <- TRUE
  counter <- 1

  locale <- Sys.getlocale("LC_MESSAGES")
  lang <- ifelse(grepl("fr", locale), "fr", "en")

  msg <- function(en, fr) {
    if (lang == "fr") fr else en
  }

  lrt <- function(model) {
    formula <- formula(model)
    response <- as.character(formula[[2]])
    terms <- attr(terms(formula), "term.labels")
    random_effects <- grep("\\|", terms, value = TRUE)
    if (length(random_effects) == 0) {
      stop(msg("The provided model has no random effects.", "Le modèle fourni n'a pas d'effets aléatoires."))
    }
    random_effects_str <- paste("(1 |", gsub(".*\\|", "", random_effects), ")", sep = "")
    formula_nul <- as.formula(paste(response, "~", paste(random_effects_str, collapse = " + ")))
    data_used <- model.frame(model)
    model_complet <- lmer(formula, data = data_used)
    model_nul <- lmer(formula_nul, data = data_used)
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
      stop(msg("The provided model has no random effects.", "Le modèle fourni n'a pas d'effets aléatoires."))
    }
    fixed_effects <- setdiff(terms, random_effects)
    if (length(fixed_effects) == 0) {
      stop(msg("The provided model has no fixed effects.", "Le modèle fourni n'a pas d'effets fixes."))
    }
    fixed_effects_str <- paste(fixed_effects, collapse = " + ")
    formula_lm <- as.formula(paste(response, "~", fixed_effects_str))
    data_used <- model.frame(model)
    model_complet <- lmer(formula, data = data_used)
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
        warning(msg(paste("Variable", var, "does not exist in the data."),
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
        warning(msg("No valid random effect variables found for crossed sorting.",
                    "Aucune variable d'effet aléatoire valide trouvée pour le tri croisé."))
      }
    }
    return(unlist(results))
  }

  bptest_lmer <- function(model) {
    if (!inherits(model, "lmerMod")) {
      stop(msg("The provided model is not an lmer model.", "Le modèle fourni n'est pas un modèle lmer."))
    }
    formula_fixed <- reformulate(attr(terms(model), "term.labels"), response = as.character(formula(model)[[2]]))
    model_lm <- lm(formula_fixed, data = model@frame)
    bp_test <- bptest(model_lm)
    return(bp_test)
  }

  count_terms_excluding_intercept <- function(model) {
    if (!inherits(model, c("lmerMod", "lm"))) {
      stop(msg("The model must be of type lmerMod (fitted with lme4) or lm.",
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
    if (any(is.na(model$coefficients))) {
      stop(msg("NA in coefficients.", "NA dans les coefficients."))
    }
    data <- model$model
    terms <- attr(model$terms, "term.labels")
    form <- formula(model)
    if (length(terms) > 1) {
      formula <- paste0(deparse(form), "+", paste(sapply(terms, function(x) paste0("I(", x, "^2)")), collapse = "+"), "+", paste(combn(terms, 2, function(x) paste0("I(", x[1], "*", x[2], ")")), collapse = "+"))
    } else {
      formula <- paste0(deparse(form), "+", paste0("I(", terms, "^2)"))
    }
    model <- lm(as.formula(formula), data = data)
    return(bptest(model))
  }

  nvar1 <- length(coef(reg))
  if (nvar1 > nvar) {
    if (verbose) cat(msg("More than", "Plus de"), nvar, msg("variables\n", "variables\n"))
    error <- FALSE
  } else if (length(fitted(reg)) < (nvar + 1)) {
    if (verbose) cat(msg("Not enough values in the subset\n", "Pas assez de valeurs dans le sous-ensemble\n"))
    error <- FALSE
  } else if (inherits(reg, "lm") && (length(summary(reg)$coefficients[, "Pr(>|t|)"])) != (length(coef(reg)))) {
    error <- FALSE
  } else {
    if (verbose) cat(counter, msg("- Analysis of the p-values of the model and its coefficients.\n", "- Analyse des p-values du modèle et de ses coefficients.\n"))
    if (inherits(reg, "lm")) {
      pval_mdl <- pf(summary(reg)$fstatistic[1], summary(reg)$fstatistic[2], summary(reg)$fstatistic[3], lower.tail = FALSE)
      if (pval_mdl > alpha) {
        if (verbose) cat(msg("\tWarning!\n\tBad significance of the model. p-value:", "\tAttention !\n\tMauvaise signification du modèle. p-value :"), pval_mdl, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tGood significance of the model. p-value:", "\tBonne signification du modèle. p-value :"), pval_mdl, "\n")
      }
    } else if (inherits(reg, "lmerMod")) {
      if (verbose) cat(msg("Please note, mixed model. Calculation of a global p-value for comparison of the model with random effects alone by Likelihood Ratio Test (LRT).\n",
                           "Veuillez noter, modèle mixte. Calcul d'une p-value globale pour la comparaison du modèle avec effets aléatoires seuls par test du rapport de vraisemblance (LRT).\n"))
      pval_mdl <- lrt(reg)
      if (pval_mdl > alpha) {
        if (verbose) cat(msg("\tWarning!\n\tThe fixed effect adds nothing to the random effect. p-value:", "\tAttention !\n\tL'effet fixe n'ajoute rien à l'effet aléatoire. p-value :"), pval_mdl, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tThe fixed effect combined with the random effect is more relevant than the random effect alone. p-value:", "\tL'effet fixe combiné à l'effet aléatoire est plus pertinent que l'effet aléatoire seul. p-value :"), pval_mdl, "\n")
      }
      aic_mdl <- lrt2(reg)
      if (aic_mdl > 0) {
        if (verbose) cat(msg("\tWarning!\n\tThe random effect does not improve the model compared to fixed effects alone.\n", "\tAttention !\n\tL'effet aléatoire n'améliore pas le modèle par rapport aux effets fixes seuls.\n"))
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tThe random effect improves the model compared to fixed effects alone.\n", "\tL'effet aléatoire améliore le modèle par rapport aux effets fixes seuls.\n"))
        pval_of_ranova <- max(suppressMessages(suppressWarnings({ranova(model)}))$`Pr(>Chisq)`, na.rm = TRUE)
        if (pval_of_ranova < alpha) {
          if (verbose) cat(msg("\t\tConfirmation by ranova().\n", "\t\tConfirmation par ranova().\n"))
        } else {
          if (verbose) cat(msg("\t\tNon-confirmation by ranova().\n", "\t\tNon-confirmation par ranova().\n"))
          error <- FALSE
        }
      }
    }

    pval_coeff <- summary(reg)$coefficients[, "Pr(>|t|)"]
    if (max(pval_coeff) > alpha) {
      if (verbose) cat(msg("\tWarning!\n\tBad significance of the coefficients. max(p.value):", "\tAttention !\n\tMauvaise signification des coefficients. max(p.value) :"), max(pval_coeff), "\n")
      error <- FALSE
    } else {
      if (verbose) cat(msg("\tGood significance of the coefficients. max(pval_coeff):", "\tBonne signification des coefficients. max(pval_coeff) :"), max(pval_coeff), "\n")
    }

    counter <- counter + 1
    if (verbose) cat(counter, msg("- Analysis of distribution of residuals.\n", "- Analyse de la distribution des résidus.\n"))
    if (length(resid(reg)) <= 50) {
      shapiro.test(resid(reg))$p.value -> pvalt
      if (pvalt < shapiro_alpha) {
        if (verbose) cat(msg("\tWarning!\n\tShapiro-Wilk test (shapiro.test()) - Non-normal distribution of residuals. p.value:", "\tAttention !\n\tTest de Shapiro-Wilk (shapiro.test()) - Distribution non normale des résidus. p.value :"), pvalt, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tShapiro-Wilk test (shapiro.test()) - Normal distribution of residuals. p.value:", "\tTest de Shapiro-Wilk (shapiro.test()) - Distribution normale des résidus. p.value :"), pvalt, "\n")
      }
    } else if (length(resid(reg)) <= 1000) {
      cat(msg("Warning! No Shapiro-Wilk test because more than 50 values.\n", "Attention ! Pas de test de Shapiro-Wilk car plus de 50 valeurs.\n"))
      jarque.bera.test(residuals(reg))$p.value -> pvalt
      if (pvalt < shapiro_alpha) {
        if (verbose) cat(msg("\tWarning!\n\tJarque-Bera test (jarque.bera.test() of {tseries}) - Non-normal distribution of residuals. p.value:", "\tAttention !\n\tTest de Jarque-Bera (jarque.bera.test() de {tseries}) - Distribution non normale des résidus. p.value :"), pvalt, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tJarque-Bera test (jarque.bera.test() of {tseries}) - Normal distribution of residuals. p.value:", "\tTest de Jarque-Bera (jarque.bera.test() de {tseries}) - Distribution normale des résidus. p.value :"), pvalt, "\n")
      }
    } else {
      cat(msg("Warning! No Shapiro-Wilk test or Jarque-Bera test because more than 1000 values.\n", "Attention ! Pas de test de Shapiro-Wilk ou de Jarque-Bera car plus de 1000 valeurs.\n"))
      ks.test(residuals(reg), "pnorm", mean = mean(residuals(reg)), sd = sd(residuals(reg)))$p.value -> pvalt
      if (pvalt < shapiro_alpha) {
        if (verbose) cat(msg("\tWarning!\n\tKolmogorov-Smirnov test (ks.test()) - Non-normal distribution of residuals. p.value:", "\tAttention !\n\tTest de Kolmogorov-Smirnov (ks.test()) - Distribution non normale des résidus. p.value :"), pvalt, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tKolmogorov-Smirnov test (ks.test()) - Normal distribution of residuals. p.value:", "\tTest de Kolmogorov-Smirnov (ks.test()) - Distribution normale des résidus. p.value :"), pvalt, "\n")
      }
    }
    if (inherits(reg, "lmerMod")) {
      pvalt <- norm_lmer(reg)
      if (pvalt < shapiro_alpha) {
        if (verbose) cat(msg("\tWarning!\n\tNon-normality of random effects for at least one random variable. p.value:", "\tAttention !\n\tNon-normalité des effets aléatoires pour au moins une variable aléatoire. p.value :"), pvalt, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tNormality of random effects for at least one random variable. p.value:", "\tNormalité des effets aléatoires pour au moins une variable aléatoire. p.value :"), pvalt, "\n")
      }
    }

    if (inherits(reg, "lm")) {
      counter <- counter + 1
      if (verbose) cat(counter, msg("- Analysis of the adequacy of model (Equivalence between the global model and the model established on the best points).\n", "- Analyse de l'adéquation du modèle (équivalence entre le modèle global et le modèle établi sur les meilleurs points).\n"))
      if ("(Intercept)" %in% names(coef(reg))) {
        nvar2 <- length(coef(reg)) - 1
      } else {
        nvar2 <- length(coef(reg))
      }
      if (nvar2 > 1) {
        raintest(reg, order.by = "mahalanobis")$p.value -> pvalt
        if (pvalt < raintest_alpha) {
          if (verbose) cat(msg("\tWarning!\n\tRainbow test ordered by mahalanobis (raintest()) - Bad adequacy. p.value:", "\tAttention !\n\tTest de Rainbow ordonné par mahalanobis (raintest()) - Mauvaise adéquation. p.value :"), pvalt, "\n")
          error <- FALSE
        } else {
          if (verbose) cat(msg("\tRainbow test (raintest()) - Good adequacy. p.value:", "\tTest de Rainbow (raintest()) - Bonne adéquation. p.value :"), pvalt, "\n")
        }
      } else {
        if (length(data) == 0) {
          data <- reg$model
        }
        z <- as.formula(paste0("~", attr(reg$terms, "term.labels")))
        raintest(reg, order.by = z, data = data)$p.value -> pvalt
        if (pvalt < raintest_alpha) {
          if (verbose) cat(msg("\tWarning!\n\tRainbow test ordered by X (raintest()) - Bad adequacy. p.value:", "\tAttention !\n\tTest de Rainbow ordonné par X (raintest()) - Mauvaise adéquation. p.value :"), pvalt, "\n")
          error <- FALSE
        } else {
          if (verbose) cat(msg("\tRainbow test (raintest()) - Good adequacy. p.value:", "\tTest de Rainbow (raintest()) - Bonne adéquation. p.value :"), pvalt, "\n")
        }
      }
    }

    counter <- counter + 1
    if (verbose) cat(counter, msg("- Analysis of independence of the residuals\n", "- Analyse de l'indépendance des résidus\n"))
    if (inherits(reg, "lm")) {
      if (length(data) > 0) {
        ypred <- predict(reg, data)
      } else {
        ypred <- predict(reg, reg$model)
      }
      try(dwtest(reg, ypred)$p.value) -> pvalt
      if (verbose) cat(msg("Warning: only by sorting the values as a function of Y. Check manually for other criteria.\n",
                           "Attention : uniquement en triant les valeurs en fonction de Y. Vérifiez manuellement pour d'autres critères.\n"))
    } else if (inherits(reg, "lmerMod")) {
      pvalt <- min(dwtest_mix(reg))
      if (verbose) cat(msg("Warning: only by sorting the values as a function of Y and according to random effect variables.\n",
                           "Attention : uniquement en triant les valeurs en fonction de Y et en fonction des variables d'effet aléatoire.\n"))
    }
    if (pvalt < dwtest_alpha) {
      if (verbose) cat(msg("\tWarning!\n\tDurbin-Watson test (dwtest()) - Bad independence of the residuals. p.value:", "\tAttention !\n\tTest de Durbin-Watson (dwtest()) - Mauvaise indépendance des résidus. p.value :"), pvalt, "\n")
      error <- FALSE
    } else {
      if (verbose) cat(msg("\tDurbin-Watson test (dwtest()) - Good independence of the residuals. p.value:", "\tTest de Durbin-Watson (dwtest()) - Bonne indépendance des résidus. p.value :"), pvalt, "\n")
    }

    counter <- counter + 1
    if (verbose) cat(counter, msg("- Analysis of variance of residuals.\n", "- Analyse de la variance des résidus.\n"))
    if (length(coef(reg)) >= 2 & length(resid(reg)) <= 1000) {
      if (inherits(reg, "lm")) {
        bptest(reg)$p.value -> pvalt
      } else if (inherits(reg, "lmerMod")) {
        bptest_lmer(reg)$p.value -> pvalt
      }
      if (pvalt < bptest_alpha) {
        if (verbose) cat(msg("\tWarning!\n\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value:", "\tAttention !\n\tTest de Breush-Pagan (bptest()) - Variance non-constante des résidus. p.value :"), pvalt, "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value:", "\tTest de Breush-Pagan (bptest()) - Variance constante des résidus. p.value :"), pvalt, "\n")
      }
    } else if (length(resid(reg)) > 1000) {
      cat(msg("\tWarning!\n\tToo many values to justify checking the constancy of the variance. Breush-Pagan test bptest() of {lmtest} may be hypersensitive.\n",
              "\tAttention !\n\tTrop de valeurs pour justifier la vérification de la constance de la variance. Le test de Breush-Pagan bptest() de {lmtest} peut être hypersensible.\n"))
    }

    counter <- counter + 1
    if (inherits(reg, "lm")) {
      cooks.distance(reg) -> cooksd
    } else if (inherits(reg, "lmerMod")) {
      cooks.distance_lmer(reg) -> cooksd
      unlist(cooksd) -> cooksd
    }
    if (verbose) cat(counter, msg("- Analysis of leverage effect.\n", "- Analyse de l'effet de levier.\n"))
    if (max(cooksd, na.rm = TRUE) > 1) {
      if (verbose) cat(msg("\tWarning!\n\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance()):", "\tAttention !\n\tDistance de Cook (cooks.distance()) - Effet de levier. max(cooks.distance()) :"), max(cooksd, na.rm = TRUE), "\n")
      error <- FALSE
    } else {
      if (verbose) cat(msg("\tCook's distance (cooks.distance()) - No leverage effect. max(cooks.distance()):", "\tDistance de Cook (cooks.distance()) - Pas d'effet de levier. max(cooks.distance()) :"), max(cooksd, na.rm = TRUE), "\n")
    }

    if (count_terms_excluding_intercept(reg) >= 2) {
      counter <- counter + 1
      car::vif(reg) -> vif_reg
      if (verbose) cat(counter, msg("- Multicollinearity test (VIF).\n", "- Test de multicolinéarité (VIF).\n"))
      if (max(vif_reg, na.rm = TRUE) > 5) {
        if (verbose) cat(msg("\tWarning!\n\tThe variance inflation factor (VIF) indicates collinear variables with car::vif(). max(vif()):", "\tAttention !\n\tLe facteur d'inflation de la variance (VIF) indique des variables collinéaires avec car::vif(). max(vif()) :"), max(vif_reg, na.rm = TRUE), "\n")
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tNo significant multicollinearity problem with car::vif(). max(vif()):", "\tPas de problème significatif de multicolinéarité avec car::vif(). max(vif()) :"), max(vif_reg, na.rm = TRUE), "\n")
      }
    }

    if (boot == TRUE & inherits(reg, "lm")) {
      counter <- counter + 1
      if (verbose) cat(counter, msg("- Analysis of solidity of model by bootstrap.\n", "- Analyse de la solidité du modèle par bootstrap.\n"))
      if (length(data) > 0) {
        bootreg(reg, verbose = FALSE, alpha = alpha, conf.level = conf.level, plot = plot, data = data) -> bootres
      } else {
        bootreg(reg, verbose = FALSE, alpha = alpha, conf.level = conf.level, plot = plot) -> bootres
      }
      if (bootres == FALSE) {
        if (verbose) cat(msg("\tWarning!\n\tBootstrap (bootreg()) - Fragility of the model in bootstrap. Please, use bootreg()\n", "\tAttention !\n\tBootstrap (bootreg()) - Fragilité du modèle en bootstrap. Veuillez utiliser bootreg()\n"))
        error <- FALSE
      } else {
        if (verbose) cat(msg("\tBootstrap (bootreg()) - Solidity of the model in bootstrap.\n", "\tBootstrap (bootreg()) - Solidité du modèle en bootstrap.\n"))
      }
    }
  }
  return(error)
}

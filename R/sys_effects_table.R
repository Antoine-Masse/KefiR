#' Convert an ANOVA-like table to the $effects table of m.test()
#'
#' @param tab Table with one row per term and a p-value column.
#' @return data.frame with columns term and p.value, or NULL if not applicable.
#' @keywords internal
.effects_from_table <- function(tab) {
  if (is.null(tab)) return(NULL)
  tab <- tryCatch(as.data.frame(tab), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  pval_col <- intersect(c("Pr(>F)", "p.value", "Pr(Prob)", "P(>|F|)", "Pr(>Chisq)"), colnames(tab))
  if (length(pval_col) == 0) return(NULL)
  terms_names <- if ("model term" %in% colnames(tab)) tab[["model term"]] else rownames(tab)
  terms_names <- trimws(as.character(terms_names))
  keep <- !grepl("^Residuals|^\\(Intercept\\)$", terms_names, ignore.case = TRUE) &
    !is.na(tab[[pval_col[1]]])
  if (!any(keep)) return(NULL)
  data.frame(term = terms_names[keep], p.value = as.numeric(tab[[pval_col[1]]][keep]),
             stringsAsFactors = FALSE)
}

#' P-value of the first factor of a formula in an $effects table
#'
#' Unique rule of m.test() (bp.log, BP-031): the global p-value is the test of the
#' first factor of the formula, which is also the first factor of the post-hocs.
#'
#' @param effects Output of .effects_from_table().
#' @param formula Model formula.
#' @return Numeric p-value, or NA.
#' @keywords internal
.first_factor_pvalue <- function(effects, formula) {
  if (is.null(effects) || nrow(effects) == 0) return(NA_real_)
  labels <- tryCatch(attr(terms(formula), "term.labels"), error = function(e) character(0))
  labels <- labels[!grepl(":", labels, fixed = TRUE) & !grepl("^Error\\(", labels)]
  found <- match(labels, effects$term)
  found <- found[!is.na(found)]
  if (length(found) > 0) effects$p.value[found[1]] else effects$p.value[1]
}

#' $effects table of a multi-factor analysis, from the test actually conducted
#'
#' @param model Fitted parametric model, or NULL.
#' @param robust_results Robust analysis results, or NULL.
#' @param formula Model formula.
#' @param data Analysis data.
#' @param alpha Significance level.
#' @return data.frame with columns term and p.value, or NULL.
#' @keywords internal
.multi_effects_table <- function(model, robust_results, formula, data, alpha = 0.05) {
  test_result <- robust_results$test_result
  method <- robust_results$method

  if (!is.null(test_result)) {
    # Modèle mixte
    if (inherits(test_result, c("lmerMod", "lmerModLmerTest"))) {
      return(.effects_from_table(tryCatch(stats::anova(test_result, type = "III"),
                                          error = function(e) NULL)))
    }
    # WRS2::t2way / t3way : une p-value par terme, nommée A, B, AB...
    if (isTRUE(method %in% c("WRS2_t2way", "WRS2_t3way"))) {
      p_names <- grep("\\.p\\.value$", names(test_result), value = TRUE)
      if (length(p_names) == 0) return(NULL)
      return(data.frame(term = sub("\\.p\\.value$", "", p_names),
                        p.value = as.numeric(unlist(test_result[p_names])),
                        stringsAsFactors = FALSE))
    }
    # Test unique (Kruskal-Wallis, Friedman)
    if (inherits(test_result, "htest")) {
      return(data.frame(term = "global", p.value = test_result$p.value, stringsAsFactors = FALSE))
    }
    # Permutations : summary() renvoie une liste de tables
    if (inherits(test_result, c("aovp", "lmp"))) {
      smry <- tryCatch(summary(test_result), error = function(e) NULL)
      if (is.list(smry) && !is.data.frame(smry) && length(smry) > 0) {
        return(.effects_from_table(smry[[1]]))
      }
      return(.effects_from_table(smry))
    }
    # Scheirer-Ray-Hare et tables similaires
    return(.effects_from_table(test_result))
  }

  if (is.null(model)) return(NULL)

  # Mesures répétées : une table par strate
  if (inherits(model, "aovlist")) {
    smry <- tryCatch(summary(model), error = function(e) NULL)
    if (is.null(smry)) return(NULL)
    tables <- lapply(smry, function(s) if (length(s) > 0) .effects_from_table(s[[1]]))
    tables <- Filter(Negate(is.null), tables)
    if (length(tables) == 0) return(NULL)
    res <- do.call(rbind, tables)
    rownames(res) <- NULL
    return(res)
  }

  if (!inherits(model, c("aov", "lm"))) return(NULL)

  # Modèle paramétrique : même type de sommes des carrés que le tableau affiché
  ss_type <- tryCatch(
    .select_ss_type(formula = formula, data = data, model = model, analysis_type = "ANOVA",
                    alpha = alpha, verbose = FALSE)$ss_type,
    error = function(e) "II"
  )
  .effects_from_table(tryCatch(
    if (identical(ss_type, "I")) anova(model) else suppressMessages(car::Anova(model, type = ss_type)),
    error = function(e) tryCatch(anova(model), error = function(e2) NULL)
  ))
}

#' Compact letters from a lower-triangular matrix of p-values
#'
#' Same lettering as the other post-hocs of m.test(), via catego().
#'
#' @param p_lower Lower-triangular p-value matrix, in the format of pairwise.t.test().
#' @param position Named statistic used to sort the groups.
#' @param alpha Significance level.
#' @return data.frame with columns categories and letters, in the order of position.
#' @keywords internal
.letters_from_pairs <- function(p_lower, position, alpha = 0.05) {
  levels_all <- names(position)
  p_full <- matrix(NA_real_, length(levels_all), length(levels_all),
                   dimnames = list(levels_all, levels_all))
  for (r in rownames(p_lower)) {
    for (cc in colnames(p_lower)) {
      value <- p_lower[r, cc]
      if (!is.na(value)) p_full[r, cc] <- p_full[cc, r] <- value
    }
  }
  ord <- order(position)
  p_sorted <- p_full[ord, ord, drop = FALSE]
  p_sorted[upper.tri(p_sorted, diag = TRUE)] <- NA
  p_sorted <- p_sorted[-1, -length(levels_all), drop = FALSE]
  letters_sorted <- catego(list(p.value = p_sorted), alpha = alpha)$groups
  data.frame(categories = levels_all,
             letters = letters_sorted[match(levels_all, letters_sorted[, "categories"]), "groups"],
             stringsAsFactors = FALSE)
}

#' Multivariate statistic and F approximation from eigenvalues
#'
#' Same formulas as summary.manova(), whose helpers are not exported by {stats}.
#' Each formula is checked against summary.manova() in tests/testthat/test-manova.R.
#'
#' @param eig Eigenvalues of solve(E) %*% H.
#' @param q Hypothesis degrees of freedom.
#' @param df_res Error degrees of freedom.
#' @param test "Wilks", "Pillai", "Hotelling-Lawley" or "Roy".
#' @return Numeric vector: statistic, approximate F, num df, den df, p-value.
#' @keywords internal
.manova_stat <- function(eig, q, df_res, test) {
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df_res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  res <- switch(test,
    "Wilks" = {
      stat <- prod(1 / (1 + eig))
      t1 <- df_res - 0.5 * (p - q + 1)
      t2 <- (p * q - 2) / 4
      t3 <- p^2 + q^2 - 5
      t3 <- if (t3 > 0) sqrt(((p * q)^2 - 4) / t3) else 1
      c(stat, ((stat^(-1 / t3) - 1) * (t1 * t3 - 2 * t2)) / p / q, p * q, t1 * t3 - 2 * t2)
    },
    "Pillai" = {
      stat <- sum(eig / (1 + eig))
      t1 <- 2 * m + s + 1
      t2 <- 2 * n + s + 1
      c(stat, (t2 / t1 * stat) / (s - stat), s * t1, s * t2)
    },
    "Hotelling-Lawley" = {
      stat <- sum(eig)
      t1 <- 2 * m + s + 1
      t2 <- 2 * (s * n + 1)
      c(stat, (t2 * stat) / s / s / t1, s * t1, t2)
    },
    "Roy" = {
      stat <- max(eig)
      t1 <- max(p, q)
      t2 <- df_res - t1 + q
      c(stat, (t2 * stat) / t1, t1, t2)
    })
  c(res, pf(res[2], res[3], res[4], lower.tail = FALSE))
}

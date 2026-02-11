.formulator_safe <- function(formula, data = NULL, mode = c("simple", "linear", "lmer"), debug=FALSE) {
  mode <- match.arg(mode)
  stopifnot(inherits(formula, "formula"))

  f <- formula

  # Auto-detection syntaxe lmer : si (1|id) ou (var|id) present, forcer mode "lmer"
  formula_str <- deparse(formula)
  if (grepl("[(][^)]+[|][^)]+[)]", formula_str) && mode != "lmer") {
    .dbg(".formulator_safe() - lmer syntax detected, switching to mode='lmer'",
         ".formulator_safe() - Syntaxe lmer d\u00e9tect\u00e9e, passage en mode='lmer'",
         debug = debug)
    mode <- "lmer"
  }

  # 0) si data= fourni, normaliser data$col -> col
  if (!is.null(data)) {
    f <- .normalize_formula_dollar(f, data)
  }

  # 1) Ne JAMAIS alterer la LHS
  lhs <- f[[2L]]

  # 2) Essayer d'utiliser .formulator(), mais on RESTAURE la LHS quoi qu'il arrive
  f2 <- try(.formulator(f, mode = mode, debug=debug), silent = TRUE)
  if (!inherits(f2, "try-error") && inherits(f2, "formula")) {
    f2[[2L]] <- lhs
    f <- f2
  }

  # 3) Si data= fourni, re-normaliser (au cas ou .formulator a remis des $)
  if (!is.null(data)) {
    f <- .normalize_formula_dollar(f, data)
  }

  f
}

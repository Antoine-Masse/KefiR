.formulator_safe <- function(formula, data = NULL, mode = c("simple", "linear", "lmer"), debug=FALSE) {
  mode <- match.arg(mode)
  stopifnot(inherits(formula, "formula"))

  f <- formula

  # Auto-détection syntaxe lmer : si (1|id) ou (var|id) présent, forcer mode "lmer"
  formula_str <- deparse(formula)
  if (grepl("[(][^)]+[|][^)]+[)]", formula_str) && mode != "lmer") {
    .dbg(".formulator_safe() - lmer syntax detected, switching to mode='lmer'",
         ".formulator_safe() - Syntaxe lmer détectée, passage en mode='lmer'",
         debug = debug)
    mode <- "lmer"
  }

  # 0) si data= fourni, normaliser data$col -> col
  if (!is.null(data)) {
    f <- .normalize_formula_dollar(f, data)
  }

  # 1) Ne JAMAIS altérer la LHS
  lhs <- f[[2L]]

  # 2) Essayer d'utiliser .formulator(), mais on RESTAURE la LHS quoi qu’il arrive
  f2 <- try(.formulator(f, mode = mode, debug=debug), silent = TRUE)
  if (!inherits(f2, "try-error") && inherits(f2, "formula")) {
    f2[[2L]] <- lhs
    f <- f2
  }

  # 3) Si data= fourni, re-normaliser (au cas où .formulator a remis des $)
  if (!is.null(data)) {
    f <- .normalize_formula_dollar(f, data)
  }

  f
}

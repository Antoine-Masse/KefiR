.normalize_formula_dollar <- function(f, data) {
  if (!inherits(f, "formula") || is.null(data)) return(f)

  data_sym <- as.name(deparse(substitute(data)))

  replace_dollar <- function(expr) {
    # 1) si ce n'est pas un appel/langage -> on rend tel quel
    if (!is.language(expr)) return(expr)
    # 2) si c'est un symbole nu (ex: B) -> tel quel (pas de expr[[i]])
    if (is.name(expr)) return(expr)

    # 3) cas data$col -> col (uniquement si LHS == data)
    if (is.call(expr) && identical(expr[[1L]], as.name("$"))) {
      if (identical(expr[[2L]], data_sym)) return(expr[[3L]])
      return(expr)  # ne pas toucher aux autres $
    }

    # 4) récursion sûre
    for (i in seq_along(expr)) expr[[i]] <- replace_dollar(expr[[i]])
    expr
  }

  f2 <- f
  f2[[2L]] <- replace_dollar(f[[2L]])  # LHS
  f2[[3L]] <- replace_dollar(f[[3L]])  # RHS
  f2
}

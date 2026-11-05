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

# Remplace factor(F) et as.factor(F) par F dans le membre droit d'une formule, et
# convertit la colonne F de data en facteur. Les appels avec d'autres arguments
# (ex. factor(F, levels = ...)) ou sur une expression sont laiss\u00e9s tels quels.
.strip_factor_calls <- function(f, data) {
  to_factor <- character(0)

  replace_factor <- function(expr) {
    if (!is.call(expr)) return(expr)
    is_factor_call <- identical(expr[[1L]], as.name("factor")) ||
      identical(expr[[1L]], as.name("as.factor"))
    if (is_factor_call && length(expr) == 2L && is.name(expr[[2L]]) &&
        as.character(expr[[2L]]) %in% names(data)) {
      to_factor <<- c(to_factor, as.character(expr[[2L]]))
      return(expr[[2L]])
    }
    for (i in seq_along(expr)[-1L]) expr[[i]] <- replace_factor(expr[[i]])
    expr
  }

  # z$b, o\u00f9 z contient la m\u00eame colonne b que data : \u00e9crire b, utilisable dans toute formule
  # et quel que soit l'environnement d'appel (dans une fonction par exemple)
  formula_env <- environment(f)
  replace_dollar <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1L]], as.name("$")) && is.name(expr[[2L]]) && is.name(expr[[3L]])) {
      col <- as.character(expr[[3L]])
      source_name <- as.character(expr[[2L]])
      if (col %in% names(data) && exists(source_name, envir = formula_env)) {
        source_obj <- get(source_name, envir = formula_env)
        if (is.data.frame(source_obj) && identical(source_obj[[col]], data[[col]])) return(expr[[3L]])
      }
      return(expr)
    }
    for (i in seq_along(expr)[-1L]) expr[[i]] <- replace_dollar(expr[[i]])
    expr
  }

  for (side in seq_along(f)[-1L]) f[[side]] <- replace_dollar(f[[side]])

  rhs <- length(f)
  f[[rhs]] <- replace_factor(f[[rhs]])
  for (var in unique(to_factor)) data[[var]] <- factor(data[[var]])
  list(formula = f, data = data)
}

.strip_data_dollar_safe <- function(f, data) {
  if (!inherits(f, "formula") || is.null(data)) return(f)
  data_sym <- as.name(deparse(substitute(data)))

  replace_dollar <- function(expr) {
    if (!is.language(expr)) return(expr)
    if (is.name(expr)) return(expr)  # <- évite expr[[i]] sur un symbole

    if (is.call(expr) && identical(expr[[1L]], as.name("$"))) {
      # On ne “déneste” que si le LHS du $ est exactement le symbole 'data'
      if (identical(expr[[2L]], data_sym)) return(expr[[3L]])
      return(expr)
    }
    for (i in seq_along(expr)) expr[[i]] <- replace_dollar(expr[[i]])
    expr
  }

  f2 <- f
  f2[[2L]] <- replace_dollar(f[[2L]])  # LHS
  f2[[3L]] <- replace_dollar(f[[3L]])  # RHS
  f2
}

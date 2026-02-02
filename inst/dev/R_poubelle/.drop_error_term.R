#.drop_error_term <- function(form) {
  #f <- deparse(form)
  #f <- gsub("\\s*\\+\\s*Error\\s*\\([^)]*\\)", "", f)
  #stats::as.formula(f)
.drop_error_term <- function(f) {
    stopifnot(inherits(f, "formula"))
    # 1) Essai par terms() avec special Error
    tt <- try(terms(f, specials = "Error"), silent = TRUE)
    if (!inherits(tt, "try-error")) {
      # Supprimer la response et reconstruire
      lhs <- as.character(f[[2L]])
      rhs_terms <- attr(delete.response(tt), "term.labels")
      if (length(rhs_terms) == 0L) rhs_terms <- "1"
      return(reformulate(rhs_terms, lhs = lhs))
    }
    # 2) Filet de sécurité : regex qui supprime tout bloc + Error(...)
    df <- deparse(f)
    df1 <- paste(df, collapse = "")
    # supprime " + Error(...)" n'importe où dans le RHS
    df1 <- sub("\\s*\\+\\s*Error\\s*\\([^)]*\\)\\s*$", "", df1)
    df1 <- sub("\\s*\\+\\s*Error\\s*\\([^)]*\\)", "", df1)
    as.formula(df1, env = environment(f))
}

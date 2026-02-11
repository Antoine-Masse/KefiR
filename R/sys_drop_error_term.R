#.drop_error_term <- function(form) {
  #f <- deparse(form)
  #f <- gsub("\\s*\\+\\s*Error\\s*\\([^)]*\\)", "", f)
  #stats::as.formula(f)
.drop_error_term <- function(f) {
    stopifnot(inherits(f, "formula"))

    # NOUVEAU: Détecter et gérer la syntaxe lmer (1|id) ou (var|id)
    # Ces termes doivent être SUPPRIMÉS pour model.frame() car ce n'est pas Error()
    formula_str <- deparse(f)
    is_lmer_syntax <- grepl("[(][^)]+[|][^)]+[)]", formula_str)

    if (is_lmer_syntax) {
      # Supprimer les termes lmer (1|id) ou (var|id) de la formule
      # Ils ne sont pas compris par model.frame()
      lhs <- deparse(f[[2L]])

      # Regex pour supprimer (1 | id) ou (var | id) avec espaces optionnels
      rhs_clean <- gsub("\\s*\\+\\s*[(][^)]+[|][^)]+[)]", "", formula_str)
      rhs_clean <- gsub("[(][^)]+[|][^)]+[)]\\s*\\+\\s*", "", rhs_clean)
      rhs_clean <- gsub("[(][^)]+[|][^)]+[)]", "", rhs_clean)

      # Extraire juste le RHS (après ~)
      rhs_part <- sub("^[^~]+~\\s*", "", rhs_clean)
      rhs_part <- trimws(rhs_part)

      # Si RHS vide, mettre 1 (intercept only)
      if (rhs_part == "" || rhs_part == lhs) rhs_part <- "1"

      return(as.formula(paste(lhs, "~", rhs_part), env = environment(f)))
    }

    # Gestion classique pour Error()
    # 1) Essai par terms() avec special Error
    tt <- try(terms(f, specials = "Error"), silent = TRUE)
    if (!inherits(tt, "try-error")) {
      # Supprimer la response et reconstruire
      lhs <- deparse(f[[2L]])
      rhs_terms <- attr(delete.response(tt), "term.labels")

      # IMPORTANT: Filtrer les termes Error(...)
      # terms() avec specials="Error" inclut encore les termes Error dans term.labels
      # Il faut les supprimer manuellement
      rhs_terms <- rhs_terms[!grepl("^Error\\s*\\(", rhs_terms)]

      if (length(rhs_terms) == 0L) rhs_terms <- "1"
      return(reformulate(rhs_terms, response = lhs))
    }
    # 2) Filet de sécurité : regex qui supprime tout bloc + Error(...)
    df <- deparse(f)
    df1 <- paste(df, collapse = "")
    # supprime " + Error(...)" n'importe où dans le RHS
    df1 <- sub("\\s*\\+\\s*Error\\s*\\([^)]*\\)\\s*$", "", df1)
    df1 <- sub("\\s*\\+\\s*Error\\s*\\([^)]*\\)", "", df1)
    as.formula(df1, env = environment(f))
}

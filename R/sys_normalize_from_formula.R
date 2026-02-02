.normalize_from_formula <- function(formula, data, id = NULL, wt_names = NULL,
                                    na_strategy = c("preserve", "listwise")) {
  na_strategy <- match.arg(na_strategy)

  # Drop Error() pour construire le model.frame proprement
  .drop_error_term <- function(f) {
    ftxt <- deparse(f)
    ftxt <- sub("\\+\\s*Error\\(.*\\)\\s*$", "", ftxt)
    as.formula(ftxt)
  }

  f_tmp <- if (!is.null(data)) .strip_data_dollar_safe(formula, data) else formula
  f_mf  <- .drop_error_term(f_tmp)

  # 1) Construire mf SANS suppression de NA (toujours)
  mf_full <- model.frame(
    f_mf,
    data = if (!is.null(data)) data else parent.frame(),
    na.action = na.pass,
    drop.unused.levels = FALSE
  )

  # 2) Préparer x/g sur la base de mf_full
  response    <- names(mf_full)[1L]
  predictors  <- names(mf_full)[-1L]
  g_df        <- mf_full[, predictors, drop = FALSE]

  # 3) Ne jamais laisser id / wt dans g
  drop_cols <- unique(na.omit(c(id, wt_names)))
  drop_cols <- intersect(drop_cols, names(g_df))
  if (length(drop_cols)) {
    g_df <- g_df[, setdiff(names(g_df), drop_cols), drop = FALSE]
  }

  # 4) Stratégie NA
  if (na_strategy == "listwise") {
    keep <- stats::complete.cases(mf_full)
    mf   <- mf_full[keep, , drop = FALSE]
  } else {
    keep <- rep(TRUE, nrow(mf_full))  # on ne supprime rien
    mf   <- mf_full
  }

  # 5) x / g (g reste un data.frame; si 1 col, à toi de le factoriser au besoin plus tard)
  x_out <- mf[[1L]]
  g_out <- mf[, -1L, drop = FALSE]

  structure(
    list(
      x = x_out,
      g = g_out,
      mf = mf,
      keep = keep,                 # masque des lignes conservées (utile si on doit synchroniser `data`)
      had_na = any(!complete.cases(mf_full)),
      response = response,
      predictors = predictors
    ),
    class = "norm_form"
  )
}
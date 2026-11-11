#' Sphericity-corrected p-values of the within-subject terms
#'
#' Always applied, whatever Mauchly's test says: Mauchly has little power with
#' few subjects, and the uncorrected F is liberal as soon as sphericity fails
#' (12 subjects x 6 times, AR(1) rho = 0.9: 8.7 % of rejections under H0).
#' Rule of Girden (1992): Greenhouse-Geisser when its epsilon is <= 0.75,
#' Huynh-Feldt (capped at 1) above.
#'
#' @param anova_result Output of rstatix::anova_test(detailed = TRUE).
#' @return data.frame (term, F, DFn, DFd, p_uncorrected, GGe, HFe, method,
#'   epsilon, p_corrected), or NULL when no within term has 3 levels or more.
#' @keywords internal
.sphericity_corrected <- function(anova_result) {
  # rstatix arrondit les epsilons a 3 decimales : on refait le calcul de car::Anova()
  # sur le meme plan pour des p-values exactes ; repli sur la table arrondie sinon.
  exact <- tryCatch({
    a <- attr(anova_result, "args")
    design <- rstatix::factorial_design(data = a$data, dv = a$dv, wid = a$wid,
                                        between = a$between, within = a$within)
    res <- suppressMessages(car::Anova(design$model, idata = design$idata,
                                       idesign = design$idesign, type = 3))
    sm <- suppressWarnings(summary(res, multivariate = FALSE))
    adj <- as.data.frame(unclass(sm$pval.adjustments))
    uni <- as.data.frame(unclass(sm$univariate.tests))
    if (nrow(adj) == 0) NULL else {
      j <- match(rownames(adj), rownames(uni))
      data.frame(Effect = rownames(adj), GGe = adj[["GG eps"]], HFe = adj[["HF eps"]],
                 F = uni[j, "F value"], DFn = uni[j, "num Df"], DFd = uni[j, "den Df"],
                 stringsAsFactors = FALSE)
    }
  }, error = function(e) NULL)
  if (!is.null(exact)) {
    corr <- exact
    tab <- exact
  } else {
    corr <- tryCatch(as.data.frame(anova_result$`Sphericity Corrections`), error = function(e) NULL)
    tab <- tryCatch(as.data.frame(anova_result$ANOVA), error = function(e) NULL)
  }
  if (is.null(corr) || is.null(tab) || nrow(corr) == 0 || !"Effect" %in% names(tab)) return(NULL)
  rows <- lapply(seq_len(nrow(corr)), function(i) {
    term <- as.character(corr$Effect[i])
    j <- match(term, as.character(tab$Effect))
    if (is.na(j)) return(NULL)
    gge <- as.numeric(corr$GGe[i])
    hfe <- as.numeric(corr$HFe[i])
    use_gg <- is.na(hfe) || gge <= 0.75
    eps <- if (use_gg) gge else min(1, hfe)
    f <- as.numeric(tab$F[j]); dfn <- as.numeric(tab$DFn[j]); dfd <- as.numeric(tab$DFd[j])
    data.frame(term = term, F = f, DFn = dfn, DFd = dfd,
               p_uncorrected = stats::pf(f, dfn, dfd, lower.tail = FALSE),
               GGe = gge, HFe = hfe,
               method = if (use_gg) "Greenhouse-Geisser" else "Huynh-Feldt",
               epsilon = eps,
               p_corrected = stats::pf(f, eps * dfn, eps * dfd, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

#' Replace uncorrected p-values by sphericity-corrected ones in an $effects table
#'
#' @param effects data.frame (term, p.value).
#' @param corrections Output of .sphericity_corrected().
#' @return The table with corrected p-values for the within-subject terms.
#' @keywords internal
.apply_sphericity <- function(effects, corrections) {
  if (is.null(effects) || is.null(corrections)) return(effects)
  norm <- function(s) gsub("\\s+", "", s)
  idx <- match(norm(effects$term), norm(corrections$term))
  hit <- !is.na(idx)
  effects$p.value[hit] <- corrections$p_corrected[idx[hit]]
  effects
}

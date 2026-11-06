# Agent 24 - fonctions communes (2026-11-06). Modele : voir PROJET_METHODO / consigne.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
RACINE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/"
for (f in c("sys_msg.R","sys_exit.R","sys_bay_core.R","sys_bay_oneway.R"))
  source(paste0(RACINE, "R/", f), encoding = "UTF-8")
SORTIES <- paste0(RACINE, "dev/mbay/scripts/resultats/")
HET <- .bay_bf_oneway_hetero; ZS <- .bay_bf_oneway_zs
ic <- function(x, n) { if (n == 0) return(c(NA, NA)); as.numeric(binom.test(x, n)$conf.int) }
fmt_ic <- function(x, n) { i <- ic(x, n); sprintf("%.3f [%.3f;%.3f]", x/n, i[1], i[2]) }
het_safe <- function(y, g, r = .5, n_draws = 1000L, seed = 1L) {
  tryCatch(suppressWarnings(HET(y, g, r, as.integer(n_draws), as.integer(seed))),
           error = function(e) list(log_bf = NA, error_pct = NA, ess = NA, n_draws = NA,
                                    low_ess = NA, msg = conditionMessage(e)))
}
md <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, TRUE); d <- df
  d[num] <- lapply(df[num], function(v) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "fg", flag = "#")))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(l) paste0("| ", paste(l, collapse = " | "), " |")))
}

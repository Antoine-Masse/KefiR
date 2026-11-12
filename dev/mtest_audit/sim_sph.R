# Faux positifs sous H0 en mesures repetees, avec et sans sphericite.
# Usage : Rscript sim_sph.R <pkg> <N> <out.rds> <workers>
args <- commandArgs(TRUE)
pkg <- args[1]; N <- as.integer(args[2]); out <- args[3]; workers <- as.integer(args[4]); scens <- if (length(args) > 4) args[-(1:4)] else c("ar1", "cs")
lib <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
gen <- function(i, scen) {
  set.seed(5000 + i + if (scen == "cs") 100000 else 0)
  n <- 12; k <- 6
  sig <- if (scen == "ar1") 0.9^abs(outer(1:k, 1:k, "-")) else matrix(0.5, k, k) + diag(0.5, k)
  y <- matrix(rnorm(n * k), n, k) %*% chol(sig)
  data.frame(id = factor(rep(1:n, k)), T = factor(rep(paste0("t", 1:k), each = n)), y = as.vector(y))
}
# Plan en tranches (audit s) : 2 groupes de 6 sujets x 6 temps, AR(1), aucun effet
one_split <- function(i) {
  d <- gen(i, "ar1"); d$G <- factor(rep(rep(c("A", "B"), each = 6), 6))
  e <- tryCatch(suppressWarnings(suppressMessages({
    utils::capture.output(v <- m.test(y ~ G * T + Error(id / T), data = d, boot = FALSE, plot = FALSE)); v })),
    error = function(e) NULL)
  pe <- function(term) {
    if (is.null(e)) return(NA_real_)
    x <- e$effects$p.value[e$effects$term == term]
    if (length(x)) x[1] else NA_real_
  }
  c(i = i, p_mtest = if (is.null(e)) NA_real_ else as.numeric(e$global_pvalue), p_aov = NA_real_,
    p_G = pe("G"), p_T = pe("T"), p_GT = pe("G:T"),
    hors_anova = as.numeric(is.null(e) || !grepl("ANOVA", paste(e$method, collapse = " "))))
}
one <- function(i, scen) {
  if (scen == "split") return(one_split(i))
  d <- gen(i, scen)
  p <- tryCatch(suppressWarnings(suppressMessages({
    utils::capture.output(v <- m.test(y ~ T + Error(id / T), data = d, return = FALSE, boot = FALSE)); v })),
    error = function(e) NA_real_)
  a <- summary(aov(y ~ T + Error(id / T), data = d))[["Error: id:T"]][[1]]["T", "Pr(>F)"]
  c(i = i, p_mtest = as.numeric(p[1]), p_aov = a)
}
# workers = 1 : un seul processus, sans cluster (machine partagee)
if (workers > 1) {
  cl <- parallel::makeCluster(workers)
  parallel::clusterExport(cl, c("lib", "pkg", "gen", "one"))
  invisible(parallel::clusterEvalQ(cl, { .libPaths(c(lib, .libPaths()))
    suppressMessages(devtools::load_all(pkg, quiet = TRUE)) }))
  run <- function(scen) parallel::parLapply(cl, seq_len(N), one, scen = scen)
} else {
  .libPaths(c(lib, .libPaths()))
  suppressMessages(devtools::load_all(pkg, quiet = TRUE))
  run <- function(scen) lapply(seq_len(N), one, scen = scen)
}
res <- list()
for (scen in scens) {
  r <- as.data.frame(do.call(rbind, run(scen)))
  r$scen <- scen
  res[[scen]] <- r
  saveRDS(res, out)
  if (scen == "split") {
    cat(scen, ": p globale", round(100 * mean(r$p_mtest < 0.05, na.rm = TRUE), 1),
        "%  G", round(100 * mean(r$p_G < 0.05, na.rm = TRUE), 1),
        "%  T", round(100 * mean(r$p_T < 0.05, na.rm = TRUE), 1),
        "%  G:T", round(100 * mean(r$p_GT < 0.05, na.rm = TRUE), 1),
        "%  hors ANOVA", sum(r$hors_anova), "  NA", sum(is.na(r$p_mtest)), "\n")
  } else cat(scen, ": m.test", round(100 * mean(r$p_mtest < 0.05, na.rm = TRUE), 1), "%  aov non corrige",
      round(100 * mean(r$p_aov < 0.05), 1), "%  NA", sum(is.na(r$p_mtest)), "\n")
}
if (workers > 1) parallel::stopCluster(cl)

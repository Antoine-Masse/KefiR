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
one <- function(i, scen) {
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
  saveRDS(do.call(rbind, res), out)
  cat(scen, ": m.test", round(100 * mean(r$p_mtest < 0.05, na.rm = TRUE), 1), "%  aov non corrige",
      round(100 * mean(r$p_aov < 0.05), 1), "%  NA", sum(is.na(r$p_mtest)), "\n")
}
if (workers > 1) parallel::stopCluster(cl)

source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
rep_ <- function(st, S = 4000L) suppressWarnings(.bay_hetero_core(st, .5, S, 1L))
zs_st <- function(n, m, ss) { # zs a partir de resumes : on reconstruit des donnees
  y <- unlist(lapply(seq_along(n), function(i) { e <- scale(qnorm(ppoints(n[i])))[, 1]; m[i] + e * sqrt(ss[i] / (n[i] - 1)) }))
  .bay_bf_oneway_zs(y, rep(seq_along(n), n))
}
cat("=== E. Ajout de groupes identiques, duplication\n")
base <- list(n = c(10, 10), m = c(0, 1.2), ss = c(9, 9))
for (extra in c(0, 1, 2, 5, 10, 20)) { n <- c(base$n, rep(10, extra)); m <- c(base$m, rep(0.6, extra)); ss <- c(base$ss, rep(9, extra))
  cat(sprintf("%2d groupes ajoutes au centre : log_bf hetero %7.3f  zs %7.3f\n", extra, rep_(st_of(n, m, ss))$log_bf, zs_st(n, m, ss))) }
cat("Doubler chaque observation (n=10->20, memes moyennes, SS*2+) :", rep_(st_of(c(20, 20), c(0, 1.2), c(18, 18)))$log_bf, "\n")
cat("=== F. k grand : troncature a gauche de la grille (F ~ 1, r = 1, n grand)\n")
wide <- function(st, r, S = 2000L, h = 1/16) {
  set.seed(1); s2 <- .bay_hetero_draws(st, S)
  a <- .bay_hetero_log_h1_inner(s2, st, r, h)
  # bornes elargies vers la gauche : memes tirages
  v <- sweep(s2, 2, st$n, "/"); tau2 <- rowMeans(s2); M <- matrix(rep(st$m, each = S), S, st$k)
  grid <- seq(log(r^2 / 1e5), log(r^2 / 120) - h, by = h); lp <- .5 * log(r^2 / 2) - .5 * log(pi); acc <- a - log(h)
  for (t in grid) { g <- exp(t); val <- .bay_hetero_log_kernel(v + g * tau2, M) + lp - .5 * t - r^2 / (2 * g)
    hi <- pmax(acc, val); acc <- hi + log1p(exp(pmin(acc, val) - hi)) }
  b <- acc + log(h); c(code = .bay_is_summary(a)$log_mean, large = .bay_is_summary(b)$log_mean)
}
for (cfg in list(c(50, 1000, .5), c(100, 1000, .5), c(100, 1000, 1), c(200, 1000, .5), c(200, 1000, 1), c(500, 1000, 1), c(500, 1e4, .5))) {
  k <- cfg[1]; ni <- cfg[2]; r <- cfg[3]; set.seed(k); n <- rep(ni, k); m <- rnorm(k, 0, 1 / sqrt(ni)); st <- st_of(n, m, rep(ni - 1, k))
  w <- wide(st, r); cat(sprintf("k=%d n_i=%g r=%g : H1 code %.3f bornes larges %.3f ecart %.3f\n", k, ni, r, w[1], w[2], w[2] - w[1]))
}
cat("=== G. temps de l'integrale H0 a grand n (limite 40 s)\n")
for (cfg in list(list(n = rep(1e5, 3), m = c(0, 0, 0)), list(n = rep(1e6, 3), m = c(0, 0, 0)), list(n = rep(1e6, 3), m = c(0, 1e-3, 2e-3)), list(n = rep(1e7, 3), m = c(0, 0, 0)))) {
  st <- st_of(cfg$n, cfg$m, (cfg$n - 1) * c(1, 4, 9))
  setTimeLimit(elapsed = 40, transient = TRUE)
  r <- tryCatch({ tt <- system.time(a <- .bay_hetero_h0_exact(st))[3]; paste("h0", round(a, 3), "en", round(tt, 2), "s") }, error = function(e) paste("ARRET:", conditionMessage(e)))
  setTimeLimit(); cat("n =", cfg$n[1], "m =", cfg$m, ":", r, "\n")
}

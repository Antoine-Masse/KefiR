source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
# inner avec tau au choix (copie de .bay_hetero_log_h1_inner, seule la ligne tau2 change)
inner_tau <- function(s2, st, r, h, tau = "arith") {
  k <- st$k; S <- nrow(s2); v <- sweep(s2, 2L, st$n, "/")
  tau2 <- if (tau == "arith") rowMeans(s2) else exp(rowMeans(log(s2)))
  M <- matrix(rep(st$m, each = S), S, k); ssm <- sum((st$m - mean(st$m))^2)
  t_grid <- seq(log(r^2 / 120), max(log(r^2), log(ssm) - log(min(tau2)), log(k)) + 25, by = h)
  lp <- 0.5 * log(r^2 / 2) - 0.5 * log(pi); acc <- rep(-Inf, S)
  for (t in t_grid) { g <- exp(t)
    val <- .bay_hetero_log_kernel(v + g * tau2, M) + lp - 0.5 * t - r^2 / (2 * g)
    hi <- pmax(acc, val); acc <- hi + log1p(exp(pmin(acc, val) - hi)) }
  acc + log(h)
}
bf_tau <- function(st, tau, r = .5, S = 4000L, h = 1/4, seed = 1) {
  set.seed(seed); s2 <- .bay_hetero_draws(st, S)
  sm <- .bay_is_summary(inner_tau(s2, st, r, h, tau))
  c(log_bf = .bay_hetero_log_const(st) + sm$log_mean - .bay_hetero_h0_exact(st), ess = sm$ess)
}
rep_ <- function(st, S = 4000L) suppressWarnings(.bay_hetero_core(st, .5, S, 1L))
zs_st <- function(n, m, ss) { # zs a partir de resumes : on reconstruit des donnees
  y <- unlist(lapply(seq_along(n), function(i) { e <- scale(qnorm(ppoints(n[i])))[, 1]; m[i] + e * sqrt(ss[i] / (n[i] - 1)) }))
  .bay_bf_oneway_zs(y, rep(seq_along(n), n))
}
cat("=== A. Groupe sans rapport : pair 1-2 a t ~ 3.5 (n=20, d=1.1 sigma), groupe 3 etale d'un facteur S\n")
for (S in c(1, 10, 100, 1e3, 1e4, 1e6)) {
  n <- c(20, 20, 20); m <- c(0, 1.1, .55); ss <- c(19, 19, 19 * S^2); st <- st_of(n, m, ss)
  a <- rep_(st); gm <- bf_tau(st, "geom")
  cat(sprintf("S=%8g  arith log_bf %8.3f (BF %.3g) err%% %.2f ess %5.0f n %6d | geom %8.3f ess %5.0f\n", S, a$log_bf, exp(a$log_bf), a$error_pct, a$ess, a$n_draws, gm[1], gm[2]))
}
cat(sprintf("Reference k=2, groupes 1 et 2 seuls : log_bf %.3f\n", rep_(st_of(c(20, 20), c(0, 1.1), c(19, 19)))$log_bf))
cat("=== B. Eloigner un petit groupe : BF en fonction de d, a SS fixes (n_small, k)\n")
for (cfg in list( 
                 list(n = c(30, 30, 2), lab = "n=(30,30,2)"), list(n = c(30, 30, 5), lab = "n=(30,30,5)"),
                 list(n = c(20, 20, 20, 20, 3), lab = "n=(20x4,3)"), list(n = c(20, 20, 20, 20, 5), lab = "n=(20x4,5)"))) {
  n <- cfg$n; k <- length(n); ss <- n - 1
  out <- sapply(c(3, 10, 100, 1e3, 1e4), function(d) { m <- c(rep(0, k - 1), d); st <- st_of(n, m, ss); a <- rep_(st); gm <- bf_tau(st, "geom")
    z <- tryCatch(zs_st(n, m, ss), error = function(e) NA)
    c(a$log_bf, gm[1], z, a$ess) })
  cat(cfg$lab, "\n"); dimnames(out) <- list(c("arith", "geom", "zs", "ess_arith"), paste0("d=", c(3, 10, 100, 1e3, 1e4))); print(round(out, 2))
}
cat("=== C. Paire n=2 quasi egale, autre groupe n=50 : BF selon eps\n")
set.seed(1); y2 <- rnorm(50); y2 <- (y2 - mean(y2)) / sd(y2)
for (d in c(0.5, 1.5)) for (eps in c(1, .1, 1e-2, 1e-3, 1e-5, 1e-8)) {
  a <- rep_(.bay_hetero_prep(c(d, d + eps, y2), rep(1:2, c(2, 50))))
  cat(sprintf("d=%.1f eps=%6g log_bf %7.3f err%% %.2f ess %5.0f\n", d, eps, a$log_bf, a$error_pct, a$ess))
}
cat("=== D. Consistance, resumes exacts\n")
for (N in c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)) {
  n <- rep(N, 3); h0 <- rep_(st_of(n, c(0, 0, 0), (n - 1) * c(1, 4, 9)), 5000L)
  h1 <- rep_(st_of(n, c(0, 5 / sqrt(N), 0), (n - 1) * c(1, 4, 9)), 5000L); h2 <- rep_(st_of(n, c(0, .3, 0), (n - 1) * c(1, 4, 9)), 5000L)
  dese <- rep_(st_of(c(10, N), c(0, 0), c(9, (N - 1))), 5000L)
  cat(sprintf("n=%8g H0 %8.3f (pente th -(k-1)/2 log n = %7.3f) | d=5/sqrt n %7.3f | d=.3 %9.2f | n=(10,N) H0 %7.3f\n", N, h0$log_bf, -log(N), h1$log_bf, h2$log_bf, dese$log_bf))
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

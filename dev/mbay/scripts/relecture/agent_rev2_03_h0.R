source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
# Reference de l'integrale de H0 : decoupage aux milieux entre moyennes triees, grille asinh locale
# sur u = mu - m_i, sans passer par la soustraction mu - m_i (precision complete pres du pic).
ref_h0 <- function(st, M = if (st$k <= 8) 100001 else 4001) {
  k <- st$k; n <- st$n; m <- st$m; ss <- st$ss
  const <- sum(-n / 2 * log(2 * pi) + lgamma(n / 2))
  o <- order(m); br <- c(-Inf, (m[o][-1] + m[o][-k]) / 2, Inf)
  span <- diff(range(m)) + max(sqrt(ss / n))
  tot <- numeric(k)
  for (ii in seq_len(k)) {
    i <- o[ii]; si <- sqrt(ss[i] / n[i]) 
    lo <- max(br[ii] - m[i], -1e7 * span); hi <- min(br[ii + 1] - m[i], 1e7 * span)
    x <- seq(asinh(lo / si), asinh(hi / si), length.out = M)
    u <- si * sinh(x); du <- si * cosh(x) * (x[2] - x[1])
    lf <- n[i] / 2 * (log(2) - log(ss[i] + n[i] * u^2))
    for (j in setdiff(seq_len(k), i)) lf <- lf + n[j] / 2 * (log(2) - log(ss[j] + n[j] * (u + (m[i] - m[j]))^2))
    tot[ii] <- .bay_log_sum_exp(lf + log(du))
  }
  const + .bay_log_sum_exp(tot)
}
show <- function(lbl, st) { 
  w <- NULL
  a <- withCallingHandlers(.bay_hetero_h0_exact(st), warning = function(c) {w <<- conditionMessage(c); invokeRestart("muffleWarning")})
  b <- ref_h0(st); cat(sprintf("%-46s code %.6f ref %.6f ecart %.2e %s\n", lbl, a, b, a - b, if (is.null(w)) "" else "[W]"))
}
show("regulier", st_of(c(12, 15, 9), c(0, .4, .2), c(11, 30, 5)))
show("3 pics eloignes, sd .01", st_of(c(40, 40, 40), c(0, 1e3, 5e3), rep(39 * 1e-4, 3)))
show("3 pics 1e6 sigmas", st_of(c(5, 5, 5), c(0, 1, 2) * 1e6, rep(4, 3)))
show("3 pics 1e12 sigmas", st_of(c(5, 5, 5), c(0, 1, 2) * 1e12, rep(4, 3)))
show("3 pics 1e15 sigmas", st_of(c(5, 5, 5), c(0, 1, 2) * 1e15, rep(4, 3)))
show("3 pics 1e18 sigmas", st_of(c(5, 5, 5), c(0, 1, 2) * 1e18, rep(4, 3)))
show("n=2, 2 ; 1e8", st_of(c(2, 2), c(0, 1e8), c(1, 1)))
show("n=2,2 ; 1e16", st_of(c(2, 2), c(0, 1e16), c(1, 1)))
show("sd tres inegaux 1e-6 / 1e6", st_of(c(10, 10), c(0, 1), c(9e-12, 9e12)))
show("sd tres inegaux 1e-12 / 1 gap 1", st_of(c(10, 10), c(0, 1), c(9e-24, 9)))
show("grand N 1e6 x3 gap 1", st_of(c(1e6, 1e6, 1e6), c(0, 1, 3), rep(1e6, 3)))
show("N 1e6 chevauchants", st_of(c(1e6, 1e6, 1e6), c(0, 1e-3, 2e-3), rep(1e6, 3)))
show("k=30 n=2", {set.seed(1); st_of(rep(2, 30), rnorm(30), rchisq(30, 1))})
show("k=200 n=5 sd inegaux", {set.seed(2); st_of(rep(5, 200), rnorm(200, 0, 3), rchisq(200, 4) * exp(rnorm(200, 0, 3)))})
show("k=500 n=1000", {set.seed(3); st_of(rep(1e3, 500), rnorm(500, 0, 1e-2), rep(1e3, 500))})
set.seed(4); worst <- 0
for (rep_ in 1:40) {
  k <- sample(2:8, 1); n <- sample(c(2, 3, 5, 20, 300, 1e4), k, TRUE)
  st <- st_of(n, rnorm(k) * 10^runif(1, -3, 6), (n - 1) * 10^runif(k, -6, 6))
  a <- suppressWarnings(.bay_hetero_h0_exact(st)); b <- ref_h0(st, 40001)
  if (abs(a - b) > 1e-6) cat("ECART", round(a - b, 6), "n", n, "\n"); worst <- max(worst, abs(a - b))
}
cat("pire ecart sur 40 tirages aleatoires", worst, "\n")
cat("--- precision : pic plus etroit que la resolution de m apres centrage\n")
for (sd_ in c(1e-6, 1e-10, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-20, 1e-30)) {
  set.seed(5); yA <- rnorm(5, 0, sd_); yB <- rnorm(5, 1, .1)
  y <- c(yA, yB); g <- rep(1:2, each = 5); st <- .bay_hetero_prep(y, g)
  a <- suppressWarnings(.bay_hetero_h0_exact(st)); b <- ref_h0(st)
  res <- suppressWarnings(.bay_hetero_core(st, .5, 2000L, 1L))
  cat(sprintf("sd_A=%.0e  H0 code %.3f ref %.3f ecart %.3f | log_bf code %.3f  log_bf ref %.3f\n", sd_, a, b, a - b, res$log_bf, res$log_ml_h1 - b))
}
cat("--- groupes (0.3 et 0.1+0.2) : variance nulle, moyennes differentes de 1 ulp\n")
print(tryCatch(.bay_bf_oneway_zs(c(.3, .3, .3, .1 + .2, .1 + .2, .1 + .2), rep(1:2, each = 3)), error = function(e) conditionMessage(e)))

source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
cat("=== temps et pas de grille\n")
for (cfg in list(list(n = rep(30, 3), lab = "k=3 n=30"), list(n = rep(10, 10), lab = "k=10"), list(n = rep(10, 30), lab = "k=30"), list(n = rep(5, 60), lab = "k=60"))) {
  set.seed(1); k <- length(cfg$n)
  for (d in c(0, 1e3)) { m <- rnorm(k) * d / sqrt(cfg$n[1]) ; st <- st_of(cfg$n, m + rnorm(k) / sqrt(cfg$n[1]), cfg$n - 1)
    tt <- system.time(a <- suppressWarnings(.bay_hetero_core(st, .5, 20000L, 1L)))[3]
    set.seed(9); s2 <- .bay_hetero_draws(st, 3000L)
    A <- .bay_hetero_log_h1_inner(s2, st, .5, a$grid_step); B <- .bay_hetero_log_h1_inner(s2, st, .5, a$grid_step / 4)
    cat(sprintf("%-10s d=%5g %6.1f s h=%.4f n=%6d grow=%d ess=%6.0f err%%=%5.2f low=%s | max|Delta grille|=%.2e\n", cfg$lab, d, tt, a$grid_step, a$n_draws, a$n_grow, a$ess, a$error_pct, a$low_ess, max(abs(A - B))))
  }
}
cat("=== calibration sous H0 et puissance : n_draws = 1000, 150 jeux, meme graine\n")
sim <- function(n, sds, eff, R = 150, seed = 100) {
  set.seed(seed); out <- matrix(NA_real_, R, 3)
  g <- rep(seq_along(n), n)
  for (i in seq_len(R)) {
    y <- rnorm(sum(n), rep(eff, n), rep(sds, n))
    a <- tryCatch(suppressWarnings(.bay_bf_oneway_hetero(y, g, n_draws = 1000L, seed = i)$log_bf), error = function(e) NA)
    out[i, ] <- c(a, .bay_bf_oneway_zs(y, g), oneway.test(y ~ factor(g))$p.value)
  }
  bf <- exp(out[, 1])
  cat(sprintf("n=(%s) sd=(%s) eff=(%s): hetero P(BF>3)=%.3f P(BF>10)=%.3f P(BF>19)=%.3f P(BF<1/3)=%.3f E[BF]=%.3f med=%.3f | zs P(BF>3)=%.3f P(BF>19)=%.3f | Welch p<.05=%.3f p<.005=%.3f\n",
    paste(n, collapse = ","), paste(sds, collapse = ","), paste(eff, collapse = ","),
    mean(bf > 3, na.rm = TRUE), mean(bf > 10, na.rm = TRUE), mean(bf > 19, na.rm = TRUE), mean(bf < 1/3, na.rm = TRUE), mean(bf, na.rm = TRUE), median(bf, na.rm = TRUE),
    mean(out[, 2] > log(3)), mean(out[, 2] > log(19)), mean(out[, 3] < .05), mean(out[, 3] < .005)))
}
sim(c(5, 5, 40), c(4, 4, 1), c(0, 0, 0))
sim(c(2, 50), c(1, 1), c(0, 0))
sim(c(2, 50), c(10, 1), c(0, 0))
sim(c(2, 50), c(.1, 1), c(0, 0))
sim(c(2, 2, 40), c(4, 4, 1), c(0, 0, 0))
sim(c(3, 3, 40), c(4, 4, 1), c(0, 0, 0))
sim(c(2, 2, 2, 2, 2, 2), rep(1, 6), rep(0, 6))
sim(c(2, 40, 40), c(.2, 1, 1), c(0, 0, 0))
sim(c(5, 5, 40), c(4, 4, 1), c(4, 0, 0))
sim(c(5, 5, 40), c(4, 4, 1), c(8, 0, 0))
sim(c(30, 30, 30), c(1, 1, 1), c(0, .8, 0))
sim(c(30, 30, 30), c(1, 1, 20), c(0, .8, 0))
sim(c(10, 10, 10, 10, 3), c(1, 1, 1, 1, 1), c(0, 0, 0, 0, 3))

source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
# 1. Reference quadrature contre IS sur un cas regulier (les constantes)
st <- st_of(c(8, 12), c(0, 1), c(7, 22))
t0 <- Sys.time(); rb <- ref_bf(st, h = 1/8, span_hi = 30, step = .3); cat("ref", rb, "temps", format(Sys.time() - t0), "\n")
a <- .bay_hetero_core(st, .5, 20000L, 1L); cat("IS", a$log_bf, "err%", a$error_pct, "ess", a$ess, "\n")
# 2. Reference totalement independante : Monte-Carlo sur le modele complet, conditionne a mu,a,s
# k = 2, n petits : quadrature brute 4D evitee ; on integre mu, delta analytiquement impossible ->
# verification par integrate imbriquee pour H0 ; H1 par s2 tires du POSTERIEUR H0 (Gibbs).
gibbs_ratio <- function(st, r = .5, S = 40000, burn = 2000, h = 1/8, seed = 3) {
  set.seed(seed)
  mu <- sum(st$n * st$m) / st$N
  out <- matrix(0, S, st$k)
  for (it in seq_len(S + burn)) {
    s2 <- (st$ss + st$n * (st$m - mu)^2) / rchisq(st$k, st$n)
    w <- st$n / s2
    mu <- rnorm(1, sum(w * st$m) / sum(w), 1 / sqrt(sum(w)))
    if (it > burn) out[it - burn, ] <- s2
  }
  v <- sweep(out, 2, st$n, "/")
  M <- matrix(rep(st$m, each = S), S, st$k)
  l1 <- .bay_hetero_log_h1_inner(out, st, r, h)
  l0 <- .bay_hetero_log_kernel(v, M)
  list(log_bf = .bay_is_summary(l1 - l0)$log_mean, ess = .bay_is_summary(l1 - l0)$ess)
}
g <- gibbs_ratio(st); cat("Gibbs H0 ratio", g$log_bf, "ess", g$ess, "\n")
st3 <- st_of(c(5, 9, 7), c(0, .8, -.3), c(5, 30, 2))
cat("k=3 IS", .bay_hetero_core(st3, .5, 20000L, 1L)$log_bf, " Gibbs", gibbs_ratio(st3)$log_bf, "
")
# 3. delta prior k=2 : quantile |a1-a2|/tau = Cauchy(sqrt2 r)
set.seed(1); gg <- 0.25 / rchisq(4e5, 1); dd <- sqrt(2 * gg) * rnorm(4e5); cat("mediane |delta|", median(abs(dd)), "\n")
# 4. invariance par permutation des etiquettes
d <- list(y = c(rnorm(10), rnorm(20, 1, 2), rnorm(5, 0, .5)), g = rep(c("a", "b", "c"), c(10, 20, 5)))
set.seed(2); d$y <- c(rnorm(10), rnorm(20, 1, 2), rnorm(5, 0, .5))
a1 <- .bay_bf_oneway_hetero(d$y, d$g)$log_bf
a2 <- .bay_bf_oneway_hetero(d$y, chartr("abc", "cab", d$g))$log_bf
cat("permutation etiquettes", a1, a2, "\n")
# 5. Plafond de P(H0|D) a F = 0, annonce au projet : 0,720 / 0,859 / 0,909 (3 groupes, r = 0,5)
for (ni in c(3, 15, 30)) { e <- scale(qnorm(ppoints(ni)))[, 1]; lb <- .bay_bf_oneway_zs(rep(e, 3), rep(1:3, each = ni)); cat("n_i", ni, "P(H0|F=0)", round(1 / (1 + exp(lb)), 3), "\n") }

source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
cat("--- plans equilibres : un groupe eloigne de d, k groupes de n = 5\n")
for (k in c(4, 6, 8)) for (d in c(10, 100, 1e3)) { n <- rep(5, k); a <- suppressWarnings(.bay_hetero_core(st_of(n, c(rep(0, k - 1), d), n - 1), .5, 4000L, 1L))
  cat(sprintf("k=%d n=5 d=%5g log_bf %8.3f err%% %.2f ess %.0f\n", k, d, a$log_bf, a$error_pct, a$ess)) }
cat("--- donnees brutes, reproductible : 4 groupes n=10 N(0,1), 5e groupe n=3 autour de 50\n")
set.seed(7); y <- c(rnorm(40), rnorm(3, 50, 1)); g <- rep(c("a","b","c","d","e"), c(10,10,10,10,3))
a <- .bay_bf_oneway_hetero(y, g); cat("log_bf", a$log_bf, "err%", a$error_pct, "ess", a$ess, " zs", .bay_bf_oneway_zs(y, g), " Welch p", oneway.test(y ~ factor(g))$p.value, "\n")
cat("--- ESS menteur : n_draws croissant, r = 0.005, 0.05, 0.1, 0.5\n")
st <- st_of(c(2, 30), c(0, 1e4), c(1e-4, 29))
for (r in c(.005, .05, .1, .5)) { rb <- ref_bf(st, r, h = 1/4, span_hi = 50, step = .3); cat("r", r, "REF", round(rb, 3), "\n")
  for (S in if (r == .005) c(2000L, 20000L, 200000L, 800000L) else c(20000L)) {
    a <- suppressWarnings(.bay_hetero_core(st, r, S, 1L, max_grow = 0L))
    cat(sprintf("   S=%7d log_bf %7.3f ecart %+6.3f err%% %5.2f ess %7.0f\n", a$n_draws, a$log_bf, a$log_bf - rb, a$error_pct, a$ess)) } }

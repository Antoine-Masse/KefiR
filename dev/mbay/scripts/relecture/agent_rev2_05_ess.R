source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
run <- function(lbl, st, r = .5, Ss = c(4000L, 20000L)) {
  rb <- ref_bf(st, r, h = 1/4, span_hi = 50, step = .3)
  cat(sprintf("%-34s r=%.3g REF %.3f\n", lbl, r, rb))
  for (S in Ss) for (sd_ in 1:2) {
    a <- suppressWarnings(.bay_hetero_core(st, r, S, sd_))
    cat(sprintf("   S=%6d seed %d n=%6d log_bf %7.3f ecart %+7.3f err%% %6.2f z %+6.1f ess %7.0f %s\n", S, sd_, a$n_draws, a$log_bf, a$log_bf - rb, a$error_pct, (a$log_bf - rb) / (a$error_pct / 100), a$ess, if (a$low_ess) "LOW" else ""))
  }
}
run("n=(2,2) d=1e3", st_of(c(2, 2), c(0, 1e3), c(1, 1)))
run("n=(2,2) d=1e6", st_of(c(2, 2), c(0, 1e6), c(1, 1)))
run("n=(3,3) d=1e3", st_of(c(3, 3), c(0, 1e3), c(2, 2)))
run("n=(2,30) d=1e4 SS1=1e-4 r=.05", st_of(c(2, 30), c(0, 1e4), c(1e-4, 29)), .05)
run("n=(2,30) d=1e4 SS1=1e-4 r=.005", st_of(c(2, 30), c(0, 1e4), c(1e-4, 29)), .005)
cat("--- k = 3 : quadrature 3D independante du tirage pour le paradoxe de l'eloignement\n")
st <- st_of(c(30, 30, 2), c(0, 0, 100), c(29, 29, 1))
cat("n=(30,30,2) d=100 REF", ref_bf(st, h = 1/4, span_hi = 30, step = .4), " IS", suppressWarnings(.bay_hetero_core(st, .5, 20000L, 1L))$log_bf, "\n")
cat("--- plans equilibres : un groupe eloigne de d, k groupes de n = 5\n")
for (k in c(4, 6, 8)) for (d in c(10, 100, 1e3)) { n <- rep(5, k); a <- suppressWarnings(.bay_hetero_core(st_of(n, c(rep(0, k - 1), d), n - 1), .5, 4000L, 1L))
  cat(sprintf("k=%d n=5 d=%5g log_bf %8.3f err%% %.2f ess %.0f\n", k, d, a$log_bf, a$error_pct, a$ess)) }

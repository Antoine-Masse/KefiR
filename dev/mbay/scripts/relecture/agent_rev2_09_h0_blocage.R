source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
for (nn in c(2e5, 4e5, 6e5, 8e5)) { st <- st_of(rep(nn, 3), c(0, 0, 0), (nn - 1) * c(1, 4, 9))
  setTimeLimit(elapsed = 20, transient = TRUE)
  r <- tryCatch({ tt <- system.time(a <- .bay_hetero_h0_exact(st))[3]; paste("h0", round(a, 3), "en", round(tt, 2), "s") }, error = function(e) "ARRET: limite 20 s atteinte")
  setTimeLimit(); cat("n_i =", nn, ":", r, "\n") }
cat("--- ulp : variance nulle exacte, moyennes a 1 ulp, y quasi constant\n")
print(tryCatch(.bay_bf_oneway_zs(c(.3, .3, .3, .1 + .2, .1 + .2, .1 + .2), rep(1:2, each = 3)), error = function(e) conditionMessage(e)))
st <- st_of(c(30, 30, 2), c(0, 0, 100), c(29, 29, 1))
cat("n=(30,30,2) d=100 REF 3D", ref_bf(st, h = 1/2, span_hi = 26, step = .5), " IS", suppressWarnings(.bay_hetero_core(st, .5, 20000L, 1L))$log_bf, "\n")

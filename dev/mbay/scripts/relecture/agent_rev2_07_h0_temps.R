source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev2_lib.R")
cat("=== G. temps de l'integrale H0 a grand n (limite 40 s)\n")
for (cfg in list(list(n = rep(1e5, 3), m = c(0, 0, 0)), list(n = rep(1e6, 3), m = c(0, 0, 0)), list(n = rep(1e6, 3), m = c(0, 1e-3, 2e-3)), list(n = rep(1e7, 3), m = c(0, 0, 0)))) {
  st <- st_of(cfg$n, cfg$m, (cfg$n - 1) * c(1, 4, 9))
  setTimeLimit(elapsed = 40, transient = TRUE)
  r <- tryCatch({ tt <- system.time(a <- .bay_hetero_h0_exact(st))[3]; paste("h0", round(a, 3), "en", round(tt, 2), "s") }, error = function(e) paste("ARRET:", conditionMessage(e)))
  setTimeLimit(); cat("n =", cfg$n[1], "m =", cfg$m, ":", r, "\n")
}

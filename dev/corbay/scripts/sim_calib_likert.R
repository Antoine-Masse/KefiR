# Calibration under the prior of Kendall (tau-b) on a 1 to 5 scale: latent normal pairs cut in 5 classes.
source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
ct <- function(v) as.integer(cut(v, c(-Inf, -1.28, -0.52, 0.52, 1.28, Inf)))
set.seed(13); out <- list()
for (n in c(20, 50, 100)) for (r in 1:1500) {
	h1 <- runif(1) < 0.5; rho <- if (h1) 2 * rbeta(1, 3, 3) - 1 else 0
	x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
	out[[length(out) + 1]] <- data.frame(n = n, h1 = h1, pK = 1 - post3(lbf_kendall(ct(x), ct(y)))[["none"]])
}
cl <- do.call(rbind, out); saveRDS(cl, "resultats/sim_calib_likert.rds")
for (n in unique(cl$n)) { s <- cl[cl$n == n, ]; b <- cut(s$pK, c(0, .2, .4, .6, .8, 1), include.lowest = TRUE)
	cat("n =", n, "\n"); print(data.frame(announced = round(tapply(s$pK, b, mean), 3), freq_h1 = round(tapply(s$h1, b, mean), 3), count = as.vector(table(b)))) }

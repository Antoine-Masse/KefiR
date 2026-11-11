# Estimate under the test prior (kappa = 1/3) against Jeffreys's estimation prior (kappa = 2):
# median and 95 % HDI, small n, weak to strong links; coverage of the true rho over 200 sets.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0", quiet = TRUE))
one <- function(m, x, y, k) unlist(.corbay_estimate(m, x, y, k, "two.sided", 4000L, 1L, 0.95)[c("median", "lo", "hi")])
set.seed(3); out <- list()
for (rho in c(0.3, 0.6, 0.9)) for (n in c(8, 15, 30)) for (r in 1:60) {
	x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
	for (m in c("pearson", "kendall")) {
		tr <- if (m == "kendall") 2 / pi * asin(rho) else rho
		a <- one(m, x, y, 1/3); b <- one(m, x, y, 2)
		obs <- if (m == "kendall") cor(x, y, method = "kendall") else cor(x, y)
		out[[length(out) + 1]] <- data.frame(m = m, rho = rho, n = n, obs = obs, med13 = a[1], cov13 = a[2] <= tr & tr <= a[3], in13 = a[2] <= obs & obs <= a[3], w13 = a[3] - a[2],
			med2 = b[1], cov2 = b[2] <= tr & tr <= b[3], in2 = b[2] <= obs & obs <= b[3], w2 = b[3] - b[2], tr = tr)
	}
}
o <- do.call(rbind, out); saveRDS(o, "C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/resultats/check_estim.rds")
f <- function(v) round(100 * mean(v), 0)
print(do.call(rbind, lapply(split(o, list(o$m, o$rho, o$n)), function(s) data.frame(m = s$m[1], rho = s$rho[1], n = s$n[1], true = round(s$tr[1], 2),
	obs = round(mean(s$obs), 2), med_k13 = round(mean(s$med13), 2), med_k2 = round(mean(s$med2), 2),
	cover_k13 = f(s$cov13), cover_k2 = f(s$cov2), obs_in_k13 = f(s$in13), obs_in_k2 = f(s$in2), width_k13 = round(mean(s$w13), 2), width_k2 = round(mean(s$w2), 2)))), row.names = FALSE)

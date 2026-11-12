# Estimate (bilan, Jeffreys) at rho = 0 and small n: centring and width; Kendall with the
# standard deviation of T* under H1 (eq. 14 bound) against sd = 1. 200 sets per cell.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0", quiet = TRUE))
E <- function(m, x, y, k, v = FALSE, s) unlist(.corbay_estimate(m, x, y, k, "two.sided", 2000L, s, 0.95, var_h1 = v)[c("median", "lo", "hi")])
set.seed(77); out <- list()
for (rho in c(0, 0.3, 0.6, 0.9)) for (n in c(10, 20, 30)) {
	for (r in 1:200) {
		x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n); tt <- 2 / pi * asin(rho)
		pj <- E("pearson", x, y, 2, s = r); pt <- E("pearson", x, y, 1/3, s = r)
		kj <- E("kendall", x, y, 2, s = r); kv <- E("kendall", x, y, 2, TRUE, s = r); kt <- E("kendall", x, y, 1/3, s = r)
		row <- function(lab, e, tr) data.frame(rho = rho, n = n, est = lab, med = e[1], lo = e[2], hi = e[3], cover = e[2] <= tr & tr <= e[3], width = e[3] - e[2])
		out[[length(out) + 1]] <- rbind(row("pearson_jeffreys", pj, rho), row("pearson_test", pt, rho), row("kendall_jeffreys", kj, tt),
			row("kendall_jeffreys_varH1", kv, tt), row("kendall_test", kt, tt))
	}
	cat(rho, n, "\n")
}
o <- do.call(rbind, out); saveRDS(o, "C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/resultats/check_estim2.rds")
f <- function(v) round(100 * mean(v), 1)
t <- do.call(rbind, lapply(split(o, list(o$est, o$rho, o$n)), function(s) data.frame(est = s$est[1], rho = s$rho[1], n = s$n[1],
	cover = f(s$cover), width = round(mean(s$width), 2), mean_med = round(mean(s$med), 3), mean_abs_med = round(mean(abs(s$med)), 3),
	q90_abs_med = round(quantile(abs(s$med), 0.9), 3))))
print(t[order(t$est, t$rho, t$n), ], row.names = FALSE)

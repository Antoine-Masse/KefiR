# Effect of kappa (1/3 here, 1 in JASP) on false "yes" and power, normal data. Output: resultats/sim_kappa.rds
source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
set.seed(21); res <- list()
for (rho in c(0, 0.3, 0.5)) for (n in c(20, 50)) for (r in 1:300) {
	x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
	res[[length(res) + 1]] <- data.frame(rho = rho, n = n,
		p13 = 1 - post3(lbf_pearson(x, y, 1/3))[["none"]], p1 = 1 - post3(lbf_pearson(x, y, 1))[["none"]],
		k13 = 1 - post3(lbf_kendall(x, y, 1/3))[["none"]], k1 = 1 - post3(lbf_kendall(x, y, 1))[["none"]])
}
r <- do.call(rbind, res); saveRDS(r, "resultats/sim_kappa.rds")
f <- function(v) round(100 * mean(v), 1)
print(do.call(rbind, lapply(split(r, list(r$rho, r$n)), function(s) data.frame(rho = s$rho[1], n = s$n[1],
	P13_95 = f(s$p13 > .95), P1_95 = f(s$p1 > .95), P13_75 = f(s$p13 > .75), P1_75 = f(s$p1 > .75),
	H0_13_75 = f(s$p13 < .25), H0_1_75 = f(s$p1 < .25), K13_95 = f(s$k13 > .95), K1_95 = f(s$k1 > .95)))))

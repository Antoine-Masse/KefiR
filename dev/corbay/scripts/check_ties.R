source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
ct <- function(v) as.integer(cut(v, c(-Inf, -1.28, -0.52, 0.52, 1.28, Inf)))
# S statistic and its variance under H0 with ties (Kendall 1970; the formula of cor.test, exact = FALSE)
s_ties <- function(x, y) {
	n <- length(x); S <- 0
	S <- sum(sign(outer(x, x, "-")) * sign(outer(y, y, "-"))) / 2
	tx <- table(x); ty <- table(y)
	v0 <- n * (n - 1) * (2 * n + 5)
	vt <- sum(tx * (tx - 1) * (2 * tx + 5)); vu <- sum(ty * (ty - 1) * (2 * ty + 5))
	v1 <- sum(tx * (tx - 1)) * sum(ty * (ty - 1)); v2 <- sum(tx * (tx - 1) * (tx - 2)) * sum(ty * (ty - 1) * (ty - 2))
	var <- (v0 - vt - vu) / 18 + v1 / (2 * n * (n - 1)) + v2 / (9 * n * (n - 1) * (n - 2))
	c(S = S, var = var, z_ties = S / sqrt(var), tstar = kendall_tstar(x, y))
}
set.seed(3)
for (n in c(20, 50, 100)) {
	z <- t(replicate(2000, s_ties(ct(rnorm(n)), ct(rnorm(n)))))
	cat(n, "sd T* (no-tie scaling):", round(sd(z[, "tstar"]), 3), " sd z with tie variance:", round(sd(z[, "z_ties"]), 3),
		" ratio tau-b scaled / S scaled:", round(mean(z[,"tstar"]/z[,"z_ties"], na.rm=TRUE), 3), "\n")
}
z <- t(replicate(2000, s_ties(rnorm(50), rnorm(50)))); cat("normal 50: sd T*", round(sd(z[, "tstar"]), 3), "\n")

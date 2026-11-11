# Numerical agreement: two-sided BF from the two halves against the direct correlationBF call;
# Kendall integral against the trapezoid rule; over 200 data sets.
source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
set.seed(41); d1 <- d2 <- numeric(0)
for (i in 1:200) {
	n <- sample(c(5, 10, 20, 50, 200), 1); rho <- runif(1, -0.8, 0.8)
	x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
	lp <- lbf_pearson(x, y); ld <- extractBF(correlationBF(y, x, rscale = 1/3), logbf = TRUE)$bf
	d1 <- c(d1, lp[["two"]] - ld)
	ts <- kendall_tstar(x, y); g <- seq(-1, 1, length.out = 400001)
	f <- exp(dnorm(ts, 1.5 * sqrt(n) * g, 1, log = TRUE) - dnorm(ts, 0, 1, log = TRUE)) * prior_tau(g, 1/3)
	d2 <- c(d2, lbf_kendall(x, y)[["two"]] - log(sum((f[-1] + f[-length(f)]) / 2) * diff(g[1:2])))
}
cat("Pearson halves vs direct, |diff log BF|: median", signif(median(abs(d1)), 2), "max", signif(max(abs(d1)), 2), "\n")
cat("Kendall integrate vs trapezoid, |diff log BF|: median", signif(median(abs(d2)), 2), "max", signif(max(abs(d2)), 2), "\n")

# Power compared with cor.test, and calibration under the prior. Output: resultats/sim_power.rds, sim_calib.rds
source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
args <- commandArgs(TRUE); only <- if (length(args) > 1) strsplit(args[2], ",")[[1]] else NULL; sfx <- if (length(args) > 2) args[3] else ""; R <- if (length(args)) as.integer(args[1]) else 300
bvn <- function(n, rho) { x <- rnorm(n); list(x = x, y = rho * x + sqrt(1 - rho^2) * rnorm(n)) }
gen <- list(
	normal = function(n, rho) bvn(n, rho),
	t3_biv = function(n, rho) { d <- bvn(n, rho); s <- sqrt(3 / rchisq(n, 3)); list(x = d$x * s, y = d$y * s) },
	monotone = function(n, rho) { d <- bvn(n, rho); list(x = d$x, y = exp(2 * d$y)) },
	likert = function(n, rho) { d <- bvn(n, rho); ct <- function(v) as.integer(cut(v, c(-Inf, -1.28, -0.52, 0.52, 1.28, Inf))); list(x = ct(d$x), y = ct(d$y)) }
)
set.seed(11); res <- list()
for (sc in if (is.null(only)) names(gen) else only) for (rho in c(0.2, 0.3, 0.5)) for (n in c(20, 50, 100)) {
	for (r in seq_len(R)) {
		d <- gen[[sc]](n, rho)
		res[[length(res) + 1]] <- data.frame(sc = sc, rho = rho, n = n,
			pP = post3(lbf_pearson(d$x, d$y))[["plus"]], pK = post3(lbf_kendall(d$x, d$y))[["plus"]], pS = post3(lbf_spearman(d$x, d$y))[["plus"]],
			cP = cor.test(d$x, d$y)$p.value, cK = suppressWarnings(cor.test(d$x, d$y, method = "kendall")$p.value),
			cS = suppressWarnings(cor.test(d$x, d$y, method = "spearman")$p.value), rP = cor(d$x, d$y))
	}
	cat(sc, rho, n, "\n")
}
saveRDS(do.call(rbind, res), paste0("resultats/sim_power", sfx, ".rds"))
# Calibration under the prior: half the sets without link, half with rho from the stretched beta(3, 3).
set.seed(12); cal <- list()
for (n in c(10, 20, 50)) for (r in seq_len(4 * R)) {
	h1 <- runif(1) < 0.5; rho <- if (h1) 2 * rbeta(1, 3, 3) - 1 else 0
	d <- bvn(n, rho)
	cal[[length(cal) + 1]] <- data.frame(n = n, h1 = h1, rho = rho,
		pP = 1 - post3(lbf_pearson(d$x, d$y))[["none"]], pK = 1 - post3(lbf_kendall(d$x, d$y))[["none"]], pS = 1 - post3(lbf_spearman(d$x, d$y))[["none"]])
}
saveRDS(do.call(rbind, cal), "resultats/sim_calib.rds")

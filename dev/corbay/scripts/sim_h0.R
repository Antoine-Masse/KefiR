# Calibration under H0 (no link) and diagnostics of the auto rule. Output: resultats/sim_h0.rds
source("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts/proto.R")
diag_flags <- function(x, y) {
	far1 <- function(v) { q <- quantile(v, c(.25, .75), names = FALSE); i <- q[2] - q[1]; i > 0 && any(v < q[1] - 3 * i | v > q[2] + 3 * i) }
	rob <- tryCatch(withr::with_preserve_seed({ set.seed(1); cr <- MASS::cov.rob(cbind(x, y), method = "mcd"); d <- mahalanobis(cbind(x, y), cr$center, cr$cov); any(d > qchisq(0.999, 2)) }), error = function(e) NA)
	c(far = far1(x) || far1(y), biv = rob, ties = max(1 - length(unique(x)) / length(x), 1 - length(unique(y)) / length(y)))
}
gen <- list(
	normal = function(n) list(x = rnorm(n), y = rnorm(n)),
	t3_indep = function(n) list(x = rt(n, 3), y = rt(n, 3)),
	t3_biv = function(n) { s <- sqrt(3 / rchisq(n, 3)); list(x = rnorm(n) * s, y = rnorm(n) * s) },
	lognorm = function(n) list(x = exp(rnorm(n)), y = exp(rnorm(n))),
	likert = function(n) list(x = sample(1:5, n, TRUE, c(.1, .2, .4, .2, .1)), y = sample(1:5, n, TRUE, c(.1, .2, .4, .2, .1))),
	outlier = function(n) { x <- rnorm(n); y <- rnorm(n); x[n] <- 6; y[n] <- 6; list(x = x, y = y) }
)
args <- commandArgs(TRUE); only <- if (length(args) > 1 && nzchar(args[2]) && args[2] != "all") strsplit(args[2], ",")[[1]] else NULL; sfx <- if (length(args) > 2) args[3] else ""; R <- if (length(args)) as.integer(args[1]) else 400
res <- list()
set.seed(2026)
for (sc in if (is.null(only)) names(gen) else only) for (n in c(10, 20, 50, 100)) {
	for (r in seq_len(R)) {
		d <- gen[[sc]](n)
		while (length(unique(d$x)) < 2 || length(unique(d$y)) < 2) d <- gen[[sc]](n)
		lp <- lbf_pearson(d$x, d$y); lk <- lbf_kendall(d$x, d$y); ls <- lbf_spearman(d$x, d$y)
		fl <- diag_flags(d$x, d$y)
		res[[length(res) + 1]] <- data.frame(sc = sc, n = n, r = r,
			pP = 1 - post3(lp)[["none"]], pK = 1 - post3(lk)[["none"]], pS = 1 - post3(ls)[["none"]],
			cP = cor.test(d$x, d$y)$p.value, cK = suppressWarnings(cor.test(d$x, d$y, method = "kendall")$p.value),
			cS = suppressWarnings(cor.test(d$x, d$y, method = "spearman")$p.value),
			far = fl[["far"]], biv = fl[["biv"]], ties = fl[["ties"]])
	}
	cat(sc, n, "\n")
}
res <- do.call(rbind, res)
dir.create("resultats", showWarnings = FALSE)
saveRDS(res, paste0("resultats/sim_h0", sfx, ".rds"))

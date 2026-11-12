# Moderate links (rho = 0.2, 0.3; n = 30, 50, 100): how often does the bilan put "no link" ahead
# with moderate or strong support while the credible interval of rho (bilan, Jeffreys) excludes 0?
# Package code: checks, rule, Bayes factors, support as in the bilan, estimate. Args: R, rho, suffix.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0", quiet = TRUE))
args <- commandArgs(TRUE); R <- as.integer(args[1]); rhos <- as.numeric(strsplit(args[2], ",")[[1]]); sfx <- args[3]
set.seed(round(100 * rhos[1]) + 51); res <- list()
for (rho in rhos) for (n in c(30, 50, 100)) {
	for (r in seq_len(R)) {
		x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
		ck <- .corbay_checks(x, y, 1L); ch <- .corbay_choose(ck)
		lb <- .corbay_lbf(ch$method, x, y, 1/3); h <- .corbay_hyp(lb, "two.sided", 0.5, "x", "y")
		li <- .corbay_lead_info(list(hypotheses = h, sensitivity = NULL, contradiction = FALSE))
		e <- .corbay_estimate(ch$method, x, y, .corbay_kappa_est(), "two.sided", 2000L, r, 0.95)
		et <- .corbay_estimate(ch$method, x, y, 1/3, "two.sided", 2000L, r, 0.95)
		res[[length(res) + 1]] <- data.frame(rho = rho, n = n, meth = ch$method, lead = h$key[1], p_lead = h$posterior[1], level = li$sup$level,
			lo = e$lo, hi = e$hi, lo_t = et$lo, hi_t = et$hi)
	}
	cat(rho, n, "\n")
}
saveRDS(do.call(rbind, res), paste0("resultats/sim_null_hdi", sfx, ".rds"))

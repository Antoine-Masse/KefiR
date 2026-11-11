# The rule of method = "auto" as coded in the package (checks, choice, Bayes factors, cross-check),
# under H0 and H1, against cor.test. Output: resultats/sim_auto.rds
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0", quiet = TRUE))
args <- commandArgs(TRUE); only <- if (length(args) > 1 && nzchar(args[2]) && args[2] != "all") strsplit(args[2], ",")[[1]] else NULL; sfx <- if (length(args) > 2) args[3] else ""; R <- if (length(args)) as.integer(args[1]) else 300
bvn <- function(n, rho) { x <- rnorm(n); list(x = x, y = rho * x + sqrt(1 - rho^2) * rnorm(n)) }
ct <- function(v) as.integer(cut(v, c(-Inf, -1.28, -0.52, 0.52, 1.28, Inf)))
gen <- list(
	normal = function(n, rho) bvn(n, rho),
	t3_biv = function(n, rho) { d <- bvn(n, rho); s <- sqrt(3 / rchisq(n, 3)); list(x = d$x * s, y = d$y * s) },
	t3_indep_margins = function(n, rho) { d <- bvn(n, rho); list(x = qt(pnorm(d$x), 3), y = qt(pnorm(d$y), 3)) },
	monotone = function(n, rho) { d <- bvn(n, rho); list(x = d$x, y = exp(2 * d$y)) },
	likert = function(n, rho) { d <- bvn(n, rho); list(x = ct(d$x), y = ct(d$y)) },
	outlier = function(n, rho) { d <- bvn(n, rho); d$x[n] <- 6; d$y[n] <- -6 * sign(rho + 1e-9) * (rho > 0) + 6 * (rho == 0); d }
)
one <- function(x, y) {
	ck <- .corbay_checks(x, y, 1L); ch <- .corbay_choose(ck)
	lb <- .corbay_lbf(ch$method, x, y, 1/3); h <- .corbay_hyp(lb, "two.sided", 0.5, "x", "y")
	om <- if (ch$method == "kendall") "pearson" else "kendall"
	oh <- .corbay_hyp(.corbay_lbf(om, x, y, 1/3), "two.sided", 0.5, "x", "y")
	expl <- ch$method == "kendall" && any(c("far", "influence") %in% ch$why)
	contra <- !expl && !identical(oh$key[1], h$key[1]) && h$posterior[1] >= 0.75 && oh$posterior[1] >= 0.75
	data.frame(meth = ch$method, why = paste(ch$why, collapse = "+"), other_lead = oh$key[1], other_p = oh$posterior[1], lead = h$key[1], p_lead = h$posterior[1], p_link = 1 - h$posterior[h$key == "none"],
		p_pos = h$posterior[h$key == "pos"], contra = contra, cP = cor.test(x, y)$p.value,
		cK = suppressWarnings(cor.test(x, y, method = "kendall")$p.value))
}
set.seed(31); res <- list()
for (sc in if (is.null(only)) names(gen) else only) for (rho in c(0, 0.3, 0.5)) for (n in c(10, 20, 50, 100)) {
	if (rho > 0 && n == 10) next
	for (r in seq_len(R)) { d <- gen[[sc]](n, rho); while (length(unique(d$x)) < 2 || length(unique(d$y)) < 2) d <- gen[[sc]](n, rho); res[[length(res) + 1]] <- cbind(sc = sc, rho = rho, n = n, one(d$x, d$y)) }
	cat(sc, rho, n, "\n")
}
saveRDS(do.call(rbind, res), paste0("resultats/sim_auto", sfx, ".rds"))

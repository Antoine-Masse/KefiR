# agent_66 : la correction D19 (lois de variance par case moyennees, sys_bay_fact_laws.R)
# rejouee sur les plans d'agent_64 (H0, variances inegales) et d'agent_65 (variances egales,
# effet de A 0,8). Memes graines. Mesures sans aucun plafond : part des jeux ou un terme est a
# P > .95 (agent_64, cible ~0,7 %), puissance P(A) > .95 (agent_65, 17,5 % avant ; perte
# toleree un tiers), P(variance par case), plafond residuel (echec du calcul), temps par appel.
args <- commandArgs(TRUE)
which_ <- if (length(args) >= 1L) args[1L] else "64"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
off <- if (length(args) >= 3L) as.integer(args[3L]) else 0L  # further datasets: seeds after the first N
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
one <- function(d) {
	t0 <- proc.time()[["elapsed"]]
	rr <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	list(rr = rr, t = proc.time()[["elapsed"]] - t0)
}
if (which_ == "64") {
	n <- c(5, 5, 5, 40); sdv <- c(4, 4, 4, 1)
	r <- t(vapply(seq_len(N), function(i) {
		set.seed(640000L + off + i)
		cell <- rep(1:4, n)
		d <- data.frame(A = factor(c("a1", "a2", "a1", "a2")[cell]), B = factor(c("b1", "b1", "b2", "b2")[cell]))
		d$y <- rnorm(length(cell), 0, sdv[cell])
		o <- one(d)
		c(any95 = any(o$rr$effects$p_h1 > 0.95), any75 = any(o$rr$effects$p_h1 > 0.75),
			het = isTRUE(o$rr$het_used), cap = isTRUE(o$rr$het_cells),
			p_het = if (is.null(o$rr$laws)) NA else o$rr$laws$posterior[2L], t = o$t)
	}, numeric(6)))
	out <- c(R.version.string, sprintf("agent_66/64 : 2 x 2, H0, n = (5,5,5,40), sd = (4,4,4,1), N = %d (graines +%d), lois moyennees (D19)", N, off),
		sprintf("un terme P>.95 : %.3f | P>.75 : %.3f | lois moyennees : %.3f | plafond residuel : %.3f | P(var. par case) mediane : %.3f | temps median %.2f s, max %.2f s",
			mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), stats::median(r[, 5], na.rm = TRUE), stats::median(r[, 6]), max(r[, 6])))
} else {
	r <- t(vapply(seq_len(N), function(i) {
		set.seed(650000L + off + i)
		d <- expand.grid(rep = 1:10, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
		d$y <- 0.8 * (d$A == "a2") + rnorm(nrow(d))
		o <- one(d)
		pA <- o$rr$effects$p_h1[o$rr$effects$term == "A"]
		c(a75 = pA > 0.75, a95 = pA > 0.95, other95 = any(o$rr$effects$p_h1[o$rr$effects$term != "A"] > 0.95),
			cap = isTRUE(o$rr$het_cells), p_het = if (is.null(o$rr$laws)) NA else o$rr$laws$posterior[2L], t = o$t)
	}, numeric(6)))
	out <- c(R.version.string, sprintf("agent_66/65 : 2 x 2, variances egales, effet A 0.8, n = 10 par case (N = 40), N = %d (graines +%d), lois moyennees (D19)", N, off),
		sprintf("P(A)>.75 : %.3f | P(A)>.95 : %.3f (avant 0.175) | B ou A:B P>.95 : %.3f | plafond residuel : %.3f | P(var. par case) mediane : %.3f | temps median %.2f s, max %.2f s",
			mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), stats::median(r[, 5], na.rm = TRUE), stats::median(r[, 6]), max(r[, 6])))
}
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_66_%s%s.txt", which_, if (off) paste0("_suite", off) else ""))); cat(out, sep = "\n")

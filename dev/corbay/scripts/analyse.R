# Tables of VALIDATION from the simulations (sim_h0.rds, sim_power.rds, sim_calib.rds).
setwd("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts")
h <- readRDS("resultats/sim_h0.rds")
h$biv[is.na(h$biv)] <- FALSE
h$kend <- h$far | h$biv | h$sc == "likert"
h$pA <- ifelse(h$kend, h$pK, h$pP)
f <- function(v) round(100 * mean(v), 1)
cat("== H0: % of sets with P(link) > 0.75 / > 0.95, and cor.test p < 0.05\n")
t <- do.call(rbind, lapply(split(h, list(h$sc, h$n), drop = TRUE), function(s) data.frame(sc = s$sc[1], n = s$n[1], R = nrow(s),
	kendall_chosen = f(s$kend),
	auto75 = f(s$pA > .75), auto95 = f(s$pA > .95),
	P75 = f(s$pP > .75), P95 = f(s$pP > .95), K75 = f(s$pK > .75), K95 = f(s$pK > .95), S75 = f(s$pS > .75), S95 = f(s$pS > .95),
	pP05 = f(s$cP < .05), pK05 = f(s$cK < .05), pS05 = f(s$cS < .05))))
t <- t[order(t$sc, t$n), ]; rownames(t) <- NULL
print(t)
if (file.exists("resultats/sim_power.rds")) {
	p <- readRDS("resultats/sim_power.rds")
	cat("\n== Power: % of sets with P(positive link) > 0.75 / > 0.95, and cor.test p < 0.05 (two-sided)\n")
	t2 <- do.call(rbind, lapply(split(p, list(p$sc, p$rho, p$n), drop = TRUE), function(s) data.frame(sc = s$sc[1], rho = s$rho[1], n = s$n[1],
		P75 = f(s$pP > .75), P95 = f(s$pP > .95), K75 = f(s$pK > .75), K95 = f(s$pK > .95), S95 = f(s$pS > .95),
		pP05 = f(s$cP < .05), pK05 = f(s$cK < .05), pS05 = f(s$cS < .05))))
	t2 <- t2[order(t2$sc, t2$rho, t2$n), ]; rownames(t2) <- NULL
	print(t2)
}
if (file.exists("resultats/sim_calib.rds")) {
	cl <- readRDS("resultats/sim_calib.rds")
	cat("\n== Calibration under the prior: announced P(link) by class against the frequency of a true link\n")
	for (m in c("pP", "pK", "pS")) for (n in unique(cl$n)) {
		s <- cl[cl$n == n, ]
		b <- cut(s[[m]], c(0, .2, .4, .6, .8, 1), include.lowest = TRUE)
		cat(m, "n =", n, "\n")
		print(data.frame(announced = round(tapply(s[[m]], b, mean), 3), freq_h1 = round(tapply(s$h1, b, mean), 3), count = as.vector(table(b))))
	}
}

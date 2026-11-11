setwd("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts")
a <- readRDS("resultats/sim_auto.rds")
f <- function(v) round(100 * mean(v), 1)
t <- do.call(rbind, lapply(split(a, list(a$sc, a$rho, a$n), drop = TRUE), function(s) {
	right <- if (s$rho[1] == 0) "none" else if (s$sc[1] == "outlier") "pos" else "pos"
	data.frame(sc = s$sc[1], rho = s$rho[1], n = s$n[1], kendall = f(s$meth == "kendall"),
		link75 = f(s$p_link > .75), link95 = f(s$p_link > .95),
		pos75 = f(s$lead == "pos" & s$p_lead > .75), pos95 = f(s$lead == "pos" & s$p_lead > .95),
		none75 = f(s$lead == "none" & s$p_lead > .75), wrongsign75 = f(s$lead == "neg" & s$p_lead > .75),
		contra = f(s$contra), pP05 = f(s$cP < .05), pK05 = f(s$cK < .05))
}))
t <- t[order(t$sc, t$rho, t$n), ]; rownames(t) <- NULL
print(t, row.names = FALSE)

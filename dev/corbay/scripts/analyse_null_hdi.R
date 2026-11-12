setwd("C:/Users/masse/Desktop/KefiR/KefiR/.claude/worktrees/feat-cor-bay-59ccb0/dev/corbay/scripts")
a <- do.call(rbind, lapply(list.files("resultats", "^sim_null_hdi_r0[23][.]rds$", full.names = TRUE), readRDS))
f <- function(v) round(100 * mean(v), 1)
ex <- function(lo, hi) lo > 0 | hi < 0
t <- do.call(rbind, lapply(split(a, list(a$rho, a$n)), function(s) {
	none <- s$lead == "none"
	data.frame(rho = s$rho[1], n = s$n[1], N = nrow(s), none_lead = f(none),
		none_mod = f(none & s$level >= 3), none_strong = f(none & s$level >= 4),
		hdi_excl0 = f(ex(s$lo, s$hi)),
		none_mod_excl0 = f(none & s$level >= 3 & ex(s$lo, s$hi)), none_strong_excl0 = f(none & s$level >= 4 & ex(s$lo, s$hi)),
		none_lead_excl0 = f(none & ex(s$lo, s$hi)), none_lead_excl0_testprior = f(none & ex(s$lo_t, s$hi_t)),
		max_p_none = round(max(c(0, s$p_lead[none])), 3))
}))
print(t[order(t$rho, t$n), ], row.names = FALSE)

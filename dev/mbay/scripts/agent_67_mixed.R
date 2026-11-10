# agent_67 : modele mixte avec repliques (D21). 10 sujets x 3 conditions x 4 repliques ;
# ecarts-types sujet 1, sujet x condition 0,7, residuel 1. Scenario H0 (aucun effet) et
# scenario d'effet (0 ; 0,5 ; 0,6). Part des jeux a P > .75 et P > .95 pour m.bay (sujet et
# sujet x condition aleatoires), pour le modele naif (sujet seul, lmBF, les repliques comptees
# comme independantes), et p < .05 pour l'ANOVA a mesures repetees sur les moyennes par sujet
# et condition (analyse classique correcte). Temps par appel.
args <- commandArgs(TRUE)
which_ <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
eff <- if (which_ == "H0") c(0, 0, 0) else c(0, 0.5, 0.6)
r <- t(vapply(seq_len(N), function(i) {
	set.seed(670000L + i + if (which_ == "H0") 0L else 50000L)
	d <- expand.grid(rep = 1:4, cond = factor(c("c1", "c2", "c3")), id = factor(1:10))
	u0 <- stats::rnorm(10); uc <- matrix(stats::rnorm(30, 0, 0.7), 10)
	d$y <- u0[d$id] + uc[cbind(as.integer(d$id), as.integer(d$cond))] + eff[d$cond] + stats::rnorm(nrow(d))
	t0 <- proc.time()[["elapsed"]]
	rr <- suppressMessages(m.bay(y ~ cond | id, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	tt <- proc.time()[["elapsed"]] - t0
	p <- rr$effects$p_h1[1L]
	nv <- BayesFactor::extractBF(BayesFactor::lmBF(y ~ cond + id, data = d, whichRandom = "id", progress = FALSE) /
		BayesFactor::lmBF(y ~ id, data = d, whichRandom = "id", progress = FALSE), logbf = TRUE)$bf
	pn <- 1 / (1 + exp(-nv))
	m <- stats::aggregate(y ~ cond + id, data = d, FUN = mean)
	pa <- summary(stats::aov(y ~ cond + Error(id / cond), data = m))[["Error: id:cond"]][[1]][1L, "Pr(>F)"]
	c(p75 = p > 0.75, p95 = p > 0.95, n75 = pn > 0.75, n95 = pn > 0.95, aov = pa < 0.05, t = tt)
}, numeric(6)))
out <- c(R.version.string, sprintf("agent_67 : mixte avec repliques, scenario %s, 10 sujets x 3 conditions x 4 repliques, N = %d", which_, N),
	sprintf("m.bay P>.75 : %.3f | P>.95 : %.3f || naif (sujet seul) P>.75 : %.3f | P>.95 : %.3f || ANOVA RM sur moyennes p<.05 : %.3f || temps median %.2f s, max %.2f s",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5]), stats::median(r[, 6]), max(r[, 6])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_67_%s.txt", which_))); cat(out, sep = "\n")

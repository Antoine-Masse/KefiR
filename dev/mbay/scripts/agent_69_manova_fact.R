# agent_69 : MANOVA a deux facteurs (D23), 2 x 2, 3 reponses (correlations 0,3), n = 12 par
# case. Scenarios : H0, et effet de A de 0,8 sur une reponse. Part des jeux ou un terme est
# a P > .75 / > .95 (m.bay, sans sensibilite), et Pillai p < .05 sur ce terme.
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 60L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(690000L + i + if (sc == "H0") 0L else 5000L)
	d <- expand.grid(rep = 1:12, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	R <- matrix(0.3, 3, 3); diag(R) <- 1
	Y <- matrix(stats::rnorm(nrow(d) * 3), ncol = 3) %*% chol(R)
	if (sc != "H0") Y[d$A == "a2", 1] <- Y[d$A == "a2", 1] + 0.8
	d$y1 <- Y[, 1]; d$y2 <- Y[, 2]; d$y3 <- Y[, 3]
	t0 <- proc.time()[["elapsed"]]
	rr <- suppressMessages(m.bay(cbind(y1, y2, y3) ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	tt <- proc.time()[["elapsed"]] - t0
	pil <- summary(stats::manova(cbind(y1, y2, y3) ~ A * B, data = d))$stats[1:3, "Pr(>F)"]
	p <- rr$effects$p_h1
	c(any75 = any(p > 0.75), any95 = any(p > 0.95), A75 = p[1] > 0.75, A95 = p[1] > 0.95, other95 = any(p[2:3] > 0.95),
		pilA = pil[1] < 0.05, pil_any = any(pil < 0.05), t = tt)
}, numeric(8)))
out <- c(R.version.string, sprintf("agent_69 : MANOVA 2 x 2, 3 reponses, n = 12 par case, scenario %s, N = %d", sc, N),
	sprintf("un terme P>.75 : %.3f | P>.95 : %.3f || A P>.75 : %.3f | A P>.95 : %.3f | B ou A:B P>.95 : %.3f || Pillai A p<.05 : %.3f | un terme p<.05 : %.3f || temps median %.1f s",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5]), mean(r[, 6]), mean(r[, 7]), stats::median(r[, 8])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_69_%s.txt", sc))); cat(out, sep = "\n")

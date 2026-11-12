# agent_75 : reponse ordonnee a plusieurs facteurs (D26), 2 x 2, n = 25 par case, notes de
# 1 a 5 (normale latente arrondie). Scenarios : H0, et effet de A de 0,5 ecart-type latent.
# Part des jeux ou un terme est a P > .75 / .95, P(A) > .75 / .95, et polr (rapport de
# vraisemblance sur A, modele additif) a p < .05.
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(750000L + i + if (sc == "H0") 0L else 10000L)
	d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$y <- pmin(5, pmax(1, round(3 + (if (sc == "A") 0.5 else 0) * (d$A == "a2") + stats::rnorm(nrow(d)))))
	rr <- suppressMessages(suppressWarnings(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, response = "ordinal")))
	p <- rr$effects$p_h1
	o <- factor(d$y, ordered = TRUE)
	l1 <- stats::logLik(MASS::polr(o ~ A + B, data = d)); l0 <- stats::logLik(MASS::polr(o ~ B, data = d))
	c(any75 = any(p > 0.75), any95 = any(p > 0.95), A75 = p[1L] > 0.75, A95 = p[1L] > 0.95,
		lrt = stats::pchisq(2 * (as.numeric(l1) - as.numeric(l0)), 1, lower.tail = FALSE) < 0.05)
}, numeric(5)))
out <- c(R.version.string, sprintf("agent_75 : reponse ordonnee a plusieurs facteurs, 2 x 2, n = 25, notes 1-5, scenario %s, N = %d", sc, N),
	sprintf("un terme P>.75 : %.3f | P>.95 : %.3f || P(A)>.75 : %.3f | P(A)>.95 : %.3f || polr, rapport de vraisemblance sur A p<.05 : %.3f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_75_%s.txt", sc))); cat(out, sep = "\n")

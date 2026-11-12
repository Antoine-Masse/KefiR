# agent_73 : reponses discretes a plusieurs facteurs (D25), 2 x 2, n = 25 par case.
#   bin_H0 / bin_A : oui/non, p = 0,4 partout, ou 0,4 -> 0,65 selon A.
#   cnt_H0 / cnt_A : comptages surdisperses (binomiale negative, taille 2), moyenne 4, ou
#                    4 -> 6 selon A.
# Part des jeux ou un terme est a P > .75 / .95, P(A) > .75 / .95, et glm (test du rapport de
# vraisemblance sur A, modele additif) a p < .05.
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "bin_H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(730000L + i + match(sc, c("bin_H0", "bin_A", "cnt_H0", "cnt_A")) * 10000L)
	d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	eff <- grepl("_A", sc) * (d$A == "a2")
	if (startsWith(sc, "bin")) {
		d$y <- factor(ifelse(stats::rbinom(nrow(d), 1, 0.4 + 0.25 * eff) == 1, "oui", "non"), levels = c("non", "oui"))
		g1 <- stats::glm(I(y == "oui") ~ A + B, data = d, family = stats::binomial()); g0 <- stats::glm(I(y == "oui") ~ B, data = d, family = stats::binomial())
	} else {
		d$y <- stats::rnbinom(nrow(d), mu = 4 + 2 * eff, size = 2)
		g1 <- MASS::glm.nb(y ~ A + B, data = d); g0 <- MASS::glm.nb(y ~ B, data = d)
	}
	# Counts: response = "counts", so that the engine is measured and not the routing of D17
	# (a count with few values and a maximum at 10 can be read as a rating scale).
	rr <- suppressMessages(suppressWarnings(if (startsWith(sc, "cnt")) m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, response = "counts") else
		m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE)))
	p <- rr$effects$p_h1
	lrt <- stats::pchisq(2 * (as.numeric(stats::logLik(g1)) - as.numeric(stats::logLik(g0))), 1, lower.tail = FALSE)
	c(any75 = any(p > 0.75), any95 = any(p > 0.95), A75 = p[1L] > 0.75, A95 = p[1L] > 0.95, lrt = lrt < 0.05)
}, numeric(5)))
out <- c(R.version.string, sprintf("agent_73 : discret a plusieurs facteurs, 2 x 2, n = 25, scenario %s, N = %d", sc, N),
	sprintf("un terme P>.75 : %.3f | P>.95 : %.3f || P(A)>.75 : %.3f | P(A)>.95 : %.3f || glm, rapport de vraisemblance sur A p<.05 : %.3f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), mean(r[, 5])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_73_%s%s.txt", sc, if (startsWith(sc, "cnt")) "_rc" else ""))); cat(out, sep = "\n")

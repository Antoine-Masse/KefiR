# agent_79 : D25 (BIC, information unitaire) sous un effet modere, demande de la session
# parente apres l'audit de D30 : le BIC dit-il « non, soutien fort » quand le facteur agit ?
#   bin_or2 : oui/non 2 x 2, n = 25 par case, rapport de cotes de 2 sur A (p = 0,4 -> 0,571)
#   cnt_rr15: comptages 2 x 2, binomiale negative taille 2, moyenne 4 -> 6 selon A (rapport 1,5)
# Mesures sur A : P(A) <= .25 (« non »), facteur de Bayes de A <= 1/20 (« non, soutien fort »),
# <= 1/3 (au moins modere), P(A) > .75, et test du rapport de vraisemblance p < .05.
args <- commandArgs(TRUE)
sc <- args[1L]
N <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(790000L + i + if (sc == "bin_or2") 0L else 50000L)
	d <- expand.grid(rep = 1:25, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	if (sc == "bin_or2") {
		p <- ifelse(d$A == "a2", 0.4 * 2 / (0.6 + 0.4 * 2), 0.4)
		d$y <- factor(ifelse(stats::rbinom(nrow(d), 1, p) == 1, "oui", "non"), levels = c("non", "oui"))
		g1 <- stats::glm(I(y == "oui") ~ A + B, data = d, family = stats::binomial()); g0 <- stats::glm(I(y == "oui") ~ B, data = d, family = stats::binomial())
		rr <- suppressMessages(suppressWarnings(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE)))
	} else {
		d$y <- stats::rnbinom(nrow(d), mu = ifelse(d$A == "a2", 6, 4), size = 2)
		g1 <- MASS::glm.nb(y ~ A + B, data = d); g0 <- MASS::glm.nb(y ~ B, data = d)
		rr <- suppressMessages(suppressWarnings(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, response = "counts")))
	}
	pA <- rr$effects$p_h1[1L]; lA <- rr$effects$log_bf[1L]
	lrt <- stats::pchisq(2 * (as.numeric(stats::logLik(g1)) - as.numeric(stats::logLik(g0))), 1, lower.tail = FALSE)
	c(non = pA <= 0.25, strong_non = lA <= -log(20), mod_non = lA <= -log(3), oui = pA > 0.75, lrt = lrt < 0.05)
}, numeric(5)))
wil <- function(k, n) { z <- 1.96; ph <- k / n
	c((ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
f <- function(j) { k <- sum(r[, j]); w <- wil(k, N); sprintf("%d/%d = %.3f [%.3f ; %.3f]", k, N, k / N, w[1], w[2]) }
out <- c(R.version.string, sprintf("agent_79 : D25 sous effet modere, scenario %s, N = %d", sc, N),
	paste0("P(A)<=.25 (non) : ", f(1), " | BF(A)<=1/20 (non, soutien fort) : ", f(2), " | BF(A)<=1/3 : ", f(3)),
	paste0("P(A)>.75 (oui) : ", f(4), " | rapport de vraisemblance p<.05 : ", f(5)))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_79_%s.txt", sc))); cat(out, sep = "\n")

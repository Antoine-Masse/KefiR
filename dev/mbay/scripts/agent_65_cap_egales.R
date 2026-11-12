# agent_65 : la prudence "variances par case" sur un 2 x 2 a variances egales (n = 10 par case,
# effet de A de 0,8 ecart-type) : part des jeux ou elle se declenche a tort, et puissance sur A
# (P(A) > .75, > .95) avec et sans le plafond du qualificatif.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
r <- t(vapply(seq_len(N), function(i) {
	set.seed(650000L + i)
	d <- expand.grid(rep = 1:10, A = factor(c("a1", "a2")), B = factor(c("b1", "b2")))
	d$y <- 0.8 * (d$A == "a2") + rnorm(nrow(d))
	rr <- suppressMessages(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE, draws = 2000))
	pA <- rr$effects$p_h1[rr$effects$term == "A"]
	c(cap = isTRUE(rr$het_cells), a75 = pA > 0.75, a95 = pA > 0.95, a95_affiche = pA > 0.95 && !isTRUE(rr$het_cells))
}, numeric(4)))
out <- c(R.version.string, sprintf("agent_65 : 2 x 2, variances egales, effet A 0.8, n = 10 par case, N = %d", N),
	sprintf("prudence declenchee a tort : %.3f | P(A)>.75 : %.3f | P(A)>.95 : %.3f | oui fort affiche sur A : %.3f", mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4])))
writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_65.txt")); cat(out, sep = "\n")

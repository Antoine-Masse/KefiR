# agent_74 : la symetrie composee supposee en mesures repetees (D15), mesuree contre une
# sphericite violee. 15 sujets x 4 temps ; covariance AR(1), rho = 0,8 (temps voisins tres
# correles, eloignes peu), ecart-type 1. Scenarios : H0, et effet sur le dernier temps (0,5),
# et H0 avec symetrie composee (rho = 0,5 partout), temoin. Part des jeux ou P(temps) > .75 /
# .95 (m.bay), p < .05 en ANOVA a mesures repetees non corrigee et corrigee de
# Greenhouse-Geisser (epsilon calcule a la main), epsilon median.
args <- commandArgs(TRUE)
sc <- if (length(args) >= 1L) args[1L] else "H0"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
k <- 4L; n <- 15L
Sig <- if (sc == "H0_cs") { m <- matrix(0.5, k, k); diag(m) <- 1; m } else 0.8^abs(outer(1:k, 1:k, "-"))
mu <- if (sc == "A") c(0, 0, 0, 0.5) else rep(0, k)
C <- stats::contr.helmert(k); C <- sweep(C, 2L, sqrt(colSums(C^2)), "/")
r <- t(vapply(seq_len(N), function(i) {
	set.seed(740000L + i + match(sc, c("H0", "A", "H0_cs")) * 10000L)
	Y <- matrix(stats::rnorm(n * k), n, k) %*% chol(Sig) + rep(mu, each = n)
	d <- data.frame(y = as.vector(Y), t = factor(rep(paste0("T", 1:k), each = n)), id = factor(rep(1:n, k)))
	rr <- suppressMessages(m.bay(y ~ t | id, data = d, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, draws = 2000))
	p <- rr$effects$p_h1[1L]
	# Uncorrected repeated-measures F and Greenhouse-Geisser.
	Z <- Y %*% C
	ssb <- n * sum((colMeans(Y) - mean(Y))^2)
	sse <- sum((Y - outer(rowMeans(Y), rep(1, k)) - outer(rep(1, n), colMeans(Y)) + mean(Y))^2)
	Fv <- (ssb / (k - 1)) / (sse / ((k - 1) * (n - 1)))
	S <- stats::cov(Z)
	eps <- sum(diag(S))^2 / ((k - 1) * sum(S^2))
	c(p75 = p > 0.75, p95 = p > 0.95, raw = stats::pf(Fv, k - 1, (k - 1) * (n - 1), lower.tail = FALSE) < 0.05,
		gg = stats::pf(Fv, eps * (k - 1), eps * (k - 1) * (n - 1), lower.tail = FALSE) < 0.05, eps = eps)
}, numeric(5)))
out <- c(R.version.string, sprintf("agent_74 : sphericite, 15 sujets x 4 temps, scenario %s, N = %d", sc, N),
	sprintf("m.bay P(temps)>.75 : %.3f | >.95 : %.3f || ANOVA MR non corrigee p<.05 : %.3f | Greenhouse-Geisser p<.05 : %.3f || epsilon median %.2f",
		mean(r[, 1]), mean(r[, 2]), mean(r[, 3]), mean(r[, 4]), stats::median(r[, 5])))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_74_%s.txt", sc))); cat(out, sep = "\n")

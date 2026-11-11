# agent_71 : calibration du post-hoc au-dela de 8 groupes (echantillonneur de Gibbs sur les
# regroupements, D11). Trois parties.
#   accord : k = 8, calcul exact contre Gibbs force (option KefiR.bay_exact_max = 7), meme jeu ;
#            ecart maximal des P(differents) par paire.
#   H0     : k = 10 et 12, n = 8 par groupe, aucun effet ; part des paires dites « differentes »
#            (P >= 0,75) et part des jeux avec au moins une telle paire.
#   paquets: k = 10, deux paquets de 5 groupes ecartes de 1 ecart-type ; part des paires bien
#            classees (differentes entre paquets, semblables dedans).
args <- commandArgs(TRUE)
part <- if (length(args) >= 1L) args[1L] else "accord"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 30L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
mk <- function(k, n, mu, seed) {
	set.seed(seed)
	g <- factor(rep(sprintf("g%02d", seq_len(k)), each = n))
	list(y = mu[as.integer(g)] + stats::rnorm(k * n), g = g)
}
out <- c(R.version.string)
if (part == "accord") {
	dev <- vapply(seq_len(N), function(i) {
		d <- mk(8, 8, c(0, 0, 0, 0.8, 0.8, 0.8, 1.6, 1.6), 710000L + i)
		ex <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, i)
		old <- options(KefiR.bay_exact_max = 7L); on.exit(options(old))
		gb <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, i)
		options(old)
		c(max(abs(ex$pairs$p_differ - gb$pairs$p_differ)), gb$mc_error, identical(gb$method, "gibbs"))
	}, numeric(3))
	out <- c(out, sprintf("agent_71 accord : k = 8, N = %d ; ecart max |P exact - P Gibbs| : mediane %.3f, max %.3f ; erreur Monte-Carlo annoncee mediane %.3f ; Gibbs utilise : %.2f",
		N, stats::median(dev[1, ]), max(dev[1, ]), stats::median(dev[2, ]), mean(dev[3, ])))
} else if (part == "H0") {
	for (k in c(10L, 12L)) {
		r <- vapply(seq_len(N), function(i) {
			d <- mk(k, 8, rep(0, k), 711000L + 100L * k + i)
			ph <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, i)
			c(mean(ph$pairs$p_differ >= 0.75), any(ph$pairs$p_differ >= 0.75), any(ph$pairs$p_differ >= 0.95))
		}, numeric(3))
		out <- c(out, sprintf("agent_71 H0 : k = %d, n = 8, N = %d ; part des paires P>=.75 : %.4f ; jeux avec une paire P>=.75 : %.3f ; P>=.95 : %.3f",
			k, N, mean(r[1, ]), mean(r[2, ]), mean(r[3, ])))
	}
} else {
	r <- vapply(seq_len(N), function(i) {
		d <- mk(10, 8, rep(c(0, 1), each = 5), 712000L + i)
		ph <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, i)
		pr <- ph$pairs
		between <- (as.integer(substr(pr$group1, 2, 3)) <= 5) != (as.integer(substr(pr$group2, 2, 3)) <= 5)
		c(mean(pr$p_differ[between] >= 0.75), mean(pr$p_differ[!between] >= 0.75), mean(pr$p_differ[!between] <= 0.25))
	}, numeric(3))
	out <- c(out, sprintf("agent_71 paquets : k = 10, deux paquets a 1 ecart-type, n = 8, N = %d ; paires entre paquets dites differentes : %.3f ; paires d'un meme paquet dites differentes : %.3f, semblables : %.3f",
		N, mean(r[1, ]), mean(r[2, ]), mean(r[3, ])))
}
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_71_%s.txt", part))); cat(out, sep = "\n")

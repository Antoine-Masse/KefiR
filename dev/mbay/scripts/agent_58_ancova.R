# agent_58 : calibration de l'analyse de covariance (D16).
# Par scenario (3 groupes de n, covariable U(5, 15), pente 0.8, residu sd 1) : part des jeux
# ou m.bay dit "difference selon g" (P > .5, > .75, > .95), a cote de car::Anova type II
# p < .05 ; couverture des intervalles d'ecart ajuste (paires) ; avec interaction ecrite
# (y ~ g * x) : part des jeux ou des pentes differentes sont dites probables (P(g:x) > .5),
# a cote de p(g:x) < .05.
# Usage : Rscript agent_58_ancova.R N coeurs
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	suppressMessages(devtools::load_all(BASE, quiet = TRUE))
	NULL
}, LIB, BASE))
SCEN <- list(
	list(nom = "paral_H0_n15", eff = c(0, 0, 0), sl = 0, n = 15, inter = FALSE),
	list(nom = "paral_un_ecart_n15", eff = c(0, 0, 1), sl = 0, n = 15, inter = FALSE),
	list(nom = "inter_ecrite_pentes_egales_n20", eff = c(0, 0, 0.5), sl = 0, n = 20, inter = TRUE),
	list(nom = "inter_ecrite_pente_differente_n20", eff = c(0, 0, 0.5), sl = 0.4, n = 20, inter = TRUE))
un <- function(i, s) {
	set.seed(580000L + 37L * nchar(s$nom) + 11L * s$n + i)
	d <- data.frame(g = factor(rep(c("a", "b", "c"), each = s$n)), x = runif(3 * s$n, 5, 15))
	d$y <- 10 + 0.8 * d$x + s$eff[as.integer(d$g)] + s$sl * (d$g == "c") * (d$x - 10) + rnorm(nrow(d))
	fo <- if (s$inter) y ~ g * x else y ~ g + x
	r <- suppressMessages(m.bay(fo, data = d, verbose = FALSE, draws = 2000L, sensitivity = FALSE, seed = i))
	av <- car::Anova(stats::lm(fo, data = d), type = 2)
	pg <- r$effects$p_h1[r$effects$term == "g"]
	pint <- if (s$inter) r$effects$p_h1[r$effects$term == "g:x"] else NA
	cover <- NA
	if (identical(r$posthoc$method, "exact")) {
		pr <- r$posthoc$pairs
		tv <- s$eff[match(pr$group2, c("a", "b", "c"))] - s$eff[match(pr$group1, c("a", "b", "c"))]
		cover <- mean(pr$diff_low <= tv & tv <= pr$diff_high)
	}
	c(pg = pg, pcl = av["g", "Pr(>F)"], pint = pint, pcl_int = if (s$inter) av["g:x", "Pr(>F)"] else NA, cover = cover)
}
clusterExport(cl, "un")
out <- c(R.version.string, sprintf("agent_58 : ANCOVA, N = %d jeux par scenario, 3 groupes, covariable U(5,15), pente 0.8, residu 1, r = 0.5, rscaleCont = sqrt(2)/4", N),
	"scenario : P(g)>.5 / >.75 / >.95 | car type II p<.05 || P(g:x)>.5 / >.75 | p(g:x)<.05 || couverture ecart ajuste")
for (s in SCEN) {
	r <- do.call(rbind, parLapplyLB(cl, seq_len(N), un, s = s))
	f <- function(x) if (all(is.na(x))) "-" else sprintf("%.3f", mean(x, na.rm = TRUE))
	out <- c(out, sprintf("%s : %s / %s / %s | %s || %s / %s | %s || %s", s$nom, f(r[, "pg"] > .5), f(r[, "pg"] > .75), f(r[, "pg"] > .95),
		f(r[, "pcl"] < .05), f(r[, "pint"] > .5), f(r[, "pint"] > .75), f(r[, "pcl_int"] < .05), f(r[, "cover"])))
	writeLines(out, file.path(BASE, "dev/mbay/scripts/resultats/agent_58.txt"))
}
stopCluster(cl)

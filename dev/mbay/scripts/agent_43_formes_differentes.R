# agent_43 : groupes de formes differentes, moyennes egales. Sous une loi t, m.bay compare
# des centres : un groupe asymetrique a un centre different de sa moyenne. Annonce-t-il
# un effet que les moyennes n'ont pas ? Reprise possible : resultats/agent_43_cas/.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_43_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(BASE)
	NULL
}, LIB, BASE))
lnorm_std <- function(n, s = 0.8) (rlnorm(n, 0, s) - exp(s^2 / 2)) / sqrt((exp(s^2) - 1) * exp(s^2))
gen <- function(i, n, k) {
	set.seed(i * 2741L + n + k)
	# Groupe 1 asymetrique ; les autres normaux ; memes moyennes et memes variances.
	y <- c(lnorm_std(n), rnorm(n * (k - 1L)))
	g <- factor(rep(seq_len(k), each = n))
	st <- .bay_laws_stats(y, g)
	res <- .bay_laws_average(.bay_laws_table(st, 0.5, 1.5, 1500L, 1L), 0.5)
	c(p = res$p_h1, pg = plogis(.bay_bf_oneway_zs(y, g, 0.5)), welch = oneway.test(y ~ g)$p.value,
		kw = kruskal.test(y ~ g)$p.value)
}
clusterExport(cl, c("gen", "lnorm_std"))
lignes <- "Plan | m.bay P>.5 P>.95 mediane | gaussien P>.5 P>.95 | Welch p<.05 p<.005 | KW p<.05"
for (cs in list(c(30, 2), c(100, 2), c(400, 2), c(30, 3), c(100, 3))) {
	r <- NULL
	for (d in seq(1L, N, by = 50L)) {
		fic <- file.path(CAS, sprintf("n%d_k%d_%05d.rds", cs[1], cs[2], d))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d:min(d + 49L, N), gen, n = cs[1], k = cs[2])), fic)
		r <- rbind(r, readRDS(fic))
	}
	lignes <- c(lignes, sprintf("n=%d k=%d | %.3f %.3f %.3f | %.3f %.3f | %.3f %.3f | %.3f", cs[1], cs[2], mean(r[, "p"] > .5), mean(r[, "p"] > .95),
		median(r[, "p"]), mean(r[, "pg"] > .5), mean(r[, "pg"] > .95), mean(r[, "welch"] < .05), mean(r[, "welch"] < .005), mean(r[, "kw"] < .05)))
	writeLines(c(R.version.string, sprintf("N = %d par plan ; groupe 1 lognormal standardise, autres normaux, moyennes et variances egales", N), lignes),
		file.path(BASE, "dev/mbay/scripts/resultats/agent_43.txt"))
}
stopCluster(cl)

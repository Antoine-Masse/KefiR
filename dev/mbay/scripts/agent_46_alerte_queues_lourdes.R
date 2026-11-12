# agent_46 : a partir de quelle probabilite des queues lourdes le bilan peut-il le dire ?
# Taux de fausse alerte sur donnees normales et taux de detection sur t3, n=10, 30, 100.
# Et frequence de la contradiction normale / queues lourdes (points 1 et 4 de l'audit).
# Reprise possible : resultats/agent_46_cas/.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_46_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(BASE)
	NULL
}, LIB, BASE))
gen <- function(i, n, k, type, d) {
	set.seed(i * 8831L + n * 3L + k + nchar(type) + round(10 * d))
	e <- function(m) switch(type, normal = rnorm(m), t3 = rt(m, 3) / sqrt(3),
		lognormal1 = rnorm(m))
	y <- unlist(lapply(seq_len(k), function(j) e(n)))
	if (type == "lognormal1") y[seq_len(n)] <- (rlnorm(n, 0, 0.8) - exp(0.32)) / sqrt((exp(0.64) - 1) * exp(0.64))
	y[(n * (k - 1L) + 1L):(n * k)] <- y[(n * (k - 1L) + 1L):(n * k)] + d
	g <- factor(rep(seq_len(k), each = n))
	a <- .bay_laws_average(.bay_laws_table(.bay_laws_stats(y, g), 0.5, 1.5, 1500L, 1L), 0.5)
	pl <- a$p_law[.bay_law_names()]
	lb <- a$log_bf_law[.bay_law_names()]
	nrm <- c("normal_eq", "normal_het")
	best <- nrm[which.max(replace(pl[nrm], is.na(pl[nrm]), -1))]
	c(p = a$p_h1, p_t = sum(pl[c("t_eq", "t_het")], na.rm = TRUE), p_best_normal = plogis(lb[best]))
}
clusterExport(cl, "gen")
lignes <- "type k n d | P(t) mediane q90 q95 | P(t)>.75 >.9 >.95 >.99 | contradiction, m.bay>.75 normale<.5 ou m.bay<.25 normale>.5 | m.bay>.75"
for (cs in list(list("normal", 0), list("normal", .5), list("t3", 0), list("t3", .5), list("lognormal1", 0))) {
	for (k in c(2L, 3L)) for (n in c(10L, 30L, 100L)) {
		r <- NULL
		for (d0 in seq(1L, N, by = 50L)) {
			fic <- file.path(CAS, sprintf("%s_%d_k%d_n%d_%05d.rds", cs[[1]], round(10 * cs[[2]]), k, n, d0))
			if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d0:min(d0 + 49L, N), gen, n = n, k = k, type = cs[[1]], d = cs[[2]])), fic)
			r <- rbind(r, readRDS(fic))
		}
		pt <- r[, "p_t"]
		contra <- (r[, "p"] > .75 & r[, 3L] < .5) | (r[, "p"] < .25 & r[, 3L] > .5)
		lignes <- c(lignes, sprintf("%s k=%d n=%d d=%.1f | %.3f %.3f %.3f | %.3f %.3f %.3f %.3f | %.3f | %.3f", cs[[1]], k, n, cs[[2]],
			median(pt), quantile(pt, .9), quantile(pt, .95), mean(pt > .75), mean(pt > .9), mean(pt > .95), mean(pt > .99), mean(contra), mean(r[, "p"] > .75)))
		writeLines(c(R.version.string, sprintf("N = %d par plan, 1500 tirages", N), lignes),
			file.path(BASE, "dev/mbay/scripts/resultats/agent_46.txt"))
	}
}
stopCluster(cl)

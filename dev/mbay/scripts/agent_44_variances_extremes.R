# agent_44 : cas signale en relecture (I1) : k=8, n=6, sans effet, ecarts-types de 0,01 a
# 100, P(effet) a 0,998 sur un jeu. Biais du calcul, hasard, ou comportement du modele ?
# Plusieurs jeux sous H0, le meme jeu sous plusieurs graines, et un plan moins extreme.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 100L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R")); .kefir_bay(BASE)
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_44_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
out <- file.path(BASE, "dev/mbay/scripts/resultats/agent_44.txt")
lignes <- c(R.version.string, "Le jeu du relecteur, par graine de calcul :")
set.seed(32); g0 <- rep(sprintf("g%02d", 1:8), each = 6)
y0 <- unlist(lapply(10^seq(-2, 2, length.out = 8), function(s) rnorm(6, 0, s)))
st0 <- .bay_laws_stats(y0, g0)
for (sd in 1:4) {
	a <- .bay_laws_average(.bay_laws_table(st0, 0.5, 1.5, 4000L, sd), 0.5)
	lignes <- c(lignes, sprintf("graine %d : log BF %.2f, P %.4f, err %.1f%%, lois %s, log BF par loi %s", sd, a$log_bf, a$p_h1, a$error_pct,
		paste(sprintf("%s=%.2f", names(a$p_law), a$p_law), collapse = " "), paste(sprintf("%.2f", a$log_bf_law), collapse = " ")))
	writeLines(lignes, out)
}
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R")); .kefir_bay(BASE); NULL
}, LIB, BASE))
gen <- function(i, k, n, lim) {
	set.seed(i * 3571L + k + n)
	sds <- 10^seq(-lim, lim, length.out = k)
	y <- unlist(lapply(sds, function(s) rnorm(n, 0, s)))
	g <- factor(rep(seq_len(k), each = n))
	a <- .bay_laws_average(.bay_laws_table(.bay_laws_stats(y, g), 0.5, 1.5, 1500L, 1L), 0.5)
	pl <- a$p_law[.bay_law_names()]
	c(p = a$p_h1, pl, lb = a$log_bf_law[.bay_law_names()], err = a$error_pct, welch = oneway.test(y ~ g)$p.value)
}
clusterExport(cl, "gen")
lignes <- c(lignes, "", "Sous H0, plans a ecarts-types etages : P>.5 P>.95 mediane | Welch p<.05 p<.005 | P(lois) moyennes | lois ecartees | err mediane")
for (cs in list(c(8, 6, 2), c(8, 6, 1), c(3, 6, 2), c(3, 20, 2), c(8, 20, 2))) {
	r <- NULL
	for (d in seq(1L, N, by = 25L)) {
		fic <- file.path(CAS, sprintf("k%d_n%d_l%d_%05d.rds", cs[1], cs[2], cs[3], d))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d:min(d + 24L, N), gen, k = cs[1], n = cs[2], lim = cs[3])), fic)
		r <- rbind(r, readRDS(fic))
	}
	lignes <- c(lignes, sprintf("k=%d n=%d sd 1e-%d a 1e%d : %.3f %.3f %.3f | %.3f %.3f | %s | %d | %.1f", cs[1], cs[2], cs[3], cs[3],
		mean(r[, "p"] > .5), mean(r[, "p"] > .95), median(r[, "p"]), mean(r[, "welch"] < .05), mean(r[, "welch"] < .005),
		paste(round(colMeans(r[, 2:5], na.rm = TRUE), 3), collapse = " "), sum(is.na(r[, 2:5])), median(r[, "err"])))
	writeLines(lignes, out)
}
stopCluster(cl)

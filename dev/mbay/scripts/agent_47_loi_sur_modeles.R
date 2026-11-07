# agent_47 : loi a priori sur les modeles a plusieurs termes, sur les memes jeux.
#   half : chaque terme eligible entre a 1/2 (facteur 1/2, interaction 1/8).
#   beta : p ~ Beta(1,1) integre sur les termes eligibles du modele
#          (Scott & Berger 2010, equations 8 et 9, lues ; hierarchie : Chipman 1996 ;
#          la combinaison est propre a KefiR).
#   unif : uniforme sur les modeles, pour memoire.
# Vraisemblances : anovaBF, whichModels = "all" filtre sur la hierarchie, loi normale a
# variance commune seulement. Plan 2x2x2 sous H0, et interaction reelle 2x2 a n=10, 25.
# Reprise possible : resultats/agent_47_cas/.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
nc <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_47_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB) {
	.libPaths(c(LIB, .libPaths()))
	suppressMessages(library(BayesFactor))
	options(BFprogress = FALSE)
	NULL
}, LIB))

espace <- function(termes) {
	fac <- strsplit(termes, ":", fixed = TRUE)
	par <- lapply(seq_along(termes), function(i) which(vapply(fac, function(f) length(f) < length(fac[[i]]) && all(f %in% fac[[i]]), TRUE)))
	mods <- list(integer(0))
	for (i in order(lengths(fac))) mods <- c(mods, lapply(Filter(function(m) all(par[[i]] %in% m), mods), function(m) c(m, i)))
	e <- vapply(mods, function(m) sum(vapply(seq_along(termes), function(i) all(par[[i]] %in% m), TRUE)), 0)
	s <- lengths(mods)
	lp <- list(half = e * log(1 / 2), beta = lbeta(1 + s, 1 + e - s), unif = rep(-log(length(mods)), length(mods)))
	list(mods = mods, lp = lp, termes = termes)
}
un_jeu <- function(i, plan, n, ef) {
	set.seed(i * 5209L + n * 7L + round(100 * ef) + nchar(plan))
	if (plan == "2x2x2") {
		d <- expand.grid(A = gl(2, 1), B = gl(2, 1), C = gl(2, 1), r = seq_len(n))
		d$y <- rnorm(nrow(d))
		fo <- y ~ A * B * C
		termes <- c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C")
		pv <- summary(aov(fo, d))[[1]][1:7, 5]
	} else {
		d <- expand.grid(A = gl(2, 1), B = gl(2, 1), r = seq_len(n))
		d$y <- rnorm(nrow(d), ef * ifelse(d$A == 1, -1, 1) * ifelse(d$B == 1, -1, 1))
		fo <- y ~ A * B
		termes <- c("A", "B", "A:B")
		pv <- summary(aov(fo, d))[[1]][1:3, 5]
	}
	bf <- anovaBF(fo, data = d, whichModels = "all", progress = FALSE)@bayesFactor
	sp <- espace(termes)
	cle <- function(m) paste(sort(termes[m]), collapse = " + ")
	nm <- vapply(strsplit(rownames(bf), " + ", fixed = TRUE), function(v) paste(sort(trimws(v)), collapse = " + "), "")
	lb <- vapply(sp$mods, function(m) if (!length(m)) 0 else bf$bf[match(cle(m), nm)], 0)
	stopifnot(!anyNA(lb))
	res <- numeric(0)
	for (lo in names(sp$lp)) {
		l <- lb + sp$lp[[lo]]
		p <- exp(l - max(l)); p <- p / sum(p)
		inc <- vapply(seq_along(termes), function(t) sum(p[vapply(sp$mods, function(m) t %in% m, TRUE)]), 0)
		res <- c(res, setNames(c(max(inc), inc[length(termes)], p[lengths(sp$mods) == 0L]), paste0(lo, c("_max", "_haut", "_nul"))))
	}
	c(res, p_min = min(pv), p_haut = pv[length(pv)])
}
clusterExport(cl, c("espace", "un_jeu"))
sp <- espace(c("A", "B", "A:B")); sp3 <- espace(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
inc <- function(s, lo, t) sum(exp(s$lp[[lo]])[vapply(s$mods, function(m) t %in% m, TRUE)])
lignes <- c(R.version.string, sprintf("N = %d par cas", N), "A priori sur les modeles :",
	sprintf("2 facteurs, half : A %.3f, A:B %.3f, nul %.3f ; beta : A %.3f, A:B %.3f, nul %.3f", inc(sp, "half", 1), inc(sp, "half", 3), exp(sp$lp$half[1]), inc(sp, "beta", 1), inc(sp, "beta", 3), exp(sp$lp$beta[1])),
	sprintf("3 facteurs, half : A %.3f, A:B %.3f, A:B:C %.4f, nul %.3f ; beta : A %.3f, A:B %.3f, A:B:C %.4f, nul %.3f", inc(sp3, "half", 1), inc(sp3, "half", 4), inc(sp3, "half", 7), exp(sp3$lp$half[1]), inc(sp3, "beta", 1), inc(sp3, "beta", 4), inc(sp3, "beta", 7), exp(sp3$lp$beta[1])),
	sprintf("sommes : %.6f %.6f %.6f %.6f", sum(exp(sp$lp$half)), sum(exp(sp$lp$beta)), sum(exp(sp3$lp$half)), sum(exp(sp3$lp$beta))), "")
for (cs in list(list("2x2x2", 5, 0), list("2x2x2", 10, 0), list("2x2", 10, 0), list("2x2", 10, .25), list("2x2", 10, .5), list("2x2", 25, 0), list("2x2", 25, .25), list("2x2", 25, .5))) {
	r <- NULL
	for (d0 in seq(1L, N, by = 50L)) {
		fic <- file.path(CAS, sprintf("%s_n%d_f%d_%05d.rds", cs[[1]], cs[[2]], round(100 * cs[[3]]), d0))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d0:min(d0 + 49L, N), un_jeu, plan = cs[[1]], n = cs[[2]], ef = cs[[3]])), fic)
		r <- rbind(r, readRDS(fic))
	}
	li <- sprintf("%s n=%d f=%.2f |", cs[[1]], cs[[2]], cs[[3]])
	for (lo in c("half", "beta", "unif")) {
		li <- paste0(li, sprintf(" %s : un terme>.5 %.3f, >.95 %.3f ; plus haut terme>.5 %.3f, >.95 %.3f ; P(nul) mediane %.3f |", lo,
			mean(r[, paste0(lo, "_max")] > .5), mean(r[, paste0(lo, "_max")] > .95), mean(r[, paste0(lo, "_haut")] > .5), mean(r[, paste0(lo, "_haut")] > .95), median(r[, paste0(lo, "_nul")])))
	}
	li <- paste0(li, sprintf(" aov : un p<.05 %.3f, p<.005 %.3f ; plus haut p<.05 %.3f", mean(r[, "p_min"] < .05), mean(r[, "p_min"] < .005), mean(r[, "p_haut"] < .05)))
	lignes <- c(lignes, li)
	writeLines(lignes, file.path(BASE, "dev/mbay/scripts/resultats/agent_47.txt"))
}
stopCluster(cl)

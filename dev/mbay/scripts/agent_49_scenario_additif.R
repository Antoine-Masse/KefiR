# agent_49 : loi a priori sur les modeles, scenario additif (objection de l'audit sur D2).
# Sachant les parents presents, Beta(1,1) a heredite fait partir l'interaction a 3/4 :
# annonce-t-elle des interactions quand les effets sont seulement additifs ?
# Six lois sur les memes BF (anovaBF, whichModels = "all" filtre sur la hierarchie,
# normale a variance commune) :
#   half       : pi = 1/2 par terme eligible (ancien reglage 1/2, 1/8) ;
#   beta       : pi ~ Beta(1,1) commune a tous les termes eligibles, p(M) = B(1+s, 1+e-s) ;
#   beta_ordre : un pi ~ Beta(1,1) par ordre d'interaction, independants ;
#   hop_ch     : comme beta_ordre, Beta(1, ch_j), ch_j = noeuds d'ordre j eligibles (penalisation de Wilson et al. 2010) ;
#   taille     : uniforme sur la taille, puis uniforme parmi les modeles hereditaires de cette taille ;
#   unif       : uniforme sur les modeles hereditaires.
# Noms de Taylor-Rodriguez, Womack & Bliznyuk (2016, JCGS 25:515, arXiv 1312.6611, sec. 3.1, lu) :
# half = HIP(1,1) = Chipman (1996) ; beta = HUP(1,1) ; beta_ordre = HOP(1,1) ; hop_ch = HOP(1,ch).
# On garde le log BF de chaque modele par jeu : les lois se recalculent sans refaire anovaBF.
# Reprise possible : resultats/agent_49_cas/.
# Usage : Rscript agent_49_scenario_additif.R N_2x2 N_2x2x2 coeurs   (ou "prior" pour les seuls a priori)
args <- commandArgs(TRUE)
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_49_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)

espace <- function(termes) {
	fac <- strsplit(termes, ":", fixed = TRUE)
	ordre <- lengths(fac)
	par <- lapply(seq_along(termes), function(i) which(vapply(fac, function(f) length(f) < length(fac[[i]]) && all(f %in% fac[[i]]), TRUE)))
	mods <- list(integer(0))
	for (i in order(ordre)) mods <- c(mods, lapply(Filter(function(m) all(par[[i]] %in% m), mods), function(m) c(m, i)))
	elig <- lapply(mods, function(m) which(vapply(seq_along(termes), function(i) all(par[[i]] %in% m), TRUE)))
	e <- lengths(elig)
	s <- lengths(mods)
	bo <- vapply(seq_along(mods), function(j) sum(vapply(unique(ordre), function(o) {
		so <- sum(ordre[mods[[j]]] == o)
		lbeta(1 + so, 1 + sum(ordre[elig[[j]]] == o) - so)
	}, 0)), 0)
	# HOP(1, ch) : b_j = nombre de noeuds d'ordre j eligibles (Taylor-Rodriguez et al. 2016, sec. 3.1).
	hc <- vapply(seq_along(mods), function(j) sum(vapply(unique(ordre), function(o) {
		so <- sum(ordre[mods[[j]]] == o); ch <- sum(ordre[elig[[j]]] == o)
		if (ch == 0) 0 else lbeta(1 + so, ch + ch - so) - lbeta(1, ch)
	}, 0)), 0)
	nt <- table(s)
	lp <- list(half = e * log(1 / 2), beta = lbeta(1 + s, 1 + e - s), beta_ordre = bo, hop_ch = hc,
		taille = -log(length(nt)) - log(as.numeric(nt[as.character(s)])),
		unif = rep(-log(length(mods)), length(mods)))
	list(mods = mods, lp = lp, termes = termes)
}
TERMES <- list("2x2" = c("A", "B", "A:B"), "2x2x2" = c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
LOIS <- c("half", "beta", "beta_ordre", "hop_ch", "taille", "unif")

a_priori <- function() {
	lignes <- "A priori : inclusion de chaque terme, et P(A:B | A et B presents)"
	for (pl in names(TERMES)) {
		sp <- espace(TERMES[[pl]])
		ab <- match("A:B", sp$termes)
		ia <- which(vapply(sp$mods, function(m) all(1:2 %in% m), TRUE))
		for (lo in LOIS) {
			pr <- exp(sp$lp[[lo]])
			inc <- vapply(seq_along(sp$termes), function(t) sum(pr[vapply(sp$mods, function(m) t %in% m, TRUE)]), 0)
			lignes <- c(lignes, sprintf("%s %-10s somme %.6f | %s | nul %.4f | P(A:B | A et B) %.3f", pl, lo, sum(pr),
				paste(sprintf("%s %.4f", sp$termes, inc), collapse = ", "), pr[1L],
				sum(pr[ia][vapply(sp$mods[ia], function(m) ab %in% m, TRUE)]) / sum(pr[ia])))
		}
	}
	lignes
}
if (length(args) && args[1L] == "prior") { writeLines(a_priori()); quit(save = "no") }
N2 <- if (length(args) >= 1L) as.integer(args[1L]) else 300L
N3 <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
nc <- if (length(args) >= 3L) as.integer(args[3L]) else 3L

# Scenarios : vrais termes et leur taille f (contrastes +-1, plan equilibre : f de Cohen).
sc <- function(p, n, vrai, f) list(p = p, n = n, vrai = vrai, f = f)
SCEN <- list(
	sc("2x2", 10, character(0), 0), sc("2x2", 25, character(0), 0),
	sc("2x2", 10, c("A", "B"), .25), sc("2x2", 10, c("A", "B"), .5),
	sc("2x2", 25, c("A", "B"), .25), sc("2x2", 25, c("A", "B"), .5),
	sc("2x2", 10, "A:B", .25), sc("2x2", 10, "A:B", .5),
	sc("2x2", 25, "A:B", .25), sc("2x2", 25, "A:B", .5),
	sc("2x2x2", 5, character(0), 0), sc("2x2x2", 10, character(0), 0),
	sc("2x2x2", 5, c("A", "B", "C"), .25), sc("2x2x2", 5, c("A", "B", "C"), .5),
	sc("2x2x2", 10, c("A", "B", "C"), .25), sc("2x2x2", 10, c("A", "B", "C"), .5),
	sc("2x2x2", 5, c("A", "B", "C", "A:B"), .5), sc("2x2x2", 10, c("A", "B", "C", "A:B"), .5),
	sc("2x2x2", 5, c("A", "B", "C", "A:B", "A:C", "B:C"), .5),
	sc("2x2x2", 10, c("A", "B", "C", "A:B", "A:C", "B:C"), .5),
	sc("2x2x2", 10, "A:B:C", .5))
nom <- function(s) sprintf("%s_n%d_%s_f%d", s$p, s$n, if (length(s$vrai)) gsub(":", "", paste(s$vrai, collapse = "-")) else "H0", round(100 * s$f))

un_jeu <- function(i, s) {
	set.seed(i * 7919L + s$n * 13L + round(100 * s$f) + 31L * length(s$vrai) + nchar(paste(s$vrai, collapse = "")) + nchar(s$p))
	termes <- TERMES[[s$p]]
	d <- if (s$p == "2x2") expand.grid(A = gl(2, 1), B = gl(2, 1), r = seq_len(s$n)) else expand.grid(A = gl(2, 1), B = gl(2, 1), C = gl(2, 1), r = seq_len(s$n))
	sg <- function(t) Reduce(`*`, lapply(strsplit(t, ":", fixed = TRUE)[[1]], function(v) ifelse(d[[v]] == 1, -1, 1)))
	mu <- if (length(s$vrai)) s$f * Reduce(`+`, lapply(s$vrai, sg)) else 0
	d$y <- rnorm(nrow(d), mu)
	fo <- if (s$p == "2x2") y ~ A * B else y ~ A * B * C
	pv <- summary(aov(fo, d))[[1]][seq_along(termes), 5]
	bf <- anovaBF(fo, data = d, whichModels = "all", progress = FALSE)@bayesFactor
	sp <- espace(termes)
	cle <- function(m) paste(sort(termes[m]), collapse = " + ")
	nm <- vapply(strsplit(rownames(bf), " + ", fixed = TRUE), function(v) paste(sort(trimws(v)), collapse = " + "), "")
	lb <- vapply(sp$mods, function(m) if (!length(m)) 0 else bf$bf[match(cle(m), nm)], 0)
	stopifnot(!anyNA(lb))
	c(lb, pv)
}
inclusion <- function(lb, sp, lo) {
	l <- lb + sp$lp[[lo]]
	p <- exp(l - max(l)); p <- p / sum(p)
	vapply(seq_along(sp$termes), function(t) sum(p[vapply(sp$mods, function(m) t %in% m, TRUE)]), 0)
}

cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB) {
	.libPaths(c(LIB, .libPaths())); suppressMessages(library(BayesFactor)); options(BFprogress = FALSE); NULL
}, LIB))
clusterExport(cl, c("espace", "un_jeu", "TERMES"))
lignes <- c(R.version.string, a_priori(), "",
	"Faux termes = sans effet reel ; vrais = a effet. Par loi : P(un faux terme > .5), P(un faux terme > .95), P(vrai > .5) moyenne.",
	"aov : P(un faux terme p<.05), P(vrai p<.05) moyenne.")
for (s in SCEN) {
	NN <- if (s$p == "2x2") N2 else N3
	r <- NULL
	for (d0 in seq(1L, NN, by = 50L)) {
		fic <- file.path(CAS, sprintf("%s_%05d.rds", nom(s), d0))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d0:min(d0 + 49L, NN), un_jeu, s = s)), fic)
		r <- rbind(r, readRDS(fic))
	}
	sp <- espace(TERMES[[s$p]]); nm <- length(sp$mods); nt <- length(sp$termes)
	faux <- which(!sp$termes %in% s$vrai); vrais <- which(sp$termes %in% s$vrai)
	li <- sprintf("%s (N=%d) |", nom(s), nrow(r))
	for (lo in LOIS) {
		inc <- t(apply(r[, seq_len(nm), drop = FALSE], 1L, inclusion, sp = sp, lo = lo))
		mf <- if (length(faux)) apply(inc[, faux, drop = FALSE], 1L, max) else rep(0, nrow(r))
		li <- paste0(li, sprintf(" %s %.3f %.3f %s |", lo, mean(mf > .5), mean(mf > .95),
			if (length(vrais)) sprintf("%.3f", mean(inc[, vrais] > .5)) else "-"))
	}
	pv <- r[, nm + seq_len(nt), drop = FALSE]
	li <- paste0(li, sprintf(" aov %.3f %s", mean(apply(pv[, faux, drop = FALSE] < .05, 1L, any)),
		if (length(vrais)) sprintf("%.3f", mean(pv[, vrais] < .05)) else "-"))
	lignes <- c(lignes, li)
	writeLines(lignes, file.path(BASE, "dev/mbay/scripts/resultats/agent_49.txt"))
}
stopCluster(cl)

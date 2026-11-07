# agent_41 : simulations de la voie 2 et k groupes de m.bay (R/sys_bay_laws.R).
#   A. Calibration sous l'a priori, generateur = l'a priori lui-meme, y compris sur
#      le choix de la loi : teste aussi les constantes entre lois.
#   A2. Controle de la tranche (0,8 ; 0,95] de k=2, n=30, autre graine, N plus grand.
#   B. Effets fixes : faux positifs et puissance, a cote d'aov, Welch, Kruskal, et de
#      l'hyper-a priori sur r (moyenne sur r/2, r, 2r).
#   C. Sensibilite aux a priori sur les memes generateurs.
#   D. Defauts B1, B2, I2 de RELECTURE_ADVERSE_2026-11-06.md.
# Usage : Rscript agent_41_lois_simulation.R <partie> <N> [coeurs]
# REPRISE : chaque tranche de 50 jeux est ecrite sur disque des qu'elle finit, dans
# resultats/agent_41_cas/ ; relancer la meme commande reprend ou l'on s'etait arrete.
# La synthese resultats/agent_41_<partie>.txt est reecrite a chaque cas termine.
args <- commandArgs(TRUE)
partie <- if (length(args) >= 1L) args[1L] else "A"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 500L
nc <- if (length(args) >= 3L) as.integer(args[3L]) else 4L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
# Chargement leger, sans les dependances de KefiR (agent_load_bay.R) : 8 Go de memoire.
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
.kefir_bay(Sys.getenv("KEFIR_SRC", BASE))
RES <- file.path(BASE, "dev/mbay/scripts/resultats")
CAS <- file.path(RES, "agent_41_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE, SRC) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(SRC)
	NULL
}, LIB, BASE, Sys.getenv("KEFIR_SRC", BASE)))
out <- file.path(RES, paste0("agent_41_", partie, ".txt"))
tete <- sprintf("%s ; partie %s ; N = %d par case ; %d coeurs ; reprise possible ; derniere ecriture %s\n",
	R.version.string, partie, N, nc, format(Sys.time()))
ecrit <- function(lignes) writeLines(c(tete, lignes), out)

# Calcule N jeux d'un cas par tranches, chacune sur disque ; reprend ce qui existe.
par_tranches <- function(cle, N, f, ...) {
	debut <- seq(1L, N, by = 50L)
	res <- list()
	for (d in debut) {
		fic <- file.path(CAS, sprintf("%s_%05d.rds", cle, d))
		if (file.exists(fic)) {
			res[[length(res) + 1L]] <- readRDS(fic)
			next
		}
		x <- do.call(rbind, parLapplyLB(cl, d:min(d + 49L, N), f, ...))
		saveRDS(x, fic)
		res[[length(res) + 1L]] <- x
	}
	do.call(rbind, res)
}

# Un appel : r = 1/4, 1/2, 1 calcules ensemble. p : r fixe a 1/2. p_mix : hyper-a priori
# uniforme sur ces trois valeurs, qui n'est qu'un a priori plus large sur l'effet : les
# vraisemblances marginales de H1 sont moyennees, sans rien regler sur les donnees.
une_analyse <- function(y, g, draws = 1500L, seed = 1L) {
	rs <- c(0.25, 0.5, 1)
	st <- .bay_laws_stats(y, g)
	tab <- .bay_laws_table(st, rs, 1.5, draws, seed)
	res <- .bay_laws_average(tab, 0.5)
	h1 <- tab$hypothesis == "H1"
	mix <- tab[tab$r == 0.5, , drop = FALSE]
	for (l in unique(tab$law)) {
		lm <- tab$log_ml[h1 & tab$law == l]
		w <- mix$law == l & mix$hypothesis == "H1"
		mix$log_ml[w] <- .bay_log_sum_exp(lm) - log(3)
	}
	res$p_mix <- .bay_laws_average(mix, 0.5)$p_h1
	res
}

if (partie %in% c("A", "A2")) {
	gen <- function(i, k, n, graine) {
		set.seed(i * 7919L + k * 13L + n + graine)
		law <- sample(.bay_law_names(), 1L)
		h1 <- runif(1) < 0.5
		het <- law %in% c("normal_het", "t_het")
		x <- rgamma(k, 1.5)
		vth <- x / sum(x)
		cc <- if (het) 1 / (k * vth) else rep(1, k)
		sig <- sqrt(cc)
		a <- if (h1) {
			g <- 1 / rgamma(1, 0.5, rate = 0.5^2 / 2)
			rnorm(k, 0, sqrt(g) * sig)
		} else rep(0, k)
		if (law %in% c("t_eq", "t_het")) {
			nu <- 2 + rexp(1)
			y <- unlist(lapply(seq_len(k), function(j) a[j] + sig[j] * sqrt((nu - 2) / nu) * rt(n, nu)))
		} else {
			y <- unlist(lapply(seq_len(k), function(j) rnorm(n, a[j], sig[j])))
		}
		g <- rep(seq_len(k), each = n)
		res <- tryCatch(une_analyse(y, g), error = function(e) NULL)
		if (is.null(res)) return(c(h1 = h1, p = NA, law = match(law, .bay_law_names()), rep(NA, 4)))
		c(h1 = h1, p = res$p_h1, law = match(law, .bay_law_names()), res$p_law[.bay_law_names()])
	}
	clusterExport(cl, c("gen", "une_analyse"))
	cas <- if (partie == "A") {
		list(c(2, 10), c(2, 30), c(3, 10), c(3, 30), c(5, 10), c(5, 30))
	} else {
		list(c(2, 30))
	}
	graine <- if (partie == "A") 0L else 500000L
	sortie <- character(0)
	for (cs in cas) {
		k <- cs[1L]; n <- cs[2L]
		r <- par_tranches(sprintf("%s_k%d_n%d", partie, k, n), N, gen, k = k, n = n, graine = graine)
		ok <- !is.na(r[, "p"])
		rr <- r[ok, , drop = FALSE]
		sortie <- c(sortie, sprintf("k=%d n=%d : %d jeux, %d echecs ; moyenne P(H1) %.3f (attendu %.3f +- %.3f) ; moyenne P(lois) %s ; frequence lois %s",
			k, n, nrow(rr), sum(!ok), mean(rr[, "p"]), mean(rr[, "h1"]), 2 * sd(rr[, "p"]) / sqrt(nrow(rr)),
			paste(round(colMeans(rr[, 4:7]), 3), collapse = " "),
			paste(round(tabulate(rr[, "law"], 4) / nrow(rr), 3), collapse = " ")))
		bins <- cut(rr[, "p"], c(0, .2, .4, .6, .8, .95, 1), include.lowest = TRUE)
		for (b in levels(bins)) {
			w <- bins == b
			if (!any(w)) next
			# z avec la variance annoncee, pas celle observee, qui s'effondre pres de 0 et 1.
			se <- sqrt(mean(rr[w, "p"] * (1 - rr[w, "p"])) / sum(w))
			sortie <- c(sortie, sprintf("   %-10s n=%4d annonce %.3f frequence %.3f z=%.2f", b, sum(w),
				mean(rr[w, "p"]), mean(rr[w, "h1"]), (mean(rr[w, "h1"]) - mean(rr[w, "p"])) / max(se, 1e-9)))
		}
		ecrit(sortie)
	}
}

if (partie == "B") {
	bruit <- function(n, type) {
		switch(type,
			normal = rnorm(n),
			t3 = rt(n, 3) / sqrt(3),
			lognormal = (rlnorm(n, 0, 0.8) - exp(0.32)) / sqrt((exp(0.64) - 1) * exp(0.64)))
	}
	gen <- function(i, sc) {
		set.seed(i * 104729L + sc$id)
		k <- length(sc$n)
		y <- unlist(lapply(seq_len(k), function(j) sc$m[j] + sc$sd[j] * bruit(sc$n[j], sc$type)))
		g <- factor(rep(seq_len(k), times = sc$n))
		res <- une_analyse(y, g)
		gauss <- .bay_bf_oneway_zs(y, g, 0.5)
		c(p = res$p_h1, pmix = res$p_mix, pg = plogis(gauss), aov = anova(lm(y ~ g))[1L, 5L],
			welch = oneway.test(y ~ g)$p.value, kw = kruskal.test(y ~ g)$p.value)
	}
	clusterExport(cl, c("gen", "une_analyse", "bruit"))
	sc <- list()
	add <- function(lab, n, m, sd, type) sc[[length(sc) + 1L]] <<- list(id = length(sc) + 1L, lab = lab, n = n, m = m, sd = sd, type = type)
	# Les cas decisifs d'abord, pour qu'une interruption ne les perde pas.
	add("hetero (5,5,40) sd(4,4,1) H0", c(5, 5, 40), c(0, 0, 0), c(4, 4, 1), "normal")
	add("hetero (5,5,40) sd(1,1,4) H0", c(5, 5, 40), c(0, 0, 0), c(1, 1, 4), "normal")
	add("hetero (10,30,50) sd(3,1,1) H0", c(10, 30, 50), c(0, 0, 0), c(3, 1, 1), "normal")
	add("hetero (2,2,40) sd(4,4,1) H0", c(2, 2, 40), c(0, 0, 0), c(4, 4, 1), "normal")
	add("(2x6) sd 1 H0", rep(2, 6), rep(0, 6), rep(1, 6), "normal")
	add("hetero (5,5,40) sd(4,4,1) effet 4 groupe 1", c(5, 5, 40), c(4, 0, 0), c(4, 4, 1), "normal")
	add("hetero (5,5,40) sd(4,4,1) effet 8 groupe 1", c(5, 5, 40), c(8, 0, 0), c(4, 4, 1), "normal")
	add("hetero (30,30,30) sd(1,1,3) d=.8 entre 1 et 2", c(30, 30, 30), c(.8, 0, 0), c(1, 1, 3), "normal")
	for (type in c("normal", "t3", "lognormal")) for (nn in c(10L, 30L, 100L)) for (k in c(2L, 3L)) for (d in c(0, 0.5, 0.8)) {
		add(sprintf("%s k=%d n=%d d=%.1f", type, k, nn, d), rep(nn, k), c(rep(0, k - 1L), d), rep(1, k), type)
	}
	add("k=5 n=10 H0", rep(10, 5), rep(0, 5), rep(1, 5), "normal")
	add("k=5 n=30 H0", rep(30, 5), rep(0, 5), rep(1, 5), "normal")
	add("k=5 n=30 d=.8", rep(30, 5), c(0, 0, 0, 0, .8), rep(1, 5), "normal")
	add("k=8 n=20 H0", rep(20, 8), rep(0, 8), rep(1, 8), "normal")
	sortie <- c("P>0,95 et p<0,05 ne sont pas des seuils equivalents : comparer aussi p<0,005.",
		"mix = hyper-a priori uniforme sur r = 1/4, 1/2, 1 ; gaussien = variance commune seule.",
		"scenario | P>.5 P>.95 P<.05 mediane | mix P>.5 P>.95 P<.05 | gaussien P>.5 P>.95 | aov p<.05 p<.005 | Welch p<.05 p<.005 | KW p<.05")
	for (s in sc) {
		r <- par_tranches(sprintf("B_%03d", s$id), N, gen, sc = s)
		sortie <- c(sortie, sprintf("%-45s | %.3f %.3f %.3f %.3f | %.3f %.3f %.3f | %.3f %.3f | %.3f %.3f | %.3f %.3f | %.3f", s$lab,
			mean(r[, "p"] > .5), mean(r[, "p"] > .95), mean(r[, "p"] < .05), median(r[, "p"]),
			mean(r[, "pmix"] > .5), mean(r[, "pmix"] > .95), mean(r[, "pmix"] < .05),
			mean(r[, "pg"] > .5), mean(r[, "pg"] > .95),
			mean(r[, "aov"] < .05), mean(r[, "aov"] < .005), mean(r[, "welch"] < .05), mean(r[, "welch"] < .005),
			mean(r[, "kw"] < .05)))
		ecrit(sortie)
	}
}

if (partie == "C") {
	gen <- function(i, n, k, d) {
		set.seed(i * 6101L + n + k * 7L + round(100 * d))
		m <- c(rep(0, k - 1L), d)
		y <- unlist(lapply(seq_len(k), function(j) rnorm(n, m[j], 1)))
		g <- factor(rep(seq_len(k), each = n))
		r <- m.bay(y, g, verbose = FALSE, draws = 1000L)
		s <- r$sensitivity
		c(p = r$prob_effect, lo = s$range_lead[1L], hi = s$range_lead[2L],
			change = s$verdict == "changes", crosses = s$verdict == "crosses", stable = s$verdict == "stable",
			r_side = length(unique(s$grid[, "1:1"] > .5)) > 1L,
			u_side = length(unique(c(r$prob_effect, s$by_u) > .5)) > 1L)
	}
	clusterExport(cl, "gen")
	sortie <- "Part des jeux : hypothese en tete qui change, seuil franchi, stable ; P change de cote selon r (a cotes 1:1), selon u"
	for (n in c(10L, 30L, 100L)) for (k in c(2L, 3L)) for (d in c(0, .5, .8)) {
		r <- par_tranches(sprintf("C_n%d_k%d_d%d", n, k, round(10 * d)), N, gen, n = n, k = k, d = d)
		sortie <- c(sortie, sprintf("n=%d k=%d d=%.1f : change %.3f, seuil %.3f, stable %.3f ; cote selon r %.3f ; selon u %.3f", n, k, d,
			mean(r[, "change"]), mean(r[, "crosses"]), mean(r[, "stable"]), mean(r[, "r_side"]), mean(r[, "u_side"])))
		ecrit(sortie)
	}
}

if (partie == "D") {
	sortie <- "B1 : petit groupe eloigne, plus il s'eloigne, plus H1 doit gagner"
	for (pl in list(c(10, 10, 10, 10, 3), c(30, 30, 2), c(20, 20, 20, 20, 3), c(rep(10, 7), 3), c(rep(5, 7), 2))) {
		for (d in c(5, 50, 500, 5000)) {
			set.seed(7)
			k <- length(pl)
			y <- c(rnorm(sum(pl[-k])), rnorm(pl[k], d, 1))
			g <- rep(seq_len(k), pl)
			r <- une_analyse(y, g, draws = 8000L)
			sortie <- c(sortie, sprintf("n=(%s) d=%5g : log BF %8.2f, P(H1) %.4f, P(variances inegales) %.3f, err %.1f%%, zs %.1f, Welch p %.2g",
				paste(pl, collapse = ","), d, r$log_bf, r$p_h1, sum(r$p_law[c("normal_het", "t_het")]), r$error_pct,
				.bay_bf_oneway_zs(y, g, .5), if (min(pl) > 1) oneway.test(y ~ factor(g))$p.value else NA))
			ecrit(sortie)
		}
	}
	sortie <- c(sortie, "", "B2 : un groupe disperse sans rapport ne doit pas changer la reponse sur les deux autres")
	set.seed(3)
	a <- rnorm(20); b <- rnorm(20); a <- (a - mean(a)) / sd(a); b <- (b - mean(b)) / sd(b) + 1.1
	r2 <- une_analyse(c(a, b), rep(1:2, each = 20), draws = 8000L)
	sortie <- c(sortie, sprintf("deux groupes seuls : log BF %.2f, P(H1) %.3f", r2$log_bf, r2$p_h1))
	for (S in c(1, 10, 100, 1e3, 1e4, 1e6)) {
		set.seed(5); c3 <- rnorm(20, 0, S)
		r <- une_analyse(c(a, b, c3), rep(1:3, each = 20), draws = 8000L)
		sortie <- c(sortie, sprintf("S=%g : log BF %.2f, P(H1) %.3f, err %.1f%%", S, r$log_bf, r$p_h1, r$error_pct))
	}
	sortie <- c(sortie, "", "Cas bruts")
	f <- function(y, g) { r <- tryCatch(une_analyse(y, g, draws = 8000L), error = function(e) conditionMessage(e)); if (is.character(r)) r else sprintf("log BF %.3f, P %.4f, err %.1f%%", r$log_bf, r$p_h1, r$error_pct) }
	set.seed(9)
	sortie <- c(sortie,
		paste("c(0,1e-160,2e-160,1,1+1e-10,1) :", f(c(0, 1e-160, 2e-160, 1, 1 + 1e-10, 1), rep(1:2, each = 3))),
		paste("donnees constantes a l'arrondi :", f(c(.3, .3, .3, .1 + .2, .1 + .2, .1 + .2), rep(1:2, each = 3))),
		paste("N = k+1, (1,1+1e-15,5) :", f(c(1, 1 + 1e-15, 5), c(1, 1, 2))),
		paste("N = k+1, (1,1,5) :", f(c(1, 1, 5), c(1, 1, 2))),
		paste("une valeur aberrante, n=20, k=2 :", f(c(rnorm(19), 40, rnorm(20)), rep(1:2, each = 20))),
		paste("groupe de n=1 :", f(c(rnorm(10), 3), rep(1:2, c(10, 1)))),
		paste("(2,2) :", f(c(0, 1, 5, 6), rep(1:2, each = 2))),
		paste("n=1e4 par groupe, k=2 :", f(rnorm(2e4), rep(1:2, each = 1e4))))
	ecrit(sortie)
}
stopCluster(cl)

# agent_50 : m.bay a plusieurs facteurs, voie normale par BayesFactor.
#   A : calibration sous l'a priori. Modele tire selon la loi HOP (une Beta(1,1) par ordre,
#       heredite forte), g_t ~ InvGamma(1/2, r^2/2) par terme, effets ~ N(0, g_t) sur le
#       codage de BayesFactor (fullDesignMatrix), sigma = 1. On compare P(terme | donnees)
#       a la frequence du terme vrai, par tranche de probabilite.
#   B : effets fixes, plans desequilibres, a cote d'aov (type II, car::Anova) :
#       faux termes a P>.5 et P>.95, puissance P(vrai > .5) ; aov p<.05.
# Calcul par les fonctions de KefiR (.bay_fact_space, .bay_fact_logbf, .bay_fact_terms).
# Reprise possible : resultats/agent_50_cas/.
# Usage : Rscript agent_50_multi_validation.R A|B N coeurs
args <- commandArgs(TRUE)
QUOI <- if (length(args) >= 1L) args[1L] else "A"
N <- if (length(args) >= 2L) as.integer(args[2L]) else 600L
nc <- if (length(args) >= 3L) as.integer(args[3L]) else 2L
LIB <- "C:/Users/masse/Desktop/KefiR/mbay-lib"
.libPaths(c(LIB, .libPaths()))
library(parallel)
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
CAS <- file.path(BASE, "dev/mbay/scripts/resultats/agent_50_cas")
dir.create(CAS, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(nc)
invisible(clusterCall(cl, function(LIB, BASE) {
	.libPaths(c(LIB, .libPaths()))
	source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R"))
	.kefir_bay(BASE)
	suppressMessages(library(BayesFactor))
	NULL
}, LIB, BASE))

plan_data <- function(plan) {
	cells <- switch(plan,
		"2x2_n5" = list(lv = c(2, 2), n = rep(5, 4)),
		"2x2_deseq" = list(lv = c(2, 2), n = c(4, 8, 12, 16)),
		"2x3_deseq" = list(lv = c(2, 3), n = c(3, 9, 5, 7, 10, 4)),
		"2x2x2_deseq" = list(lv = c(2, 2, 2), n = c(3, 6, 4, 8, 5, 7, 4, 6)))
	nm <- LETTERS[seq_along(cells$lv)]
	g <- expand.grid(lapply(setNames(cells$lv, nm), function(l) factor(seq_len(l))))
	d <- g[rep(seq_len(nrow(g)), cells$n), , drop = FALSE]
	rownames(d) <- NULL
	d
}
termes_de <- function(k) {
	nm <- LETTERS[seq_len(k)]
	unlist(lapply(seq_len(k), function(o) apply(utils::combn(nm, o), 2L, paste, collapse = ":")))
}

# Model drawn from the HOP prior, order by order.
tire_modele <- function(sp) {
	m <- integer(0)
	for (o in sort(unique(sp$order))) {
		can <- which(sp$order == o & vapply(seq_along(sp$keys), function(i) all(sp$parents[[i]] %in% m), logical(1)))
		if (!length(can)) next
		p <- stats::runif(1)
		m <- c(m, can[stats::runif(length(can)) < p])
	}
	m
}

un_A <- function(i, plan, r = 0.5) {
	set.seed(500000L + 1000L * nchar(plan) + i)
	d <- plan_data(plan)
	k <- ncol(d)
	termes <- termes_de(k)
	sp <- .bay_fact_space(termes)
	vrai <- tire_modele(sp)
	X <- as.matrix(BayesFactor:::fullDesignMatrix(stats::as.formula(paste("~", paste(termes, collapse = "+"))),
		d, setNames(rep("fixed", k), names(d))))
	col_terme <- vapply(strsplit(colnames(X), ".&.", fixed = TRUE), function(p) {
		paste(sort(sub("_redu_[0-9]+$", "", p)), collapse = ":")
	}, "")
	mu <- numeric(nrow(d))
	for (t in vrai) {
		cols <- which(col_terme == sp$keys[t])
		g <- 1 / stats::rgamma(1, 1 / 2, rate = r^2 / 2)
		mu <- mu + X[, cols, drop = FALSE] %*% stats::rnorm(length(cols), 0, sqrt(g))
	}
	y <- as.numeric(mu) + stats::rnorm(nrow(d))
	bf <- .bay_fact_logbf(y, d, sp, r, 10000L, i)
	e <- .bay_fact_terms(sp, bf$log_bf[, 1L], .bay_fact_prior(sp, "hop"))
	p <- e$p_h1[match(sp$keys, .bay_term_key(e$term))]
	p0 <- .bay_fact_null(bf$log_bf[, 1L], .bay_fact_prior(sp, "hop"))
	c(p, as.numeric(seq_along(termes) %in% vrai), p_nul = p0, nul_vrai = as.numeric(!length(vrai)),
		err = max(bf$error[, 1L]))
}

# Effets fixes : vrais termes de taille f, codage +-f pour 2 niveaux, (-f, 0, f) pour 3.
un_B <- function(i, plan, vrai, f) {
	set.seed(700000L + 1000L * nchar(plan) + 37L * length(vrai) + round(100 * f) + i)
	d <- plan_data(plan)
	k <- ncol(d)
	termes <- termes_de(k)
	sp <- .bay_fact_space(termes)
	code <- function(v) { l <- nlevels(v); c(-1, 0, 1)[if (l == 2) c(1, 3) else 1:3][as.integer(v)] }
	mu <- numeric(nrow(d))
	for (t in vrai) mu <- mu + f * Reduce(`*`, lapply(strsplit(t, ":", fixed = TRUE)[[1]], function(v) code(d[[v]])))
	y <- mu + stats::rnorm(nrow(d))
	bf <- .bay_fact_logbf(y, d, sp, 0.5, 10000L, i)
	e <- .bay_fact_terms(sp, bf$log_bf[, 1L], .bay_fact_prior(sp, "hop"))
	p <- e$p_h1[match(sp$keys, .bay_term_key(e$term))]
	fit <- stats::lm(stats::as.formula(paste("y ~", paste(termes, collapse = "+"))), data = cbind(y = y, d),
		contrasts = setNames(rep(list("contr.sum"), k), names(d)))
	pv <- car::Anova(fit, type = 2)[termes, "Pr(>F)"]
	c(p, pv)
}
clusterExport(cl, c("plan_data", "termes_de", "tire_modele", "un_A", "un_B"))

lot <- function(nom, fun, NN, ...) {
	r <- NULL
	for (d0 in seq(1L, NN, by = 50L)) {
		fic <- file.path(CAS, sprintf("%s_%05d.rds", nom, d0))
		if (!file.exists(fic)) saveRDS(do.call(rbind, parLapplyLB(cl, d0:min(d0 + 49L, NN), fun, ...)), fic)
		r <- rbind(r, readRDS(fic))
	}
	r
}

if (QUOI == "A") {
	lignes <- c(R.version.string, sprintf("A. Calibration sous l'a priori HOP, r = 0.5, N = %d par plan", N),
		"plan : par tranche de P(terme) : moyenne annoncee / frequence du terme vrai (nombre)")
	for (plan in c("2x2_n5", "2x3_deseq", "2x2_deseq")) {
		r <- lot(paste0("A_", plan), un_A, N, plan = plan)
		nt <- (ncol(r) - 3L) / 2L
		p <- as.vector(r[, seq_len(nt)]); v <- as.vector(r[, nt + seq_len(nt)])
		b <- cut(p, c(0, .1, .3, .5, .7, .9, 1), include.lowest = TRUE)
		tb <- vapply(split(seq_along(p), b), function(ix) c(mean(p[ix]), mean(v[ix]), length(ix)), numeric(3))
		z <- (sum(v) - sum(p)) / sqrt(sum(p * (1 - p)))
		lignes <- c(lignes, sprintf("%s : %s | ensemble : annonce %.3f, vrai %.3f, z = %.2f | nul : annonce %.3f, vrai %.3f | erreur BF max %.3f",
			plan, paste(sprintf("[%s] %.3f/%.3f (%d)", colnames(tb), tb[1, ], tb[2, ], tb[3, ]), collapse = " "),
			mean(p), mean(v), z, mean(r[, "p_nul"]), mean(r[, "nul_vrai"]), max(r[, "err"])))
		writeLines(lignes, file.path(BASE, "dev/mbay/scripts/resultats/agent_50_A.txt"))
	}
} else {
	SCEN <- list(
		list("2x2_deseq", character(0), 0), list("2x2_deseq", c("A", "B"), .5), list("2x2_deseq", "A:B", .5),
		list("2x2_deseq", c("A", "B"), .25),
		list("2x3_deseq", character(0), 0), list("2x3_deseq", c("A", "B"), .5), list("2x3_deseq", c("A", "B", "A:B"), .5),
		list("2x2x2_deseq", character(0), 0), list("2x2x2_deseq", c("A", "B", "C"), .5))
	lignes <- c(R.version.string, sprintf("B. Effets fixes, plans desequilibres, N = %d par cas (2x2x2 : N/2)", N),
		"cas : m.bay P(un faux terme)>.5, >.95, P(vrai)>.5 moyenne | aov type II : un faux p<.05, vrai p<.05 moyenne")
	for (s in SCEN) {
		NN <- if (s[[1]] == "2x2x2_deseq") N %/% 2L else N
		nom <- sprintf("B_%s_%s_f%d", s[[1]], if (length(s[[2]])) gsub(":", "", paste(s[[2]], collapse = "-")) else "H0", round(100 * s[[3]]))
		r <- lot(nom, un_B, NN, plan = s[[1]], vrai = s[[2]], f = s[[3]])
		termes <- termes_de(if (s[[1]] == "2x2x2_deseq") 3L else 2L)
		nt <- length(termes)
		faux <- which(!termes %in% s[[2]]); vrais <- which(termes %in% s[[2]])
		pf <- if (length(faux)) apply(r[, faux, drop = FALSE], 1L, max) else rep(0, nrow(r))
		lignes <- c(lignes, sprintf("%s (N=%d) : %.3f %.3f %s | %.3f %s", nom, nrow(r), mean(pf > .5), mean(pf > .95),
			if (length(vrais)) sprintf("%.3f", mean(r[, vrais] > .5)) else "-",
			mean(apply(r[, nt + faux, drop = FALSE] < .05, 1L, any)),
			if (length(vrais)) sprintf("%.3f", mean(r[, nt + vrais] < .05)) else "-"))
		writeLines(lignes, file.path(BASE, "dev/mbay/scripts/resultats/agent_50_B.txt"))
	}
}
stopCluster(cl)

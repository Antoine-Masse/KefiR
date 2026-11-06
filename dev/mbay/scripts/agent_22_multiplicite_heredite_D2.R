# agent_22 : variante d agent_20, heredite CONDITIONNELLE a 1/2 (loi annoncee dans D2 : A:B a 1/8) ; multiplicite en plan 2x2x2 sous H0 globale, recoupee sur anovaBF.
# Usage : Rscript agent_20_multiplicite_2x2x2.R <N_total> <n_jeux> <graine> [budget_s]
# Sortie : resultats/agent_20_N<N_total>.rds (resume par jeu) et .txt
#
# LOI A HEREDITE FORTE (definition exacte, prior_prob de .bay_inclusion, type "all") :
#   pi = 1/2.
#   Chaque facteur principal entre independamment avec la probabilite pi.
#   Chaque interaction (k = 2 ou 3 facteurs) entre avec la probabilite
#   pi = 1/2 SI tous ses termes parents (tous les sous-ensembles stricts de ses facteurs)
#   sont dans le modele, et avec la probabilite 0 sinon. Tirages independants,
#   conditionnellement au niveau inferieur.
#   Donc A:B : 1/2 si A et B presents (P(A:B)=1/8) ; A:B:C : 1/2 si les trois doubles presentes (P=1/128).
#   Poids d'un modele = produit des probabilites d'entree et de non-entree de tous les
#   termes eligibles. La somme sur les 19 modeles vaut 1 (verifie ci-dessous).
#
# Espace : 19 modeles respectant la marginalite (nul compris). anovaBF("withmain") donne
# 25 modeles non nuls ; on retire les 7 ou A:B:C figure sans ses trois interactions doubles.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress = FALSE)
args <- commandArgs(TRUE)
Ntot <- as.integer(args[1]); nrep <- as.integer(args[2]); seed <- as.integer(args[3])
budget <- if (length(args) >= 4) as.numeric(args[4]) else Inf
racine <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
for (f in c("sys_msg.R", "sys_exit.R", "sys_bay_core.R")) source(file.path(racine, "R", f), encoding = "UTF-8")
sortie <- file.path(racine, "dev/mbay/scripts/resultats")

rep_cell <- Ntot / 8; stopifnot(rep_cell == round(rep_cell))
des <- expand.grid(A = gl(2, 1), B = gl(2, 1), C = gl(2, 1), r = seq_len(rep_cell))
des$r <- NULL
cn <- c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C")
ord <- c(1, 1, 1, 2, 2, 2, 3)

# Espace des modeles et poids a heredite forte
marg <- as.matrix(expand.grid(rep(list(0:1), 7))); colnames(marg) <- cn
ok <- apply(marg, 1, function(z) z["A:B"] <= z["A"] * z["B"] && z["A:C"] <= z["A"] * z["C"] &&
	z["B:C"] <= z["B"] * z["C"] && z["A:B:C"] <= z["A:B"] * z["A:C"] * z["B:C"])
marg <- marg[ok, , drop = FALSE]; stopifnot(nrow(marg) == 19)
poids_hered <- function(z, pi = 1/2) {
	p <- 1
	for (t in cn) {
		k <- ord[match(t, cn)]
		elig <- if (k == 1) TRUE else {
			f <- strsplit(t, ":", fixed = TRUE)[[1]]
			par <- unlist(lapply(seq_len(k - 1), function(j) combn(f, j, paste, collapse = ":")))
			all(z[par] == 1)
		}
		if (elig) p <- p * (if (z[t] == 1) pi else 1 - pi) else stopifnot(z[t] == 0)
	}
	p
}
w_h <- apply(marg, 1, poids_hered); stopifnot(abs(sum(w_h) - 1) < 1e-12)
cle <- function(z) paste(cn[z == 1], collapse = " + ")
cles <- apply(marg, 1, cle)
termes <- as.list(cles); termes <- lapply(termes, function(s) if (s == "") character(0) else strsplit(s, " + ", fixed = TRUE)[[1]])
i_nul <- which(cles == "")

extrait <- function(y) {
	d <- des; d$y <- y
	b <- anovaBF(y ~ A * B * C, data = d, whichModels = "withmain", rscaleFixed = "medium", progress = FALSE)
	tab <- b@bayesFactor
	norm_nom <- function(s) {
		v <- trimws(strsplit(s, "+", fixed = TRUE)[[1]])
		paste(v[order(ord[match(v, cn)], match(v, cn))], collapse = " + ")
	}
	nn <- vapply(rownames(tab), norm_nom, character(1))
	lb <- tab$bf[match(cles, nn)]; err <- tab$error[match(cles, nn)]
	lb[i_nul] <- 0; err[i_nul] <- 0
	stopifnot(sum(is.na(match(cles[-i_nul], nn))) == 0)
	list(lb = lb, err = err)
}

P_mat <- function(lb) {
	tt <- function(type, pp = NULL) {
		r <- .bay_inclusion(lb, termes, type = type, prior_prob = pp)
		r$p_h1[match(cn, r$term)]
	}
	rbind(apparie = tt("matched"), uniforme = tt("all"), heredite = tt("all", w_h))
}

set.seed(seed)
P <- array(NA_real_, c(nrep, 3, 7), list(NULL, c("apparie", "uniforme", "heredite"), cn))
pv <- matrix(NA_real_, nrep, 7, dimnames = list(NULL, cn))
post_nul <- matrix(NA_real_, nrep, 2, dimnames = list(NULL, c("uniforme", "heredite")))
LB <- matrix(NA_real_, nrep, 19, dimnames = list(NULL, cles))
err_med <- err_max <- rep(NA_real_, nrep); statut <- rep("non fait", nrep); msg_err <- character(0)
t0 <- Sys.time(); fait <- 0L; arret <- FALSE
for (i in seq_len(nrep)) {
	if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > budget) { arret <- TRUE; break }
	y <- rnorm(Ntot)
	d <- des; d$y <- y
	pv[i, ] <- summary(aov(y ~ A * B * C, data = d))[[1]][1:7, 5]
	res <- tryCatch(extrait(y), error = function(e) e)
	fait <- i
	if (inherits(res, "error")) { statut[i] <- "erreur"; msg_err <- c(msg_err, conditionMessage(res)); next }
	statut[i] <- if (anyNA(res$lb)) "NA" else "ok"
	err_med[i] <- median(res$err[-i_nul], na.rm = TRUE); err_max[i] <- max(res$err, na.rm = TRUE)
	if (anyNA(res$lb)) next
	LB[i, ] <- res$lb
	P[i, , ] <- P_mat(res$lb)
	lp_u <- res$lb - .bay_log_sum_exp(res$lb)
	lp_h <- res$lb + log(w_h); lp_h <- lp_h - .bay_log_sum_exp(lp_h)
	post_nul[i, ] <- exp(c(lp_u[i_nul], lp_h[i_nul]))
	if (i %% 200 == 0) {
		cat(format(Sys.time(), "%H:%M:%S"), "jeu", i, "\n")
		saveRDS(list(P = P, pv = pv, post_nul = post_nul, err_med = err_med, err_max = err_max, statut = statut, fait = i),
			file.path(sortie, sprintf("agent_22_N%d.partiel.rds", Ntot)))
	}
}
duree <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

saveRDS(list(Ntot = Ntot, seed = seed, nrep_demande = nrep, fait = fait, arret_budget = arret, duree_s = duree,
	LB = LB[seq_len(fait), , drop = FALSE], 	P = P[seq_len(fait), , , drop = FALSE], pv = pv[seq_len(fait), , drop = FALSE],
	post_nul = post_nul[seq_len(fait), , drop = FALSE], err_med = err_med[seq_len(fait)], err_max = err_max[seq_len(fait)],
	statut = statut[seq_len(fait)], msg_err = msg_err, poids_hered = setNames(w_h, cles), R = R.version.string,
	BayesFactor = as.character(packageVersion("BayesFactor"))),
	file.path(sortie, sprintf("agent_22_N%d.rds", Ntot)))
unlink(file.path(sortie, sprintf("agent_22_N%d.partiel.rds", Ntot)))

zz <- file(file.path(sortie, sprintf("agent_22_N%d.txt", Ntot)), "wt", encoding = "UTF-8"); sink(zz)
cat(R.version.string, "; BayesFactor", as.character(packageVersion("BayesFactor")), "; graine", seed, "\n")
cat("N total", Ntot, "; jeux faits", fait, "sur", nrep, "; arret budget :", arret, "; duree", round(duree), "s ;",
	round(duree / max(fait, 1), 3), "s/jeu\n")
cat("Statut des jeux :\n"); print(table(statut[seq_len(fait)]))
if (length(msg_err)) { cat("Messages d'erreur :\n"); print(table(msg_err)) }
cat("Erreur BayesFactor (colonne error, proportionnelle) : mediane par jeu des medianes", signif(median(err_med, na.rm = TRUE), 3),
	"; maximum sur tous jeux et modeles", signif(max(err_max, na.rm = TRUE), 3),
	"; q99 du max par jeu", signif(quantile(err_max, .99, na.rm = TRUE), 3), "\n")
ic <- function(x) { x <- x[!is.na(x)]; k <- sum(x); bt <- binom.test(k, length(x))$conf.int; sprintf("%.3f [%.3f ; %.3f] (%d/%d)", k / length(x), bt[1], bt[2], k, length(x)) }
idx <- seq_len(fait)
cat("\nTaux par jeu d'au moins un terme au-dessus du seuil\n")
for (loi in dimnames(P)[[2]]) for (s in c(.5, .75, .95)) {
	a <- apply(P[idx, loi, , drop = FALSE], 1, function(v) if (anyNA(v)) NA else any(v > s))
	cat(sprintf(" %-9s P>%.2f : %s\n", loi, s, ic(a)))
}
a <- apply(pv[idx, , drop = FALSE], 1, function(v) any(v < .05))
cat(" aov p<0,05 au moins une fois :", ic(a), "; par terme", round(mean(pv[idx, ] < .05, na.rm = TRUE), 4), "\n")
cat("\nP(terme|D) par terme : moyenne et taux P>0,5 ; a priori de la loi\n")
pri <- function(pp) { r <- .bay_inclusion(rep(0, 19), termes, type = "all", prior_prob = pp); r$prior_h1[match(cn, r$term)] }
for (loi in dimnames(P)[[2]]) {
	pr <- switch(loi, apparie = rep(.5, 7), uniforme = pri(NULL), heredite = pri(w_h))
	print(round(rbind(apriori = pr, moyenne = apply(P[idx, loi, , drop = FALSE], 3, mean, na.rm = TRUE),
		mediane = apply(P[idx, loi, , drop = FALSE], 3, median, na.rm = TRUE),
		taux_P_sup_05 = apply(P[idx, loi, , drop = FALSE] > .5, 3, mean, na.rm = TRUE)), 4)); cat(" ^", loi, "\n")
}
cat("\nNombre moyen de termes P>0,5 (apparie) :", round(mean(rowSums(P[idx, "apparie", ] > .5), na.rm = TRUE), 3), "\n")
cat("\nP(modele nul | D) : a priori uniforme", round(1/19, 4), "; heredite", round(w_h[i_nul], 4), "\n")
for (l in colnames(post_nul)) cat(sprintf(" %-9s mediane %.4f ; decile bas %.4f ; moyenne %.4f\n", l,
	median(post_nul[idx, l], na.rm = TRUE), quantile(post_nul[idx, l], .1, na.rm = TRUE), mean(post_nul[idx, l], na.rm = TRUE)))
sink(); close(zz)
cat("FIN", duree, "s\n")

# agent_31 : cout de la loi retenue (D2, A et B a 1/2, A:B a 1/8) quand l'interaction existe vraiment.
# Plan 2x2, interaction pure : moyenne = f * sA * sB, sA et sB = +-1, sd = 1 (f : taille de Cohen).
# Compare P(A:B | D) et P(A joue un role | D) sous heredite et sous uniforme, a cote de aov p<0,05.
# Lois de la table : nul, A, B, A+B, A+B+A:B ; heredite 1/4, 1/4, 1/4, 1/8, 1/8 ; uniforme 1/5 chacun.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress = FALSE)
set.seed(20261131)
cles <- c("", "A", "B", "A + B", "A + B + A:B")
w_h <- c(1/4, 1/4, 1/4, 1/8, 1/8); w_u <- rep(1/5, 5)
stopifnot(abs(sum(w_h) - 1) < 1e-12)
post <- function(lb, w) { l <- lb + log(w); exp(l - max(l)) / sum(exp(l - max(l))) }
un_jeu <- function(n, f) {
	d <- expand.grid(A = gl(2, 1), B = gl(2, 1), r = seq_len(n))
	sA <- ifelse(d$A == 1, -1, 1); sB <- ifelse(d$B == 1, -1, 1)
	d$y <- rnorm(nrow(d), f * sA * sB, 1)
	tab <- anovaBF(y ~ A * B, data = d, whichModels = "withmain", progress = FALSE)@bayesFactor
	lb <- c(0, tab$bf[match(cles[-1], rownames(tab))])
	stopifnot(!anyNA(lb))
	ph <- post(lb, w_h); pu <- post(lb, w_u)
	c(AB_h = ph[5], AB_u = pu[5], A_h = sum(ph[c(2, 4, 5)]), A_u = sum(pu[c(2, 4, 5)]),
	  p_ab = summary(aov(y ~ A * B, data = d))[[1]][3, 5])
}
N <- 200
cat(R.version.string, "; BayesFactor", as.character(packageVersion("BayesFactor")), "; N =", N, "par case ; graine 20261131\n")
cat("A priori : heredite P(A:B)=0,125, P(A)=0,5 ; uniforme P(A:B)=0,2, P(A)=0,6\n\n")
for (n in c(10, 25, 50)) for (f in c(0, 0.25, 0.5)) {
	r <- t(replicate(N, un_jeu(n, f)))
	cat(sprintf("n=%d/case f=%.2f | P(A:B)>0,5 : heredite %.3f, uniforme %.3f | >0,95 : %.3f, %.3f | mediane : %.3f, %.3f | P(A role)>0,5 : %.3f, %.3f | aov A:B p<0,05 : %.3f\n",
		n, f, mean(r[, "AB_h"] > .5), mean(r[, "AB_u"] > .5), mean(r[, "AB_h"] > .95), mean(r[, "AB_u"] > .95),
		median(r[, "AB_h"]), median(r[, "AB_u"]), mean(r[, "A_h"] > .5), mean(r[, "A_u"] > .5), mean(r[, "p_ab"] < .05)))
}

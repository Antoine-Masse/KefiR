# G12 : comparaison de quatre moteurs pour l'independance (audit de 01effef par la session parente).
#   Rscript dev/gbay/scripts/sim_engines.R
# Moteurs, tous a marges libres ou neutres sous H0 et H1 :
#   bic     : -BIC/2 avec n = N, log BF = (G - ddl log N) / 2 (G11, moteur de 01effef) ;
#   robust  : BAS::bas.glm, Poisson, interaction contre effets principaux, a priori robuste
#             (Bayarri et al. 2012), n = N ;
#   hgn     : idem, a priori hyper-g/n (Li et Clyde 2018), alpha = 3, n = N ;
#   neutral : Dirichlet-multinomiale, a = 1, ou H0 recoit sur ses marges la loi que H1 leur donne
#             (Dirichlet(c a) sur les lignes, Dirichlet(r a) sur les colonnes) : forme close de
#             l'expert (chk_gbay3.R).
# Criteres (session parente) : sous H0, P(lien) > 0,95 dans moins de 1 % des tables ; sous un
# lien modere (V = 0,2 et 0,3), P(independance) > 0,95 dans au plus 5 % des tables, et une
# puissance raisonnable face au test G.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BAS))
set.seed(20261112)
R <- as.integer(Sys.getenv("GBAY_REPS", "1000"))
t0 <- Sys.time()

.G <- function(t) { E <- outer(rowSums(t), colSums(t)) / sum(t); 2 * sum(ifelse(t > 0, t * log(t / E), 0)) }
ldm <- function(n, a) lgamma(sum(a)) - lgamma(sum(a) + sum(n)) + sum(lgamma(a + n) - lgamma(a))
lb_neutral <- function(t, a = 1) {
	r <- nrow(t); cc <- ncol(t)
	ldm(as.vector(t), rep(a, r * cc)) - ldm(rowSums(t), rep(cc * a, r)) - ldm(colSums(t), rep(r * a, cc))
}
lb_bas <- function(t, prior) {
	d <- data.frame(n = as.numeric(t), A = factor(rep(seq_len(nrow(t)), ncol(t))), B = factor(rep(seq_len(ncol(t)), each = nrow(t))))
	N <- sum(t)
	bp <- if (prior == "robust") BAS::robust(n = N) else BAS::hyper.g.n(alpha = 3, n = N)
	f <- tryCatch(suppressWarnings(BAS::bas.glm(n ~ A * B, data = d, family = stats::poisson(), betaprior = bp, modelprior = BAS::uniform(),
		include.always = ~ A + B, force.heredity = TRUE, method = "BAS")), error = function(e) NULL)
	if (is.null(f)) return(NA_real_)
	f$logmarg[which.max(f$size)] - f$logmarg[which.min(f$size)]
}
use_bas <- Sys.getenv("GBAY_BAS", "1") != "0"
part <- Sys.getenv("GBAY_PART", "all")
engines <- function(t) {
	t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
	if (nrow(t) < 2L || ncol(t) < 2L) return(rep(NA_real_, 5))
	df <- (nrow(t) - 1) * (ncol(t) - 1)
	G <- .G(t)
	c(bic = (G - df * log(sum(t))) / 2, robust = if (use_bas) lb_bas(t, "robust") else NA_real_, hgn = if (use_bas) lb_bas(t, "hgn") else NA_real_, neutral = lb_neutral(t),
		gp = stats::pchisq(G, df, lower.tail = FALSE))
}
pr <- function(l) 1 / (1 + exp(-l))
nm <- c("bic", "robust", "hgn", "neutral")

# Cramer's V of a probability table.
cv <- function(P) { E <- outer(rowSums(P), colSums(P)); sqrt(sum((P - E)^2 / E) / (min(dim(P)) - 1)) }
# Diagonal link of strength V around margins r, c.
mk <- function(r, cc, V) {
	P0 <- outer(r, cc); D <- matrix(-1 / (length(cc) - 1), length(r), length(cc))
	for (i in seq_len(min(length(r), length(cc)))) D[i, i] <- 1
	f <- function(d) { P <- P0 * (1 + d * D); if (any(P <= 0)) return(1e6); cv(P / sum(P)) - V }
	d <- stats::uniroot(f, c(0, 0.999 * (length(cc) - 1)))$root
	P <- P0 * (1 + d * D); P / sum(P)
}

run <- function(P, N, rows_fixed = NULL) {
	r <- nrow(P)
	m <- vapply(seq_len(R), function(i) {
		t <- if (is.null(rows_fixed)) matrix(stats::rmultinom(1, N, as.numeric(P)), r) else
			t(vapply(seq_len(r), function(k) stats::rmultinom(1, rows_fixed, P[k, ] / sum(P[k, ]))[, 1], numeric(ncol(P))))
		engines(t)
	}, numeric(5))
	m[, !is.na(m[1L, ]), drop = FALSE]
}
# BAS fails on some tables (internal error, or NA): its rates are over the tables where it ran.
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summ <- function(m, lab) {
	out <- data.frame(lab, n_ok = ncol(m), gtest05 = mean(m[5L, ] < 0.05))
	for (e in seq_along(nm)) {
		p <- pr(m[e, ])
		out[[paste0(nm[e], "_p75")]] <- mean_na(p > 0.75)
		out[[paste0(nm[e], "_p95")]] <- mean_na(p > 0.95)
		out[[paste0(nm[e], "_h0strong")]] <- mean_na(p < 0.05)
		out[[paste0(nm[e], "_meanP")]] <- mean_na(p)
	}
	out
}

h0 <- NULL
if (part != "h1") for (dm in list(c(2, 2), c(3, 3), c(3, 4), c(4, 5))) for (mg in c("egales", "inegales")) for (N in c(30, 60, 200)) {
	rr <- if (mg == "egales") rep(1, dm[1]) else seq(1, 3, length.out = dm[1])
	cc <- if (mg == "egales") rep(1, dm[2]) else seq(1, 4, length.out = dm[2])
	P <- outer(rr / sum(rr), cc / sum(cc))
	h0 <- rbind(h0, summ(run(P, N), data.frame(table = paste0(dm[1], "x", dm[2]), margins = mg, N = N)))
	cat(".")
}
# Lignes fixees et equilibrees, le cas type du TP (5 groupes de 10).
if (part != "h1") for (k in c(2, 3)) {
	P <- outer(rep(0.2, 5), rep(1 / k, k))
	h0 <- rbind(h0, summ(run(P, 50, rows_fixed = 10), data.frame(table = paste0("5x", k), margins = "lignes fixees (10)", N = 50)))
}
cat("\n== Sous H0 ==\n"); print(h0, digits = 3, width = 250)

h1 <- NULL
if (part != "h0") for (dm in list(c(2, 2), c(3, 3), c(3, 4))) for (V in c(0.2, 0.3)) for (N in c(50, 100, 200)) {
	P <- mk(rep(1 / dm[1], dm[1]), rep(1 / dm[2], dm[2]), V)
	h1 <- rbind(h1, summ(run(P, N), data.frame(table = paste0(dm[1], "x", dm[2]), V = V, N = N)))
	cat(".")
}
cat("\n== Sous un lien (diagonal, V de Cramer donne) ==\n"); print(h1, digits = 3, width = 250)
saveRDS(list(h0 = h0, h1 = h1, R = R, seed = 20261112), paste0("dev/gbay/scripts/sim_engines_", part, ".rds"))
cat("\nDuree :", format(Sys.time() - t0), "\n")

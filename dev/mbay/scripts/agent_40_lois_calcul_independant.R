# agent_40 : les 8 modeles de R/sys_bay_laws.R contre un calcul independant.
# Reference ecrite separement, sans aucune fonction du paquet : quadrature directe
# de la vraisemblance brute (dnorm, dt) sur tous les parametres, 2 groupes, avec
# integrale de Cauchy sur l'effet standardise.
#   parametres k = 2 : mu, a = theta1 - theta2 (ou mu seul), log taubar, rho, nu
# Sous H1 : a_i | g ~ N(0, g sigma_i^2), donc avec mu plat la difference
# a = theta1 - theta2 suit, sachant g, N(0, g (sigma1^2 + sigma2^2)), et
# marginalement une Cauchy d'echelle r * sqrt(sigma1^2 + sigma2^2).
# Sous cette forme, mu plat se reparametre (mu, a) -> (theta1, theta2), jacobien 1,
# si mu est le point milieu pondere ; on l'ecrit directement sur (theta1, theta2).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(pkgload::load_all(".", quiet = TRUE, helpers = FALSE, attach_testthat = FALSE))

gl <- function(n) {
	i <- seq_len(n - 1L); b <- i / sqrt(4 * i^2 - 1)
	J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
	e <- eigen(J, symmetric = TRUE); o <- order(e$values)
	list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}
# grille produit sur [lo, hi]
nodes <- function(lo, hi, n) { g <- gl(n); list(x = (hi + lo) / 2 + (hi - lo) / 2 * g$x, w = (hi - lo) / 2 * g$w) }

lmean <- function(l) { t <- max(l); t + log(sum(exp(l - t))) }

ref_model <- function(y1, y2, h1, het, law, r = 0.5, u = 1.5, n_nodes = if (law == "t") 24L else 36L) {
	# Bornes larges autour des estimations ; donnees deja standardisees.
	s <- sqrt((sum((y1 - mean(y1))^2) + sum((y2 - mean(y2))^2)) / (length(y1) + length(y2) - 2))
	lik_group <- function(y, loc, ls2, nu) {
		# loc, ls2, nu : vecteurs de meme longueur (points de grille)
		if (law == "normal") {
			rowSums(matrix(dnorm(rep(y, each = length(loc)), rep(loc, length(y)),
				exp(0.5 * rep(ls2, length(y))), log = TRUE), length(loc), length(y)))
		} else {
			sc <- exp(0.5 * ls2) * sqrt((nu - 2) / nu)
			rowSums(matrix(dt((rep(y, each = length(loc)) - rep(loc, length(y))) / rep(sc, length(y)),
				rep(nu, length(y)), log = TRUE) - log(rep(sc, length(y))), length(loc), length(y)))
		}
	}
	lt_nodes <- nodes(-2.2 - 2 * log(s), 2.2 - 2 * log(s), n_nodes)
	if (het) {
		rho_nodes <- nodes(-5, 5, n_nodes)  # logit de rho
	} else {
		rho_nodes <- list(x = 0, w = 1)
	}
	if (law == "t") {
		nu_nodes <- nodes(-9, 2.6, 30L)   # log(nu - 2)
	} else {
		nu_nodes <- list(x = 0, w = 1)
	}
	span <- 4.5 * s * max(1 / sqrt(length(y1)), 1 / sqrt(length(y2)), 0.35)
	if (law == "t") span <- span * 1.4
	loc1 <- nodes(mean(y1) - span, mean(y1) + span, n_nodes)
	loc2 <- nodes(mean(y2) - span, mean(y2) + span, n_nodes)
	pool <- mean(c(y1, y2))
	vals <- numeric(0)
	for (irho in seq_along(rho_nodes$x)) for (inu in seq_along(nu_nodes$x)) {
		if (het) {
			rho <- plogis(rho_nodes$x[irho])
			ld <- dbeta(rho, u, u, log = TRUE) + log(rho) + log(1 - rho)
			c1 <- 1 / (2 * rho); c2 <- 1 / (2 * (1 - rho))
		} else { ld <- 0; c1 <- 1; c2 <- 1 }
		if (law == "t") {
			xv <- nu_nodes$x[inu]; nu <- 2 + exp(xv)
			lv <- -exp(xv) + xv
		} else { nu <- NA; lv <- 0 }
		base <- ld + lv + log(rho_nodes$w[irho]) + log(nu_nodes$w[inu])
		for (it in seq_along(lt_nodes$x)) {
			# Pour t : fenetre centree sur l'echelle s^2 = sigma^2 (nu - 2) / nu, sinon la grille
			# rate la region nu proche de 2, ou taubar est tres petit (premier essai : ecart de 0,08).
			lt <- lt_nodes$x[it] + if (law == "t") log((nu - 2) / nu) else 0
			ls1 <- log(c1) - lt; ls2 <- log(c2) - lt
			if (h1) {
				g1 <- rep(loc1$x, times = n_nodes); g2 <- rep(loc2$x, each = n_nodes)
				w <- log(rep(loc1$w, times = n_nodes)) + log(rep(loc2$w, each = n_nodes))
				# Cauchy sur (theta1 - theta2), echelle r sqrt(sigma1^2 + sigma2^2)
				scl <- r * sqrt(exp(ls1) + exp(ls2))
				prior <- dcauchy(g1 - g2, 0, scl, log = TRUE)
				ll <- lik_group(y1, g1, rep(ls1, length(g1)), rep(nu, length(g1))) +
					lik_group(y2, g2, rep(ls2, length(g2)), rep(nu, length(g2)))
				vals <- c(vals, lmean(ll + prior + w) + log(lt_nodes$w[it]) + base)
			} else {
				m_nodes <- nodes(pool - 1.2 * span - abs(mean(y1) - mean(y2)), pool + 1.2 * span + abs(mean(y1) - mean(y2)), 90L)
				ll <- lik_group(y1, m_nodes$x, rep(ls1, 90L), rep(nu, 90L)) +
					lik_group(y2, m_nodes$x, rep(ls2, 90L), rep(nu, 90L))
				vals <- c(vals, lmean(ll + log(m_nodes$w)) + log(lt_nodes$w[it]) + base)
			}
		}
	}
	lmean(vals)
}

ALL <- !identical(commandArgs(TRUE), "t")
jeu <- list(
	list(n = c(12, 12), d = 0.2, s = c(1, 1)),
	list(n = c(10, 18), d = 1.0, s = c(1, 1)),
	list(n = c(8, 20), d = 0.5, s = c(2.5, 1))
)
cat(R.version.string, "\n")
cat("Colonnes : jeu, loi, hypothese, reference (quadrature), KefiR (IS), ecart, erreur annoncee (rel_se)\n")
for (i in seq_along(jeu)) {
	set.seed(1000 + i)
	y1 <- rnorm(jeu[[i]]$n[1], jeu[[i]]$d, jeu[[i]]$s[1]); y2 <- rnorm(jeu[[i]]$n[2], 0, jeu[[i]]$s[2])
	y <- c(y1, y2); g <- rep(c("a", "b"), jeu[[i]]$n)
	st <- .bay_laws_stats(y, g)
	tab <- .bay_laws_table(st, 0.5, 1.5, 20000L, 10L)
	lab <- list(normal_eq = c("normal", FALSE), normal_het = c("normal", TRUE),
		t_eq = c("t", FALSE), t_het = c("t", TRUE))
	for (j in seq_len(nrow(tab))) {
		v <- lab[[tab$law[j]]]
		if (v[1] != "t" && !ALL) next
		ref <- ref_model(st$y[[1]], st$y[[2]], tab$hypothesis[j] == "H1", as.logical(v[2]), v[1])
		cat(sprintf("%d %-10s %s  ref %10.4f  kefir %10.4f  ecart %8.4f  err %.4f  ess %s\n", i, tab$law[j], tab$hypothesis[j],
			ref, tab$log_ml[j], tab$log_ml[j] - ref, tab$rel_se[j], format(round(tab$ess[j]))))
	}
}

# agent_42 : la seconde forme d'effet ("un seul groupe s'ecarte") sous loi t, k = 3,
# contre une quadrature directe ecrite ici avec dt() et dcauchy() : mu, a_s, log taubar,
# log(nu - 2). Variance commune seulement : a variance par groupe, deux dimensions de
# plus rendent la quadrature trop lente ; cette partie repose sur la calibration (agent_41 A).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
source(file.path(BASE, "dev/mbay/scripts/agent_load_bay.R")); .kefir_bay(BASE)
gl <- function(n) {
	i <- seq_len(n - 1L); b <- i / sqrt(4 * i^2 - 1)
	J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
	e <- eigen(J, symmetric = TRUE); o <- order(e$values)
	list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}
nodes <- function(lo, hi, n) { g <- gl(n); list(x = (hi + lo) / 2 + (hi - lo) / 2 * g$x, w = (hi - lo) / 2 * g$w) }
lmean <- function(l) { t <- max(l); t + log(sum(exp(l - t))) }
ref <- function(st, s, r = 0.5) {
	y_s <- st$y[[s]]
	y_o <- unlist(st$y[-s], use.names = FALSE)
	vals <- numeric(0)
	nu_n <- nodes(-9, 2.6, 30L)
	lt_n <- nodes(-2.6, 2.6, 24L)
	mu_c <- mean(y_o); as_c <- mean(y_s) - mu_c
	span_mu <- 3; span_a <- 4
	mu_n <- nodes(mu_c - span_mu, mu_c + span_mu, 40L)
	a_n <- nodes(as_c - span_a, as_c + span_a, 40L)
	MU <- rep(mu_n$x, 40L); A <- rep(a_n$x, each = 40L)
	W <- log(rep(mu_n$w, 40L)) + log(rep(a_n$w, each = 40L))
	for (iv in seq_along(nu_n$x)) {
		x <- nu_n$x[iv]; nu <- 2 + exp(x)
		for (it in seq_along(lt_n$x)) {
			lt <- lt_n$x[it] + log((nu - 2) / nu)       # fenetre sur l'echelle
			sc <- exp(-lt / 2) * sqrt((nu - 2) / nu)
			ll <- rowSums(dt(outer(MU, y_o, function(m, v) (v - m) / sc), nu, log = TRUE)) +
				rowSums(dt(outer(MU + A, y_s, function(m, v) (v - m) / sc), nu, log = TRUE)) -
				st$N * log(sc)
			pri <- dcauchy(A, 0, r * sqrt(2) * exp(-lt / 2), log = TRUE)
			vals <- c(vals, lmean(ll + pri + W) + log(lt_n$w[it]) + log(nu_n$w[iv]) - exp(x) + x)
		}
	}
	lmean(vals)
}
cat(R.version.string, "\nJeu, groupe s, reference, KefiR, ecart, erreur annoncee (rel_se), ESS\n")
for (i in 1:2) {
	set.seed(4200 + i)
	y <- c(rnorm(9), rnorm(11, c(0, 1)[i]), rnorm(10))
	st <- .bay_laws_stats(y, rep(1:3, c(9, 11, 10)))
	for (s in 1:3) {
		k <- .bay_with_seed(10L + s, .bay_t_ml(st, TRUE, FALSE, 1.5, 0.5, 16000L, s = s))
		rr <- ref(st, s)
		cat(sprintf("%d %d  ref %10.4f  kefir %10.4f  ecart %8.4f  err %.4f  ess %.0f\n", i, s, rr, k$log_ml, k$log_ml - rr, k$rel_se, k$ess))
	}
}

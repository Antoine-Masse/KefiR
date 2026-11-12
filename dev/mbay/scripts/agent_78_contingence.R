# agent_78 : moteur des tableaux de contingence (D30). Ancien (contingencyTableBF, Gunel-Dickey,
# lignes fixees, a = 1) contre nouveau (interaction log-lineaire, -BIC/2 avec n = N), sur les
# memes tables, 1000 par scenario. Sous independance : part des tables a BF > 20 et BF > 3
# (P > .95 et > .75 a cotes egales), intervalle de Wilson. Avec un effet : memes parts.
#   m45_alea : 4 x 5, N = 50, marges tirees au hasard (Dirichlet(1)) a chaque table
#   m45_egal : 4 x 5, N = 50, marges egales
#   m22      : 2 x 2, N = 40, marges egales (20 par groupe), p = 0,4
#   m32      : 3 x 2, 30 par groupe, p = 0,5
#   m32_eff  : 3 x 2, 30 par groupe, p = (0,5 ; 0,5 ; 0,8)
#   m45_eff  : 4 x 5, 25 par groupe, une ligne decalee
args <- commandArgs(TRUE)
sc <- args[1L]
N <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
BASE <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay"
suppressMessages(devtools::load_all(BASE, quiet = TRUE))
rdir <- function(k) { x <- stats::rgamma(k, 1); x / sum(x) }
tab_of <- function(i) {
	set.seed(780000L + i + 10000L * match(sc, c("m45_alea", "m45_egal", "m22", "m32", "m32_eff", "m45_eff")))
	switch(sc,
		m45_alea = matrix(stats::rmultinom(1, 50, outer(rdir(4), rdir(5))), 4, 5),
		m45_egal = matrix(stats::rmultinom(1, 50, outer(rep(1 / 4, 4), rep(1 / 5, 5))), 4, 5),
		m22 = t(vapply(1:2, function(r) as.numeric(stats::rmultinom(1, 20, c(0.6, 0.4))), numeric(2))),
		m32 = t(vapply(1:3, function(r) as.numeric(stats::rmultinom(1, 30, c(0.5, 0.5))), numeric(2))),
		m32_eff = t(vapply(1:3, function(r) as.numeric(stats::rmultinom(1, 30, c(0.5, 0.5) + c(0, 0, 0.3)[r] * c(-1, 1))), numeric(2))),
		m45_eff = t(vapply(1:4, function(r) as.numeric(stats::rmultinom(1, 25, if (r == 4) c(.4, .3, .1, .1, .1) else rep(.2, 5))), numeric(5))))
}
r <- t(vapply(seq_len(N), function(i) {
	tab <- tab_of(i)
	tab <- tab[rowSums(tab) > 0, , drop = FALSE]
	old <- tryCatch(.bay_cont_bf_pkg(tab, 1), error = function(e) NA_real_)
	new <- .bay_cont_lbf(tab)
	chi <- suppressWarnings(stats::chisq.test(tab)$p.value)
	c(old20 = old > log(20), old3 = old > log(3), new20 = new > log(20), new3 = new > log(3), chi = isTRUE(chi < 0.05))
}, numeric(5)))
wil <- function(k, n) { z <- 1.96; ph <- k / n
	c((ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) }
f <- function(j) { k <- sum(r[, j], na.rm = TRUE); n <- sum(!is.na(r[, j])); w <- wil(k, n); sprintf("%d/%d = %.4f [%.4f ; %.4f]", k, n, k / n, w[1], w[2]) }
out <- c(R.version.string, sprintf("agent_78 : contingence, scenario %s, N = %d tables", sc, N),
	paste0("ancien (Gunel-Dickey) BF>20 : ", f(1), " | BF>3 : ", f(2)),
	paste0("nouveau (log-lineaire BIC) BF>20 : ", f(3), " | BF>3 : ", f(4)),
	paste0("khi-deux p<.05 : ", f(5)))
writeLines(out, file.path(BASE, sprintf("dev/mbay/scripts/resultats/agent_78_%s.txt", sc))); cat(out, sep = "\n")

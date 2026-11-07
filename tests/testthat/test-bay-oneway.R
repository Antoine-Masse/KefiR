# Facteur de Bayes a un facteur, variance commune (zs).
# Les messages de .exit() sortent en francais quelle que soit la locale (voir
# R/sys_exit.R) ; ceux de warning(.msg()) suivent la locale : motifs a deux langues.

sim_oneway <- function(n, eff, sds = 1, seed = 1) {
	set.seed(seed)
	k <- length(n)
	g <- factor(rep(seq_len(k), times = n))
	y <- rnorm(sum(n), mean = rep(eff, times = n), sd = rep(rep_len(sds, k), times = n))
	list(y = y, g = g)
}


# ---- zs contre BayesFactor --------------------------------------------------

test_that("zs : anovaBF, equilibre et desequilibre, k = 2, 3, 5, r = 0.5 et 1", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	designs <- list(c(15, 15), c(8, 21), c(12, 12, 12), c(7, 15, 26),
		rep(10, 5), c(5, 9, 14, 20, 6))
	for (i in seq_along(designs)) {
		n <- designs[[i]]
		for (r in c(0.5, 1)) {
			d <- sim_oneway(n, seq(0, 0.7, length.out = length(n)), seed = 100 + i)
			ab <- BayesFactor::extractBF(BayesFactor::anovaBF(y ~ g,
				data.frame(y = d$y, g = d$g), rscaleFixed = r, progress = FALSE),
				logbf = TRUE)
			expect_lt(abs(.bay_bf_oneway_zs(d$y, d$g, r) - ab$bf), 1e-5)
		}
	}
})

test_that("zs : ttestBF(rscale = sqrt(2)/2) a k = 2", {
	skip_if_not_installed("BayesFactor")
	for (n in list(c(20, 20), c(11, 30))) {
		d <- sim_oneway(n, c(0, 0.6), seed = sum(n))
		x1 <- d$y[d$g == 1]
		x2 <- d$y[d$g == 2]
		tb <- BayesFactor::extractBF(BayesFactor::ttestBF(x1, x2, rscale = sqrt(2) / 2),
			logbf = TRUE)
		expect_lt(abs(.bay_bf_oneway_zs(d$y, d$g, 1/2) - tb$bf), 1e-5)
	}
})

# ---- zs : stabilite ---------------------------------------------------------

# Groupes aux residus centres : F est fixe exactement, a l'arrondi pres.
data_with_F <- function(Fval, n = 10, k = 3, seed = 3) {
	set.seed(seed)
	e <- rnorm(n * k)
	g <- factor(rep(seq_len(k), each = n))
	e <- e - ave(e, g)
	e <- e / sqrt(sum(e^2) / (n * k - k))
	mus <- seq(-1, 1, length.out = k)
	mus <- mus * sqrt(Fval * (k - 1) / (n * sum(mus^2)))
	list(y = mus[g] + e, g = g)
}

test_that("zs : F enorme, valeur finie et croissante, jamais NA", {
	Fs <- c(1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14)
	vals <- vapply(Fs, function(Fv) {
		d <- data_with_F(Fv)
		.bay_bf_oneway_zs(d$y, d$g)
	}, numeric(1))
	expect_false(anyNA(vals))
	expect_true(all(is.finite(vals)))
	expect_true(all(diff(vals) > 0))
})

test_that("zs : pente asymptotique en log F, valeur theorique (N - 1 - k)/2", {
	# Pour F >> 1 : (SST/q(g))^{(N-1)/2} apporte F^{(N-1)/2} ; l'integrale en t de
	# exp(-k t/2) jusqu'au pic, situe en g ~ F, retire F^{k/2}. N = 30, k = 3 : 13.
	d1 <- data_with_F(1e10)
	d2 <- data_with_F(1e12)
	slope <- (.bay_bf_oneway_zs(d2$y, d2$g) - .bay_bf_oneway_zs(d1$y, d1$g)) / log(100)
	expect_equal(slope, (30 - 1 - 3) / 2, tolerance = 1e-6)
})

test_that("zs : la ou l'integrale de BayesFactor echoue, nous restons finis", {
	skip_if_not_installed("BayesFactor")
	old <- options(BFprogress = FALSE)
	on.exit(options(old), add = TRUE)
	# anovaBF a F = 1e4 : ancrage, la route est encore commune.
	bf_of <- function(d) {
		out <- NULL
		utils::capture.output(out <- suppressMessages(suppressWarnings(tryCatch(
			BayesFactor::extractBF(BayesFactor::anovaBF(y ~ g, data.frame(y = d$y, g = d$g),
				progress = FALSE), logbf = TRUE)$bf, error = function(e) NA_real_))),
			type = "message")
		out
	}
	d4 <- data_with_F(1e4)
	expect_lt(abs(.bay_bf_oneway_zs(d4$y, d4$g) - bf_of(d4)), 1e-5)
	d6 <- data_with_F(1e6)
	ab6 <- bf_of(d6)
	# BayesFactor 0.9.12 : NA a F = 1e6 (integrate, "roundoff error"). Une version
	# corrigee doit alors concorder avec nous, pas casser ce test.
	v6 <- .bay_bf_oneway_zs(d6$y, d6$g)
	expect_true(is.finite(v6))
	if (!is.na(ab6)) expect_lt(abs(v6 - ab6), 1e-5)
})

test_that("zs : residu nul, Inf avec signal ; y constant et N = k refuses", {
	g <- factor(rep(1:3, each = 4))
	y <- c(1, 2, 4)[g]
	expect_warning(res <- .bay_bf_oneway_zs(y, g), "zero|nulle")
	expect_identical(as.numeric(res), Inf)
	expect_identical(attr(res, "signal"), "zero_residual")
	expect_error(.bay_bf_oneway_zs(rep(2, 12), g), "constant")
	expect_error(.bay_bf_oneway_zs(c(1, 2, 3), factor(1:3)), "N <= k")
})

test_that("zs : F = 0 exact, valeur finie et negative", {
	g <- factor(rep(1:3, each = 5))
	y <- rep(c(-2, -1, 0, 1, 2), 3)
	v <- .bay_bf_oneway_zs(y, g)
	expect_true(is.finite(v))
	expect_lt(v, 0)
})

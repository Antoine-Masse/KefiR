# cor.bay : pendant bayesien de cor.test (dev/corbay/DECISIONS.md, C1 a C9).

d <- local({
	set.seed(1)
	x <- rnorm(30)
	data.frame(x = x, y = 0.5 * x + rnorm(30))
})

test_that("accord numerique avec correlationBF (Pearson, deux moities et bilateral)", {
	r <- cor.bay(d$x, d$y, method = "pearson", verbose = FALSE, sensitivity = FALSE, draws = 1000)
	b2 <- BayesFactor::extractBF(BayesFactor::correlationBF(d$y, d$x, rscale = 1/3), logbf = TRUE)$bf
	expect_equal(r$bf$log_bf10, b2, tolerance = 1e-3)
	bh <- BayesFactor::extractBF(BayesFactor::correlationBF(d$y, d$x, rscale = 1/3, nullInterval = c(0, 1)), logbf = TRUE)$bf
	expect_equal(c(r$bf$log_bf_plus, r$bf$log_bf_minus), bh, tolerance = 1e-10)
	# probabilites : a priori 1/2 aucun lien, 1/4 chaque signe
	w <- c(none = 0.5, pos = 0.25 * exp(bh[1]), neg = 0.25 * exp(bh[2]))
	w <- w / sum(w)
	expect_equal(r$hypotheses$posterior[match(c("none", "pos", "neg"), r$hypotheses$key)], unname(w), tolerance = 1e-10)
	expect_equal(r$prob_link, sum(w[c("pos", "neg")]), tolerance = 1e-10)
	expect_equal(sum(r$hypotheses$posterior), 1)
})

test_that("Kendall : formule de van Doorn et al. 2018, a priori propre, cas limites", {
	expect_equal(stats::integrate(.corbay_prior_tau, -1, 1, kappa = 1/3)$value, 1, tolerance = 1e-8)
	expect_equal(stats::integrate(.corbay_prior_tau, -1, 1, kappa = 1)$value, 1, tolerance = 1e-8)
	# kappa = 1 : a priori uniforme sur rho, donc pi/4 cos(pi tau / 2) sur tau (eq. 10)
	expect_equal(.corbay_prior_tau(0.3, 1), pi / 4 * cos(pi * 0.3 / 2), tolerance = 1e-12)
	# integrale par un autre chemin : grille fine, methode des trapezes
	ts <- .corbay_tstar(d$x, d$y)
	n <- nrow(d)
	g <- seq(-1, 1, length.out = 200001)
	f <- exp(dnorm(ts, 1.5 * sqrt(n) * g, 1, log = TRUE) - dnorm(ts, 0, 1, log = TRUE)) * .corbay_prior_tau(g, 1/3)
	trap <- sum((f[-1] + f[-length(f)]) / 2) * diff(g[1:2])
	expect_equal(.corbay_lbf_kendall(d$x, d$y, 1/3)[["two"]], log(trap), tolerance = 1e-6)
	# n grand : l'integrale reste finie
	set.seed(4)
	x <- rnorm(4000)
	l <- .corbay_lbf_kendall(x, 0.2 * x + rnorm(4000), 1/3)
	expect_true(all(is.finite(l[c("plus", "minus", "two")])))
	expect_gt(l[["plus"]], 10)
})

test_that("unilateral : deux hypotheses, signe teste seulement", {
	r <- cor.bay(d$x, d$y, alternative = "greater", verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_setequal(r$hypotheses$key, c("pos", "none"))
	expect_equal(r$bf$log_bf10, r$bf$log_bf_plus)
	expect_gte(r$estimate$lo, 0)
	r2 <- cor.bay(d$x, d$y, alternative = "less", verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_setequal(r2$hypotheses$key, c("neg", "none"))
	expect_identical(r2$hypotheses$key[1], "none")
	expect_true(any(grepl("autre sens|other way", r2$notes$fr)))
})

test_that("method = auto : Pearson sans valeur a part, Kendall avec une valeur eloignee ou une echelle courte", {
	r <- cor.bay(d$x, d$y, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_identical(r$method, "pearson")
	expect_true(r$auto)
	expect_identical(r$other$method, "kendall")
	y <- d$y; y[30] <- 15
	r2 <- cor.bay(d$x, y, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_identical(r2$method, "kendall")
	expect_true("far" %in% r2$why)
	set.seed(2)
	l <- data.frame(a = sample(1:5, 40, TRUE), b = sample(1:5, 40, TRUE))
	r3 <- cor.bay(l$a, l$b, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_identical(r3$method, "kendall")
	expect_true("scale" %in% r3$why)
	r4 <- cor.bay(d$x, y, method = "spearman", verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_identical(r4$method, "spearman")
	expect_null(r4$other)
	expect_equal(r4$bf$log_bf_plus, .corbay_lbf_pearson(rank(d$x), rank(y), 1/3)[["plus"]])
})

test_that("entrees comme cor.test : vecteurs, formule, noms de colonnes ; refus clairs", {
	a <- cor.bay(d$x, d$y, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	b <- cor.bay(~ x + y, data = d, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	cc <- cor.bay("x", "y", data = d, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_equal(a$hypotheses, b$hypotheses)
	expect_equal(b$hypotheses, cc$hypotheses)
	expect_match(b$label, "x, y")
	dm <- d; dm$y[1:3] <- NA
	r <- cor.bay(~ x + y, data = dm, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_equal(r$n, 27L)
	expect_equal(r$n_removed, 3L)
	expect_error(cor.bay(d$x, rep(1, 30), verbose = FALSE), "constant")
	expect_error(cor.bay(d$x[1:2], d$y[1:2], verbose = FALSE), "3")
	expect_error(cor.bay(d$x, letters[1:30], verbose = FALSE), "num")
	expect_error(cor.bay(d$x, d$y[1:10], verbose = FALSE))
	expect_error(cor.bay(d$x, d$y, kappa = -1, verbose = FALSE), "kappa")
})

test_that("sensibilite : grille 3 x 3, plage et verdict ; graine sans effet sur l'etat global", {
	set.seed(99); before <- .Random.seed
	r <- cor.bay(d$x, d$y, verbose = FALSE, draws = 1000)
	expect_identical(.Random.seed, before)
	expect_equal(dim(r$sensitivity$grid), c(3L, 3L))
	expect_true(r$sensitivity$range_lead[1] <= r$hypotheses$posterior[1] && r$hypotheses$posterior[1] <= r$sensitivity$range_lead[2])
	expect_true(r$sensitivity$verdict %in% c("stable", "crosses", "changes"))
	r2 <- cor.bay(d$x, d$y, verbose = FALSE, draws = 1000)
	expect_identical(r$estimate, r2$estimate)
})

test_that("recit : six etapes, bilan visible, marques, pas de p-value, rejeu identique", {
	txt <- capture.output(r <- cor.bay(d$x, d$y, draws = 1000))
	titres <- grep("^[0-9]+\\) ", txt, value = TRUE)
	expect_equal(as.integer(sub("\\).*", "", titres)), 1:7)
	i7 <- grep("^7\\) Bilan", txt)
	expect_identical(txt[i7 - 1L], "")
	body <- txt[nzchar(txt) & !grepl("^[0-9]+\\) |^   ", txt)]
	expect_true(all(grepl("^\t(==> |--> |Note : |Note: | {4,})", body)))
	expect_false(any(grepl("p[ -]?(value)? ?[=<>] ?[0-9]", txt)))
	expect_identical(capture.output(print(r, story = TRUE)), txt)
	court <- capture.output(print(r))
	expect_true(any(grepl("Bilan de cor.bay|Bilan of cor.bay", court)))
	det <- paste(capture.output(print(r, detail = TRUE)), collapse = "\n")
	expect_match(det, "Ly et al. 2016")
	out <- capture.output(r0 <- cor.bay(d$x, d$y, verbose = FALSE, draws = 1000))
	expect_length(out, 0L)
})

test_that("80 colonnes au plus, tabulation comptee 8, francais et anglais, et width = 60", {
	vis <- function(x) nchar(gsub("\t", "        ", x, fixed = TRUE), type = "width")
	old_loc <- Sys.getlocale("LC_CTYPE"); old_w <- getOption("width")
	on.exit({ Sys.setlocale("LC_CTYPE", old_loc); options(width = old_w) }, add = TRUE)
	set.seed(5)
	lg <- data.frame(longueur_de_la_tige_principale = rnorm(25), masse_seche_des_racines = rnorm(25))
	lg$masse_seche_des_racines[25] <- 9
	cas <- list(
		function() cor.bay(d$x, d$y, draws = 1000),
		function() cor.bay(~ longueur_de_la_tige_principale + masse_seche_des_racines, data = lg, draws = 1000, detail = TRUE),
		function() cor.bay(d$x, d$y, alternative = "less", draws = 1000),
		function() cor.bay(d$x[1:6], d$y[1:6], method = "kendall", draws = 1000))
	for (lc in c("French_France.utf8", "English_United States.utf8")) {
		if (identical(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)), "")) next
		for (w in c(60, 80, 120)) {
			options(width = w)
			for (f in cas) {
				txt <- capture.output(f())
				expect_true(all(vis(txt) <= min(w, 80)), info = paste(lc, w, txt[vis(txt) > min(w, 80)], collapse = " | "))
			}
		}
	}
})

test_that("paire influente (distance de Cook > 1) : Kendall, et le desaccord de Pearson n'est pas une contradiction", {
	set.seed(7)
	x <- rnorm(30); y <- 0.6 * x + 0.8 * rnorm(30)
	x[30] <- 3; y[30] <- -3
	r <- cor.bay(x, y, verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_identical(r$method, "kendall")
	expect_true("influence" %in% r$why)
	expect_true(isTRUE(r$other$explained))
	expect_false(r$contradiction)
	expect_gt(r$checks$influential, 0L)
	# sous 20 paires, non cherchees
	r2 <- cor.bay(x[c(1:15, 30)], y[c(1:15, 30)], verbose = FALSE, sensitivity = FALSE, draws = 1000)
	expect_true(is.na(r2$checks$influential))
})

# m.bay : post-hoc a un facteur (partitions d'egalite, D11).

test_that("partitions : nombres de Bell, a priori de somme 1, aucune difference a 1 - prior_h1", {
	expect_equal(vapply(1:6, function(k) length(.bay_partitions(k)), integer(1)), c(1L, 2L, 5L, 15L, 52L, 203L))
	for (k in 2:5) {
		pa <- .bay_partitions(k)
		lp <- .bay_partition_prior(pa, 0.3)
		expect_equal(sum(exp(lp)), 1)
		expect_equal(exp(lp[vapply(pa, function(p) max(p) == 1L, logical(1))]), 0.7)
	}
	# k = 3, Ewens alpha = 1 : les 4 partitions non nulles ont le meme poids
	lp <- exp(.bay_partition_prior(.bay_partitions(3), 0.5))
	expect_equal(sort(lp), c(rep(0.125, 4), 0.5))
	# k = 4 : {1,2,3}{4} pese 2 fois {1,2}{3}{4} (Ewens : (3-1)! contre (2-1)!)
	pa <- .bay_partitions(4)
	w <- exp(.bay_partition_prior(pa, 0.5))
	key <- vapply(pa, paste, character(1), collapse = "")
	expect_equal(w[key == "1112"] / w[key == "1123"], 2)
})

test_that("probabilites de paires : recalcul direct par lmBF et la loi d'Ewens", {
	set.seed(4)
	d <- data.frame(y = c(rnorm(10, 10), rnorm(10, 11.2), rnorm(10, 10.1)), g = factor(rep(c("a", "b", "c"), each = 10)))
	ph <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, 1L)
	pa <- .bay_partitions(3)
	lb <- vapply(pa, function(p) {
		if (max(p) == 1L) return(0)
		d$h <- factor(p[as.integer(d$g)])
		BayesFactor::extractBF(BayesFactor::lmBF(y ~ h, data = d, rscaleFixed = 0.5, progress = FALSE), logbf = TRUE)$bf
	}, numeric(1))
	pri <- ifelse(vapply(pa, function(p) max(p) == 1L, logical(1)), 0.5, 0.125)
	post <- pri * exp(lb); post <- post / sum(post)
	hand <- c(ab = sum(post[vapply(pa, function(p) p[1] != p[2], logical(1))]),
		ac = sum(post[vapply(pa, function(p) p[1] != p[3], logical(1))]),
		bc = sum(post[vapply(pa, function(p) p[2] != p[3], logical(1))]))
	expect_equal(ph$pairs$p_differ, unname(hand), tolerance = 1e-4)
	# coherence : P(a != c) <= P(a != b) + P(b != c)
	expect_lte(hand[["ac"]], hand[["ab"]] + hand[["bc"]] + 1e-12)
	expect_equal(sum(ph$partitions$posterior), 1)
})

test_that("taille d'effet : delta proche du d de Cohen a grand n, intervalle comparable", {
	set.seed(9)
	n <- 400
	d <- data.frame(y = c(rnorm(n, 0), rnorm(n, 0.5)), g = factor(rep(c("a", "b"), each = n)))
	ph <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 4000, 1L)
	sp <- sqrt((var(d$y[d$g == "a"]) + var(d$y[d$g == "b"])) / 2)
	dc <- diff(tapply(d$y, d$g, mean)) / sp
	expect_lt(abs(ph$pairs$delta - dc), 0.03)
	se <- sqrt(2 / n + dc^2 / (4 * n))
	expect_lt(abs((ph$pairs$delta_high - ph$pairs$delta_low) - 2 * 1.96 * se), 0.03)
	tt <- t.test(y ~ g, data = d, var.equal = TRUE)
	expect_lt(abs(ph$pairs$diff_low - (-tt$conf.int[2])), 0.03)
})

test_that("au-dela du calcul exact : tirage de Gibbs, accord avec l'enumeration, recit", {
	expect_equal(vapply(1:9, .bay_bell, numeric(1)), c(1, 2, 5, 15, 52, 203, 877, 4140, 21147))
	# constante d'Ewens : somme directe sur les partitions non nulles
	pa <- .bay_partitions(5)
	for (a in c(0.5, 1, 2)) {
		lw <- vapply(pa[vapply(pa, max, integer(1)) > 1L], function(p) { sz <- tabulate(p); length(sz) * log(a) + sum(lgamma(sz)) }, numeric(1))
		expect_equal(.bay_ewens_lz(5, a), .bay_log_sum_exp(lw))
	}
	set.seed(21)
	d <- data.frame(y = c(rnorm(10, 0), rnorm(10, 0), rnorm(10, 1), rnorm(10, 1.2)), g = factor(rep(letters[1:4], each = 10)))
	ex <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 1000, 1L)
	expect_identical(ex$method, "exact")
	old <- options(KefiR.bay_exact_max = 3L)
	on.exit(options(old))
	gb <- .bay_posthoc_oneway(d$y, d$g, 0.5, 0.5, 6000, 1L)
	expect_identical(gb$method, "gibbs")
	expect_equal(gb$sweeps, 3000L)
	expect_lt(max(abs(gb$pairs$p_differ - ex$pairs$p_differ)), 0.03)
	expect_lt(max(abs(gb$pairs$p_min - ex$pairs$p_min)), 0.04)
	expect_lt(max(abs(gb$pairs$p_max - ex$pairs$p_max)), 0.04)
	expect_lt(gb$mc_error, 0.02)
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_match(paste(r$story$fr, collapse = " "), "par tirage")
	expect_match(paste(r$story$en, collapse = " "), "by draws")
	expect_true(any(grepl("Gibbs", r$story$detail_fr)))
})

test_that("ecart observe et ecart estime sont nommes", {
	set.seed(4)
	d <- data.frame(y = c(rnorm(10, 10), rnorm(10, 11.2), rnorm(10, 10.1)), g = factor(rep(c("a", "b", "c"), each = 10)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	# l'ecart estime est dit ramene vers 0 ; le bilan donne les moyennes observees
	expect_match(paste(r$story$fr, collapse = " "), "ramen")
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_match(gsub("\\s+", " ", court), "(Moyennes|Médianes) observ|Observed (means|medians)")
})

test_that("m.bay : posthoc = FALSE retire l'etape, le resultat garde les paires sinon", {
	set.seed(4)
	d <- data.frame(y = c(rnorm(10, 10), rnorm(10, 11.2), rnorm(10, 10.1)), g = factor(rep(c("a", "b", "c"), each = 10)))
	r <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_equal(nrow(r$posthoc$pairs), 3L)
	expect_true(any(grepl("Quels groupes", r$story$fr)))
	r0 <- m.bay(y ~ g, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE, posthoc = FALSE)
	expect_null(r0$posthoc)
	expect_true(any(grepl("^6\\) .*non demand", r0$story$fr)))
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_match(court, "Différents :|Semblables :|Non tranchés :|Different:|Alike:|Not settled:")
	expect_false(grepl("ni quels groupes", court))
})

test_that("variances inegales : le recit et la conclusion disent de lire les paires avec reserve", {
	set.seed(12)
	d <- data.frame(dosage = c(rnorm(15, 20, 1), rnorm(15, 21.5, 1), rnorm(15, 20.5, 3)), labo = factor(rep(c("L1", "L2", "L3"), each = 15)))
	r <- m.bay(dosage ~ labo, data = d, verbose = FALSE, draws = 2000, sensitivity = FALSE)
	expect_true(isTRUE(r$posthoc$caution))
	expect_match(paste(r$story$fr, collapse = " "), "avec r\u00e9serve")
	court <- paste(capture.output(print(r)), collapse = " ")
	expect_match(gsub("\\s+", " ", court), "lire avec r\u00e9serve|read them with caution")
})

test_that("echelle extreme : les intervalles restent finis et suivent l'unite de y", {
	set.seed(3)
	y <- c(rnorm(10), rnorm(10, 1))
	g <- rep(1:2, each = 10)
	a <- .bay_posthoc_oneway(y, g, 0.5, 0.5, 2000, 1L)
	b <- .bay_posthoc_oneway(y * 1e200, g, 0.5, 0.5, 2000, 1L)
	expect_false(anyNA(b$pairs[, c("diff_low", "diff_high", "delta_low", "delta_high")]))
	expect_equal(b$pairs$diff / 1e200, a$pairs$diff, tolerance = 1e-6)
	expect_equal(b$pairs$delta, a$pairs$delta, tolerance = 1e-6)
})

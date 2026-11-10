# m.bay : donnees appariees et mesures repetees (D15).

.rm_data <- function(n = 15, k = 3, mu = c(0, 0.5, 1)[seq_len(k)], seed = 1, sd_sub = 2) {
	set.seed(seed)
	d <- data.frame(id = factor(rep(seq_len(n), k)), temps = factor(rep(paste0("t", seq_len(k)), each = n)))
	d$y <- rnorm(n, 0, sd_sub)[d$id] + mu[d$temps] + rnorm(n * k)
	d
}

test_that("apparie : BF egal a ttestBF(paired) sur les differences, effet dz par posterior", {
	d <- .rm_data(n = 20, k = 2, mu = c(0, 0.6), seed = 3)
	r <- m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 2000, sensitivity = FALSE)
	expect_identical(r$plan, "paired")
	ref <- BayesFactor::extractBF(BayesFactor::ttestBF(x = d$y[d$temps == "t2"], y = d$y[d$temps == "t1"],
		paired = TRUE, rscale = 0.5), logbf = TRUE)$bf
	expect_equal(r$effects$log_bf, ref, tolerance = 1e-8)
	expect_equal(r$prob_effect, 1 / (1 + exp(-ref)), tolerance = 1e-8)
	dd <- d$y[d$temps == "t2"] - d$y[d$temps == "t1"]
	expect_lt(abs(r$posthoc$pairs$diff - mean(dd)), 0.1)
	expect_lt(abs(r$posthoc$pairs$delta - mean(dd) / sd(dd)), 0.1)
	expect_true(r$posthoc$pairs$diff_low < r$posthoc$pairs$diff && r$posthoc$pairs$diff < r$posthoc$pairs$diff_high)
	expect_true(any(grepl("ttestBF", r$story$fr)))
	expect_true(any(grepl("dz", r$story$fr)))
})

test_that("apparie : les formes d'appel de m.test donnent le meme resultat", {
	d <- .rm_data(n = 12, k = 2, mu = c(0, 0.8), seed = 5)
	a <- m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)$effects$log_bf
	b <- m.bay(y ~ temps + Error(id / temps), data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)$effects$log_bf
	c <- m.bay(d$y, d$temps, id = d$id, verbose = FALSE, draws = 1000, sensitivity = FALSE)$effects$log_bf
	e <- m.bay(d$y, d$temps, paired = TRUE, verbose = FALSE, draws = 1000, sensitivity = FALSE)$effects$log_bf
	expect_equal(b, a); expect_equal(c, a); expect_equal(e, a)
})

test_that("apparie : paire incomplete retiree et dite, ordre des lignes sans importance avec id", {
	d <- .rm_data(n = 12, k = 2, mu = c(0, 0.8), seed = 6)
	full <- m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	sh <- d[sample(nrow(d)), ]
	expect_equal(m.bay(y ~ temps | id, data = sh, verbose = FALSE, draws = 1000, sensitivity = FALSE)$effects$log_bf, full$effects$log_bf)
	d2 <- d[-1L, ]
	r2 <- m.bay(y ~ temps | id, data = d2, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_equal(r2$n_subjects, 11L)
	expect_true(any(grepl("retir", r2$story$fr)))
})

test_that("mesures repetees : BF contre lmBF(y ~ temps + id, whichRandom) / lmBF(y ~ id), regroupements", {
	d <- .rm_data(n = 15, k = 3, mu = c(0, 0, 1), seed = 7)
	r <- m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 2000, sensitivity = FALSE)
	expect_identical(r$plan, "repeated")
	ref <- BayesFactor::lmBF(y ~ temps + id, data = d, whichRandom = "id", rscaleFixed = 0.5, iterations = 50000, progress = FALSE) /
		BayesFactor::lmBF(y ~ id, data = d, whichRandom = "id", iterations = 50000, progress = FALSE)
	expect_lt(abs(r$effects$log_bf - BayesFactor::extractBF(ref, logbf = TRUE)$bf), 0.1)
	ph <- r$posthoc
	expect_equal(nrow(ph$pairs), 3L)
	expect_equal(sum(ph$partitions$posterior), 1)
	expect_identical(ph$partitions$code[1L], "112")
	expect_true(ph$pairs$p_differ[ph$pairs$group1 == "t1" & ph$pairs$group2 == "t3"] > 0.9)
	expect_true(any(grepl("sujets", r$story$fr)))
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_match(court, "sujets|subjects")
	expect_match(court, "Différents :|Different:")
})

test_that("mesures repetees : l'effet sujet change la reponse par rapport aux groupes independants", {
	d <- .rm_data(n = 15, k = 3, mu = c(0, 0.4, 0.8), seed = 8, sd_sub = 4)
	rm <- m.bay(y ~ temps | id, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	ind <- m.bay(y ~ temps, data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE, laws = "normal_eq")
	expect_gt(rm$prob_effect, ind$prob_effect)
})

test_that("plan en tranches (sujets dans des groupes) : moteur commun, termes du plan", {
	set.seed(9)
	n <- 10
	d <- expand.grid(id = factor(1:(2 * n)), temps = factor(c("t1", "t2", "t3")))
	d$groupe <- factor(ifelse(as.integer(d$id) <= n, "A", "B"))
	d$y <- rnorm(2 * n, 0, 2)[d$id] + 0.8 * (d$temps == "t3") + rnorm(nrow(d))
	r <- m.bay(y ~ groupe * temps + Error(id / temps), data = d, verbose = FALSE, draws = 1000, sensitivity = FALSE)
	expect_identical(r$plan, "repeated")
	expect_setequal(r$effects$term, c("groupe", "temps", "groupe:temps"))
	# Facteurs de Bayes des modeles contre "sujets seuls" : ceux d'anovaBF(whichRandom = "id").
	a <- BayesFactor::anovaBF(y ~ groupe * temps + id, data = d, whichRandom = "id", rscaleFixed = 0.5,
		iterations = 20000, progress = FALSE)
	ref <- BayesFactor::extractBF(a, logbf = TRUE)
	key <- function(x) vapply(strsplit(gsub(" ", "", x), "+", fixed = TRUE), function(v) {
		v <- vapply(strsplit(setdiff(v, "id"), ":", fixed = TRUE), function(w) paste(sort(w), collapse = ":"), character(1))
		paste(sort(v), collapse = "+")
	}, character(1))
	mine <- r$models[-1L, ]
	expect_lt(max(abs(mine$log_bf - ref$bf[match(key(mine$model), key(rownames(ref)))])), 0.1)
	# Post-hoc a plusieurs facteurs (D20) : familles des termes probables, sujets en nuisance.
	expect_identical(r$posthoc$method, "factorial")
})

test_that("mesures repetees : mixte avec repliques refuse clairement", {
	d <- .rm_data(n = 8, k = 2, mu = c(0, 1), seed = 10)
	d <- rbind(d, d)
	expect_error(m.bay(y ~ temps | id, data = d, verbose = FALSE), "mixte|mixed")
})

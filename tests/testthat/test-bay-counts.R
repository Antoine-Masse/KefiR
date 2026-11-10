# m.bay : comptages a un facteur, Poisson et binomiale negative comparees (D17, voie 2).

.cnt <- function(mu, theta = Inf, n = 20, seed = 1) {
	set.seed(seed)
	g <- factor(rep(seq_along(mu), each = n), labels = LETTERS[seq_along(mu)])
	m <- mu[as.integer(g)]
	y <- if (is.finite(theta)) MASS::rnegbin(length(m), mu = m, theta = theta) else stats::rpois(length(m), m)
	data.frame(y = y, g = g)
}

test_that("facteur de Bayes : -BIC/2 des quatre modeles, moyenne sur les deux lois", {
	d <- .cnt(c(4, 4, 8), theta = 2, seed = 51)
	r <- m.bay(y ~ g, data = d, response = "counts", verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)
	expect_identical(r$plan, "counts")
	lm <- c(-BIC(glm(y ~ 1, data = d, family = poisson())) / 2, -BIC(glm(y ~ g, data = d, family = poisson())) / 2,
		-BIC(suppressWarnings(MASS::glm.nb(y ~ 1, data = d))) / 2, -BIC(suppressWarnings(MASS::glm.nb(y ~ g, data = d))) / 2)
	w <- exp(lm - max(lm)); w <- w / sum(w)
	expect_equal(r$prob_effect, w[2] + w[4], tolerance = 1e-8)
	expect_equal(r$laws$posterior[2], w[3] + w[4], tolerance = 1e-8)
})

test_that("surdispersion : la binomiale negative prend la tete et le recit le dit ; Poisson sinon", {
	d <- .cnt(c(5, 5, 5), theta = 1, seed = 2)
	r <- m.bay(y ~ g, data = d, response = "counts", verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)
	expect_gt(r$laws$posterior[2], 0.9)
	expect_gt(r$dispersion, 1.5)
	expect_match(gsub("\\s+", " ", paste(r$story$fr, collapse = " ")), "Changement de stratégie : surdispersion")
	# sous H0 surdispersee, pas de difference affirmee, contrairement a Poisson seul
	expect_lt(r$prob_effect, 0.75)
	pp <- 1 / (1 + exp(-(-BIC(glm(y ~ g, data = d, family = poisson())) + BIC(glm(y ~ 1, data = d, family = poisson()))) / 2))
	expect_lt(r$prob_effect, pp)
	d2 <- .cnt(c(5, 5, 5), seed = 3)
	r2 <- m.bay(y ~ g, data = d2, response = "counts", verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)
	expect_gt(r2$laws$posterior[1], 0.5)
})

test_that("rapports des moyennes proches de glm.nb, paires par la loi unique", {
	d <- .cnt(c(4, 4, 10), theta = 3, n = 60, seed = 4)
	r <- m.bay(y ~ g, data = d, response = "counts", verbose = FALSE, sensitivity = FALSE)
	ref <- exp(stats::coef(suppressWarnings(MASS::glm.nb(y ~ g, data = d)))[c("gB", "gC")])
	pr <- r$posthoc$pairs
	expect_equal(pr$diff[pr$group1 == "A" & pr$group2 == "B"], unname(ref[1]), tolerance = 1e-6)
	expect_equal(pr$diff[pr$group1 == "A" & pr$group2 == "C"], unname(ref[2]), tolerance = 1e-6)
	expect_gt(pr$p_differ[pr$group1 == "A" & pr$group2 == "C"], 0.95)
	expect_equal(sum(r$posthoc$partitions$posterior), 1)
})

test_that("contre-verification : meme sens que BAS (Poisson, a priori robuste) sur des comptages de Poisson", {
	skip_if_not_installed("BAS")
	for (s in 1:3) {
		d <- .cnt(list(c(5, 5, 5), c(4, 4, 7), c(3, 6, 9))[[s]], seed = 10 + s)
		r <- m.bay(y ~ g, data = d, response = "counts", verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)
		b <- BAS::bas.glm(y ~ g, data = d, family = poisson(), betaprior = BAS::robust(), modelprior = BAS::uniform(), force.heredity = TRUE)
		lb_bas <- diff(b$logmarg[order(b$size)])
		expect_identical(sign(r$effects$log_bf), sign(lb_bas))
	}
})

test_that("comptages : refus clairs, 80 colonnes", {
	d <- .cnt(c(3, 3), seed = 5)
	expect_error(m.bay(y ~ g, data = transform(d, y = y + 0.5), response = "counts", verbose = FALSE), "entiers|whole")
	expect_error(m.bay(y ~ g, data = transform(d, y = -y), response = "counts", verbose = FALSE), "entiers|whole")
	d$B <- factor(rep(1:2, 20))
	expect_error(m.bay(y ~ g * B, data = d, response = "counts", verbose = FALSE), "un facteur seulement|one factor only")
	old <- options(width = 80); on.exit(options(old))
	txt <- capture.output(m.bay(y ~ g, data = .cnt(c(4, 4, 8), theta = 2, seed = 51), response = "counts"))
	expect_lte(max(nchar(gsub("\t", "        ", txt), type = "width")), 80)
	expect_true(any(grepl("^7\\) Bilan", txt)))
})

test_that("comptages choisis seuls : bascule dite, nom de la reponse, echelles de notation exclues, console propre", {
	set.seed(7)
	g <- factor(rep(c("a", "b", "c"), each = 15)); nb_poissons <- rpois(45, c(1.5, 4, 9)[g])
	expect_silent(r <- m.bay(nb_poissons, g, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE))
	expect_identical(r$plan, "counts")
	fr <- gsub("\\s+", " ", paste(r$story$fr, collapse = " "))
	expect_match(fr, "Changement de strat\u00e9gie : r\u00e9ponse de comptage")
	expect_match(fr, "nb_poissons")
	expect_false(grepl("comptages de x ", fr))
	# echelle de notation 1 a 5 (pas de zero, 5 valeurs) : voie continue
	note5 <- sample(1:5, 45, TRUE)
	expect_false(identical(m.bay(note5, g, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, draws = 500)$plan, "counts"))
	# moyennes au-dela de 10 : voie continue
	big <- rpois(45, 30)
	expect_false(identical(m.bay(big, g, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, draws = 500)$plan, "counts"))
	# forcer la voie continue
	expect_false(identical(m.bay(nb_poissons, g, response = "continuous", verbose = FALSE, sensitivity = FALSE, posthoc = FALSE, draws = 500)$plan, "counts"))
})

test_that("facteur de Bayes ecrit selon une seule regle", {
	expect_identical(.bay_bf_txt(log(12345), ","), "12\u00a0300")
	expect_identical(.bay_bf_txt(log(2e10), "."), "> 10^9")
	expect_identical(.bay_bf_txt(log(2.5), ","), "2,5")
})

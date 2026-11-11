# m.bay : reponse binaire ou categorielle a un facteur, tableau de contingence (D17, voie 1).

test_that("forme close = contingencyTableBF (multinomiales independantes, lignes fixees)", {
	for (s in 1:4) {
		set.seed(s)
		R <- c(2, 3, 3, 4)[s]; C <- c(2, 2, 3, 3)[s]
		tab <- matrix(rpois(R * C, 12), R, C)
		for (a in c(1, 2, 3)) expect_equal(.bay_cont_lbf(tab, a), .bay_cont_bf_pkg(tab, a), tolerance = 1e-6)
	}
})

test_that("m.bay : BF et probabilite egaux a contingencyTableBF ; 0/1, facteur, texte et logique donnent le meme resultat", {
	set.seed(41)
	d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)))
	d$y01 <- rbinom(90, 1, c(.5, .5, .8)[d$g])
	r <- m.bay(y01 ~ g, data = d, verbose = FALSE, sensitivity = FALSE)
	expect_identical(r$plan, "contingency")
	ref <- .bay_cont_bf_pkg(table(d$g, d$y01), 1)
	expect_equal(r$effects$log_bf, ref, tolerance = 1e-10)
	expect_equal(r$prob_effect, 1 / (1 + exp(-ref)), tolerance = 1e-10)
	d$yf <- factor(d$y01, labels = c("non", "oui"))
	d$yt <- as.character(d$yf)
	d$yl <- d$y01 == 1
	for (v in c("yf", "yt", "yl")) {
		expect_equal(m.bay(stats::as.formula(paste(v, "~ g")), data = d, verbose = FALSE, sensitivity = FALSE, posthoc = FALSE)$effects$log_bf, ref, tolerance = 1e-10)
	}
	expect_true(any(grepl("contingencyTableBF", r$story$fr)))
	expect_true(any(grepl("khi-deux", r$story$fr)))
})

test_that("post-hoc : BF d'un regroupement = contingencyTableBF du tableau fusionne ; paires par la loi unique", {
	set.seed(5)
	tab <- matrix(c(12, 18, 12, 18, 6, 24), 3, 2, byrow = TRUE, dimnames = list(c("a", "b", "c"), c("n", "o")))
	ph <- .bay_cont_posthoc(tab, 1, 0.5, 4000, 1L)
	pa <- .bay_partitions(3)
	lb <- vapply(pa, function(p) if (max(p) == 1L) 0 else .bay_cont_bf_pkg(rowsum(tab, p), 1), numeric(1))
	pri <- ifelse(vapply(pa, max, integer(1)) == 1L, 0.5, 0.125)
	post <- pri * exp(lb); post <- post / sum(post)
	hand <- c(sum(post[vapply(pa, function(p) p[1] != p[2], logical(1))]), sum(post[vapply(pa, function(p) p[1] != p[3], logical(1))]),
		sum(post[vapply(pa, function(p) p[2] != p[3], logical(1))]))
	expect_equal(ph$pairs$p_differ, hand, tolerance = 1e-8)
	expect_identical(ph$delta_unit, "logor")
})

test_that("taille d'effet : difference de proportions et log du rapport de cotes proches de glm a grand n", {
	set.seed(7)
	n <- 2000
	d <- data.frame(g = factor(rep(c("a", "b"), each = n)))
	d$y <- rbinom(2 * n, 1, c(.3, .45)[d$g])
	r <- m.bay(y ~ g, data = d, verbose = FALSE, sensitivity = FALSE)
	obs <- diff(tapply(d$y, d$g, mean))
	expect_lt(abs(r$posthoc$pairs$diff - obs), 0.005)
	lor <- unname(stats::coef(stats::glm(y ~ g, data = d, family = stats::binomial()))[2])
	expect_lt(abs(r$posthoc$pairs$delta - lor), 0.02)
})

test_that("plus de deux modalites : ecart total des proportions, sans delta ; recit et bilan lisibles", {
	set.seed(43)
	d <- data.frame(p = factor(rep(c("P1", "P2"), each = 40)))
	d$c <- factor(c(sample(c("blanc", "rose", "rouge"), 40, TRUE, c(.5, .3, .2)), sample(c("blanc", "rose", "rouge"), 40, TRUE, c(.2, .3, .5))))
	r <- m.bay(c ~ p, data = d, verbose = FALSE)
	expect_identical(r$posthoc$delta_unit, "tvd")
	tvd <- 0.5 * sum(abs(prop.table(table(d$p, d$c), 1)[2, ] - prop.table(table(d$p, d$c), 1)[1, ]))
	expect_lt(abs(r$posthoc$pairs$diff - tvd), 0.06)
	expect_true(is.na(r$posthoc$pairs$delta))
	court <- paste(capture.output(print(r)), collapse = "\n")
	expect_false(grepl("NA", court))
	expect_false(grepl("% %", court, fixed = TRUE))
})

test_that("reponse binaire avec plusieurs facteurs : voie logistique (D25), plus de deux modalites refusees", {
	d <- data.frame(y = rbinom(40, 1, .5), A = factor(rep(1:2, 20)), B = factor(rep(1:2, each = 20)))
	expect_identical(m.bay(y ~ A * B, data = d, verbose = FALSE, sensitivity = FALSE)$plan, "discrete_fact")
	d$z <- factor(sample(c("a", "b", "c"), 40, replace = TRUE))
	expect_error(m.bay(z ~ A * B, data = d, verbose = FALSE), "un facteur seulement|one factor only")
})

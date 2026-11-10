# m.bay a plusieurs facteurs : espace hereditaire, loi HOP, BayesFactor, bilan.

test_that("espace hereditaire : 5 modeles a 2 facteurs, 19 a 3, parents exiges", {
	expect_length(.bay_fact_space(c("A", "B", "A:B"))$models, 5L)
	expect_length(.bay_fact_space(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))$models, 19L)
	expect_length(.bay_fact_space(c("A", "B", "C"))$models, 8L)
	expect_error(.bay_fact_space(c("A", "A:B")), "ordre inf|lower-order")
	expect_error(.bay_fact_space(c("A", "B", "C", "A:B:C")), "ordre inf|lower-order")
})

test_that("loi HOP : somme 1, valeurs exactes a 2 facteurs, interaction a 1/2 sachant ses parents", {
	sp <- .bay_fact_space(c("A", "B", "A:B"))
	for (ty in c("hop", "hip", "unif")) expect_equal(sum(exp(.bay_fact_prior(sp, ty))), 1)
	p <- exp(.bay_fact_prior(sp, "hop"))
	keys <- vapply(sp$models, function(m) paste(sp$keys[m], collapse = "+"), "")
	expect_equal(unname(p[match(c("", "A", "B", "A+B", "A+B+A:B"), keys)]), c(1/3, 1/6, 1/6, 1/6, 1/6))
	sp3 <- .bay_fact_space(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C"))
	expect_equal(sum(exp(.bay_fact_prior(sp3, "hop"))), 1)
	inc <- .bay_fact_terms(sp3, rep(0, 19), .bay_fact_prior(sp3, "hop"))
	expect_equal(inc$p_h1[inc$term == "A"], 1/2)
	expect_equal(inc$p_h1[inc$term == "A:B"], 1/6)
	expect_equal(inc$p_h1[inc$term == "A:B:C"], 1/32)
})

test_that("a un facteur, BayesFactor et la voie interne normal_eq donnent le meme log BF", {
	skip_if_not_installed("BayesFactor")
	set.seed(5)
	for (k in c(2L, 3L, 5L)) {
		g <- factor(rep(seq_len(k), times = c(7, 9, 8, 6, 10)[seq_len(k)]))
		y <- rnorm(length(g), as.integer(g) / 3)
		for (r in c(0.25, 0.5, 1)) {
			bf <- BayesFactor::extractBF(BayesFactor::lmBF(y ~ g, data = data.frame(y = y, g = g),
				rscaleFixed = r, progress = FALSE), logbf = TRUE)$bf
			expect_equal(bf, .bay_bf_oneway_zs(y, g, r), tolerance = 1e-4)
		}
	}
})

set.seed(8)
d2 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 1:6)
d2$y <- rnorm(nrow(d2), 1.5 * (d2$A == "a2"))
gof <- function(...) m.bay(..., verbose = FALSE, sensitivity = FALSE, draws = 2000L)

test_that("plan 2x3 : resultat complet, reproductible, sans p-value", {
	skip_if_not_installed("BayesFactor")
	r <- gof(y ~ A * B, data = d2)
	expect_s3_class(r, "mbay")
	expect_identical(r$plan, "factorial")
	expect_equal(r$effects$term, c("A", "B", "A:B"))
	expect_true(all(r$effects$p_h1 > 0 & r$effects$p_h1 < 1))
	expect_equal(r$effects$prior_h1, c(1/2, 1/2, 1/6))
	expect_equal(sum(r$hypotheses$posterior), 1)
	expect_false(any(grepl("p.value|pvalue", c(names(r), names(r$effects)))))
	expect_identical(gof(y ~ A * B, data = d2)$effects$p_h1, r$effects$p_h1)
	expect_true(r$effects$p_h1[1L] > 0.5)
	# Les memes donnees sous d'autres formes.
	expect_equal(gof(d2$y, d2[c("A", "B")])$n_models, 4L)
	expect_equal(gof(y ~ A + B, data = d2)$effects$term, c("A", "B"))
})

test_that("bilan multi-facteurs : court, loi annoncee, sensibilite, detail", {
	skip_if_not_installed("BayesFactor")
	r <- m.bay(y ~ A * B, data = d2, verbose = FALSE, draws = 2000L)
	court <- capture.output(print(r))
	expect_lte(length(court), 40L)
	txt <- gsub("\\s+", " ", paste(court, collapse = " "))
	expect_match(txt, "variance commune|common variance")
	expect_match(txt, "p-value")
	expect_match(txt, "seul ou en interaction|alone or in interaction")
	long <- paste(capture.output(print(r, detail = TRUE)), collapse = " ")
	expect_match(long, "HOP")
	expect_match(long, "BayesFactor")
	expect_equal(dim(r$sensitivity$p_terms), c(3L, 9L))
})

test_that("refus et avis a plusieurs facteurs", {
	skip_if_not_installed("BayesFactor")
	e <- d2[!(d2$A == "a1" & d2$B == "b3"), ]
	expect_error(gof(y ~ A * B, data = e), "Case vide|Empty cell")
	expect_silent(r <- gof(y ~ A + B, data = e))
	r <- gof(y ~ A * B, data = d2, prior_h1 = 0.3, laws = "t_eq")
	expect_true(any(grepl("prior_h1", r$notes$en)))
	expect_true(any(grepl("laws", r$notes$en)))
	expect_false(any(grepl("laws", gof(y ~ A * B, data = d2)$notes$en)))
})

test_that("effet enorme : pas de NA (integrate de lmBF en echec, agent_50)", {
	skip_if_not_installed("BayesFactor")
	set.seed(13)
	e <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), r = 1:5)
	e$y <- rnorm(nrow(e), 150 * c(-1, 0, 1)[as.integer(e$B)])
	r <- gof(y ~ A * B, data = e)
	expect_false(anyNA(r$effects$p_h1))
	expect_true(r$effects$p_h1[2L] > 0.999)
})

test_that("relecture du 9 novembre : noms, facteurs confondus, graine, ecart-type, avis", {
	skip_if_not_installed("BayesFactor")
	# 1. Noms non syntaxiques.
	d3 <- d2
	names(d3)[1:2] <- c("my A", "B-2")
	r1 <- gof(y ~ `my A` * `B-2`, data = d3)
	expect_equal(r1$effects$p_h1, gof(y ~ A * B, data = d2)$effects$p_h1)
	# 2. Facteurs confondus ou emboites : refus.
	da <- d2
	da$C <- factor(ifelse(da$A == "a1", "c1", "c2"))
	expect_error(gof(y ~ A + C, data = da), "confondus|confounded")
	da$N <- factor(paste(da$A, da$B))
	expect_error(gof(y ~ A + N, data = da), "confondus|confounded")
	# 6. La reponse centrale ne depend pas de sensitivity.
	expect_equal(m.bay(y ~ A * B, data = d2, verbose = FALSE, sensitivity = TRUE)$effects$p_h1,
		gof(y ~ A * B, data = d2)$effects$p_h1)
	# 4. Ecart-type au sein des cases.
	expect_equal(gof(y ~ A * B, data = d2)$sd_within, sigma(lm(y ~ A * B, data = d2)))
	expect_match(paste(capture.output(print(gof(y ~ A * B, data = d2))), collapse = " "), "au sein des cases|Within-cell")
	# 10. Aucune interaction testee : dit.
	expect_true(any(grepl("No interaction", gof(y ~ A + B, data = d2)$notes$en)))
	expect_false(any(grepl("No interaction", gof(y ~ A * B, data = d2)$notes$en)))
	# 3. Precision mesuree par deux graines, en points de probabilite.
	expect_true(all(gof(y ~ A * B, data = d2)$effects$error_pct < 5))
})

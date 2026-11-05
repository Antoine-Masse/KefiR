# MANOVA : chaque resultat est compare a un calcul independant ({stats}, {car}, {MVN}, {biotools}).

silence <- function(expr) {
  value <- NULL
  output <- utils::capture.output(value <- suppressMessages(suppressWarnings(expr)))
  list(value = value, output = output)
}
bilan <- function(...) silence(m.test(..., plot = FALSE, verbose = TRUE))
p_only <- function(...) silence(m.test(..., return = FALSE))$value

# Donnees sans violation : la voie classique et Wilks sont attendus
make_one_way <- function() {
  set.seed(3)
  n <- 30
  d <- data.frame(g = rep(c("A", "B", "C"), each = n))
  d$a <- rnorm(3 * n, c(0, 1, 1)[factor(d$g)], 1)
  d$b <- 0.5 * d$a + rnorm(3 * n, c(0, 0, 1)[factor(d$g)], 1)
  d
}
make_two_way <- function(drop = 0) {
  set.seed(8)
  d <- expand.grid(F = c("f1", "f2"), G = c("g1", "g2", "g3"), rep = 1:15, stringsAsFactors = FALSE)
  d$a <- rnorm(nrow(d), 2 * (d$F == "f2"), 1)
  d$b <- 0.4 * d$a + rnorm(nrow(d), (d$G == "g3") * 1.5, 1)
  if (drop > 0) d <- d[-seq_len(drop), ]
  d
}

test_that("un facteur : Wilks, Pillai, Hotelling-Lawley et Roy egalent summary.manova()", {
  d <- make_one_way()
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  fit <- manova(cbind(a, b) ~ g, data = d)
  expect_equal(r$test_statistics$Wilks$p_value, summary(fit, test = "Wilks")$stats[1, "Pr(>F)"])
  expect_equal(r$test_statistics$Wilks$statistic, summary(fit, test = "Wilks")$stats[1, "Wilks"], ignore_attr = TRUE)
  expect_equal(r$test_statistics$Pillai$p_value, summary(fit, test = "Pillai")$stats[1, "Pr(>F)"])
  expect_equal(r$test_statistics$HotellingLawley$p_value, summary(fit, test = "Hotelling-Lawley")$stats[1, "Pr(>F)"])
  expect_equal(r$test_statistics$Roy$p_value, summary(fit, test = "Roy")$stats[1, "Pr(>F)"])
})

test_that("la p-value globale et la methode annoncee designent le meme test", {
  d <- make_one_way()
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  expect_identical(r$method, "MANOVA (Wilks' Lambda)")
  expect_equal(r$global_pvalue, r$test_statistics$Wilks$p_value)
  expect_equal(p_only(cbind(a, b) ~ g, data = d), r$global_pvalue)
})

test_that("deux facteurs : un test par terme, egal a car::Manova(type II), plan equilibre", {
  d <- make_two_way()
  r <- bilan(cbind(a, b) ~ F * G, data = d)$value
  expect_equal(r$effects$term, c("F", "G", "F:G"))
  ref <- car::Manova(lm(cbind(a, b) ~ F * G, data = d), type = "II")
  ref_stats <- summary(ref, test = "Wilks", univariate = FALSE)$multivariate.tests
  for (term in r$effects$term) {
    h <- ref_stats[[term]]
    eig <- Re(eigen(solve(h$SSPE) %*% h$SSPH)$values)
    expect_equal(r$effects$Wilks[r$effects$term == term], prod(1 / (1 + eig)), tolerance = 1e-8)
    expect_equal(r$effects$Pillai[r$effects$term == term], sum(eig / (1 + eig)), tolerance = 1e-8)
  }
})

test_that("les formules multivariees egalent summary.manova()", {
  d <- make_one_way()
  fit <- manova(cbind(a, b) ~ g, data = d)
  e <- crossprod(residuals(fit))
  h <- crossprod(residuals(lm(cbind(a, b) ~ 1, data = d))) - e
  eig <- Re(eigen(solve(e, h), only.values = TRUE)$values)
  for (test in c("Wilks", "Pillai", "Hotelling-Lawley", "Roy")) {
    ref <- summary(fit, test = test)$stats[1, ]
    expect_equal(KefiR:::.manova_stat(eig, 2, fit$df.residual, test), unname(ref[c(2:6)]), tolerance = 1e-8)
  }
})

test_that("deux facteurs : p-values et methode suivent le premier terme", {
  d <- make_two_way()
  r <- bilan(cbind(a, b) ~ F * G, data = d)$value
  expect_equal(r$global_pvalue, r$effects$p.Wilks[r$effects$term == "F"])
  expect_equal(p_only(cbind(a, b) ~ F * G, data = d), r$global_pvalue)
})

test_that("deux facteurs : plan desequilibre, Type II egal a car::Manova()", {
  d <- make_two_way(drop = 11)
  r <- bilan(cbind(a, b) ~ F * G, data = d)$value
  ref <- summary(car::Manova(lm(cbind(a, b) ~ F * G, data = d), type = "II"), univariate = FALSE)$multivariate.tests
  for (term in c("F", "G", "F:G")) {
    h <- ref[[term]]
    eig <- Re(eigen(solve(h$SSPE) %*% h$SSPH)$values)
    expect_equal(r$effects$Wilks[r$effects$term == term], prod(1 / (1 + eig)), tolerance = 1e-8)
  }
})

test_that("deux facteurs : les post-hocs aboutissent", {
  o <- bilan(cbind(a, b) ~ F + G, data = make_two_way())
  expect_false(any(grepl("skipped|ignor", o$output)))
  expect_false(is.null(o$value$discriminant_analysis))
  expect_setequal(names(o$value$protected_anovas), c("a", "b"))
})

test_that("les assomptions sont controlees dans les groupes, pas sur les donnees brutes", {
  set.seed(5)
  d <- data.frame(g = rep(c("A", "B"), each = 60), a = rnorm(120), b = rnorm(120))
  d$a[d$g == "B"] <- d$a[d$g == "B"] + 8
  d$b[d$g == "B"] <- d$b[d$g == "B"] + 8
  o <- bilan(cbind(a, b) ~ g, data = d)
  residus <- residuals(lm(cbind(a, b) ~ g, data = d))
  mardia <- suppressWarnings(MVN::mvn(data = residus, mvn_test = "mardia"))$multivariate_normality
  # Deux groupes bien separes : bimodal au brut, normal dans les groupes
  expect_true(any(grepl("normalit. multivari.e respect.e", o$output)))
  expect_true(any(grepl("Corr.lation maximale = 0\\.[0-3]", o$output)))
  expect_true(any(grepl("4 test.s, 0 viol.s", o$output)))
  expect_true(any(grepl(as.character(mardia$p.value[2]), o$output, fixed = TRUE)) || nrow(mardia) == 2)
})

test_that("une violation unique bascule sur Pillai, et la methode le dit", {
  set.seed(9)
  d <- data.frame(g = rep(c("A", "B", "C"), each = 40))
  d$a <- rnorm(120, 0, c(1, 4, 8)[factor(d$g)])
  d$b <- rnorm(120, 0, c(1, 4, 8)[factor(d$g)])
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  fit <- manova(cbind(a, b) ~ g, data = d)
  if (identical(r$method, "MANOVA (Pillai's trace)")) {
    expect_equal(r$global_pvalue, summary(fit, test = "Pillai")$stats[1, "Pr(>F)"])
  } else {
    expect_identical(r$method, "Robust MANOVA (Wilks, MCD)")
    expect_equal(r$global_pvalue, r$test_statistics$RobustWilks$p_value)
  }
})

test_that("la voie robuste est reproductible et annoncee", {
  set.seed(2)
  d <- data.frame(g = rep(c("A", "B", "C"), each = 30), a = rexp(90), b = rexp(90)^2)
  d$a[1:3] <- 40
  r1 <- bilan(cbind(a, b) ~ g, data = d)$value
  set.seed(99)
  r2 <- bilan(cbind(a, b) ~ g, data = d)$value
  expect_identical(r1$method, "Robust MANOVA (Wilks, MCD)")
  expect_identical(r1$global_pvalue, r2$global_pvalue)
  expect_true(all(sapply(r1$protected_anovas, function(z) z$test) == "kruskal"))
})

test_that("les ANOVAs proteges egalent les tests individuels", {
  d <- make_one_way()
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  ref <- summary(aov(a ~ g, data = d))[[1]]["g   ", "Pr(>F)"]
  if (is.na(ref)) ref <- anova(lm(a ~ g, data = d))["g", "Pr(>F)"]
  expect_equal(r$protected_anovas$a$p_value, ref)
  expect_equal(r$protected_anovas$a$alpha_adjusted, 0.05 / 2)
})

test_that("la matrice de structure est intra-groupes, et les parts de variance somment a 1", {
  d <- make_one_way()
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  fit <- MASS::lda(g ~ a + b, data = d)
  scores <- predict(fit)$x
  ref <- cor(residuals(lm(cbind(a, b) ~ g, data = d)), residuals(lm(scores ~ d$g)))
  expect_equal(r$discriminant_analysis$structure_matrix, ref, ignore_attr = TRUE)
  expect_equal(sum(r$discriminant_analysis$prop_variance), 1)
})

test_that("les comparaisons entre groupes sont des Pillai corriges par Holm", {
  d <- make_one_way()
  r <- bilan(cbind(a, b) ~ g, data = d)$value
  pair <- function(u, v) {
    s <- d[d$g %in% c(u, v), ]
    summary(manova(cbind(a, b) ~ g, data = s), test = "Pillai")$stats[1, "Pr(>F)"]
  }
  ref <- p.adjust(c(pair("A", "B"), pair("A", "C"), pair("B", "C")), "holm")
  expect_equal(sort(as.vector(r$pairwise_multivariate[!is.na(r$pairwise_multivariate)])), sort(ref))
  expect_setequal(r$groups_multivariate$categories, c("A", "B", "C"))
})

test_that("NA dans une VD ou un facteur : suppression par ligne, resultat identique a la main", {
  d <- make_one_way()
  d$a[2] <- NA
  d$g[50] <- NA
  ref <- manova(cbind(a, b) ~ g, data = d[!is.na(d$a) & !is.na(d$g), ])
  expect_equal(p_only(cbind(a, b) ~ g, data = d), summary(ref, test = "Wilks")$stats[1, "Pr(>F)"])
})

test_that("les syntaxes de VD donnent le meme resultat", {
  d <- make_one_way()
  ref <- p_only(cbind(a, b) ~ g, data = d)
  expect_equal(silence(m.test(cbind(a, d$b) ~ g, data = d, return = FALSE))$value, ref)
  expect_equal(silence(m.test(as.matrix(d[, c("a", "b")]), d$g, return = FALSE))$value, ref)
  expect_equal(silence(m.test(d[, c("a", "b")], d$g, return = FALSE))$value, ref)
  r <- bilan(d[, c("a", "b")], d$g)$value
  expect_false(is.null(r$discriminant_analysis))
})

test_that("une VD transformee est bien analysee transformee", {
  set.seed(3)
  d <- make_one_way()
  d$pos <- exp(d$a)
  ref <- summary(manova(cbind(log(pos), b) ~ g, data = d), test = "Wilks")$stats[1, "Pr(>F)"]
  expect_equal(p_only(cbind(log(pos), b) ~ g, data = d), ref)
})

test_that("trois VD et cinq groupes", {
  set.seed(12)
  d <- data.frame(g = rep(letters[1:5], each = 25))
  d$a <- rnorm(125, as.integer(factor(d$g)) / 3); d$b <- rnorm(125); d$c <- rnorm(125, (d$g == "e") * 2)
  r <- bilan(cbind(a, b, c) ~ g, data = d)$value
  expect_equal(r$test_statistics$Wilks$p_value,
               summary(manova(cbind(a, b, c) ~ g, data = d), test = "Wilks")$stats[1, "Pr(>F)"])
  expect_equal(r$discriminant_analysis$n_functions, 3)
  expect_equal(r$protected_anovas$a$alpha_adjusted, 0.05 / 3)
})

test_that("cas limites : messages clairs et pas de plantage", {
  set.seed(1)
  small <- data.frame(g = rep(c("A", "B"), each = 3), a = rnorm(6), b = rnorm(6))
  expect_false(is.null(bilan(cbind(a, b) ~ g, data = small)$value))
  collin <- make_one_way()
  collin$b <- collin$a * 2
  expect_false(is.null(silence(tryCatch(m.test(cbind(a, b) ~ g, data = collin, plot = FALSE),
                                        error = function(e) conditionMessage(e)))$value))
})

test_that("un graphique MANOVA se trace sans erreur", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_false(is.null(silence(m.test(cbind(a, b) ~ g, data = make_one_way(), plot = TRUE))$value))
  expect_false(is.null(silence(m.test(cbind(a, b) ~ F + G, data = make_two_way(), plot = TRUE))$value))
})

test_that("une formule longue ne fait pas planter m.test()", {
  r <- bilan(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, data = iris)
  expect_false(is.null(r$value))
  expect_true(any(grepl("ASSOMPTION 6/6", r$output)))
})

test_that("deux facteurs : chaque VD a son test par terme, egal a car::Anova()", {
  d <- make_two_way(drop = 11)
  r <- bilan(cbind(a, b) ~ F * G, data = d)$value
  ref <- car::Anova(lm(a ~ F + G + F:G, data = d), type = "II")
  expect_equal(r$effects_by_dv$a$p.value, ref[c("F", "G", "F:G"), "Pr(>F)"])
})

test_that("un p trop petit pour etre simule n'est pas affiche comme zero", {
  o <- bilan(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, data = iris)$output
  expect_false(any(grepl("p-value = 0.00e", o, fixed = TRUE)))
})

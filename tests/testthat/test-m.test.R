# Chaque p-value de m.test() est comparee a un calcul independant en R de base ou {car}.

silence <- function(expr) {
  value <- NULL
  output <- utils::capture.output(value <- suppressMessages(suppressWarnings(expr)))
  list(value = value, output = output)
}
bilan <- function(...) silence(m.test(..., plot = FALSE, verbose = TRUE))
p_only <- function(...) silence(m.test(..., return = FALSE))$value
type3 <- function(formula, data, factor_name = "F") {
  contrasts_list <- stats::setNames(list("contr.sum"), factor_name)
  car::Anova(lm(formula, data = data, contrasts = contrasts_list), type = "III")
}

# Graines retenues pour rester sur la voie parametrique
make_ancova <- function() {
  set.seed(2)
  d <- data.frame(F = rep(c("A", "B"), each = 20), x = runif(40, 200, 220))
  d$y <- 10 + 2 * (d$F == "B") + 0.05 * d$x + rnorm(40, 0, 1)
  d
}
make_slopes <- function() {
  set.seed(1)
  d <- data.frame(F = rep(c("A", "B", "C"), each = 30), x = runif(90, 100, 120))
  d$y <- 50 + c(0, 4, 8)[factor(d$F)] + c(1, 3, -2)[factor(d$F)] * (d$x - 110) + rnorm(90, 0, 2)
  d
}

test_that("une interaction non significative quitte reellement le modele", {
  d <- make_ancova()
  b <- bilan(y ~ F * x, data = d)
  ref <- type3(y ~ F + x, d)["F", "Pr(>F)"]
  expect_true(any(grepl("Formule : y ~ F \\+ x", b$output)))
  expect_false(any(grepl("^F:x", b$output)))
  expect_equal(b$value$global_pvalue, ref)
  expect_equal(b$value$p.value["A", "B"], ref)
  expect_equal(p_only(y ~ F * x, data = d), ref)
})

test_that("une interaction significative est gardee, covariable centree", {
  d <- make_slopes()
  b <- bilan(y ~ F * x, data = d)
  centered <- d
  centered$x <- d$x - mean(d$x)
  expect_true(any(grepl("^F:x", b$output)))
  expect_true(any(grepl("centr.e\\(s\\) sur leur moyenne", b$output)))
  expect_equal(b$value$global_pvalue, type3(y ~ F * x, centered)["F", "Pr(>F)"])
})

test_that("les moyennes ajustees ne dependent pas du centrage", {
  d <- make_slopes()
  b <- bilan(y ~ F * x, data = d)
  ref <- summary(emmeans::emmeans(lm(y ~ F * x, data = d), "F"))
  expect_equal(sort(b$value$descriptive_stats$adjusted_mean), sort(ref$emmean))
})

test_that("la p-value globale est celle du premier facteur", {
  set.seed(1)
  d <- data.frame(F = rep(c("A", "B", "C"), each = 15), x = runif(45, 20, 30))
  d$y <- c(0, 1, 2)[factor(d$F)] + 0.2 * d$x + rnorm(45, 0, 1)
  ref <- anova(lm(y ~ x + F, data = d))["F", "Pr(>F)"]
  expect_equal(p_only(y ~ x + F, data = d), ref)
  expect_equal(bilan(y ~ x + F, data = d)$value$global_pvalue, ref)
})

test_that("ANOVA a deux facteurs : p globale et $effects suivent le tableau affiche", {
  set.seed(2)
  d <- data.frame(F = rep(c("A", "B"), each = 20), G = rep(c("u", "v"), 20))
  d$y <- 10 + (d$F == "B") * 2 + rnorm(40)
  ref <- anova(lm(y ~ F * G, data = d))
  b <- bilan(y ~ F * G, data = d)
  expect_equal(p_only(y ~ F * G, data = d), ref["F", "Pr(>F)"])
  expect_equal(b$value$effects$term, c("F", "G", "F:G"))
  expect_equal(b$value$effects$p.value, ref[c("F", "G", "F:G"), "Pr(>F)"])
})

test_that("une reponse transformee n'est transformee qu'une fois", {
  set.seed(2)
  d <- data.frame(F = rep(c("A", "B"), each = 20), x = runif(40, 200, 220))
  d$y <- exp(1 + 0.3 * (d$F == "B") + 0.005 * d$x + rnorm(40, 0, 0.2))
  ref <- type3(log(y) ~ F + x, d)["F", "Pr(>F)"]
  expect_equal(p_only(log(y) ~ F + x, data = d), ref)
  expect_equal(p_only(I(log(y)) ~ F + x, data = d), ref)
  d$G <- rep(c("u", "v"), 20)
  expect_equal(p_only(log(y) ~ F * G, data = d), anova(lm(log(y) ~ F * G, data = d))["F", "Pr(>F)"])
  expect_equal(p_only(log(y) ~ F, data = d), t.test(log(y) ~ F, data = d, var.equal = TRUE)$p.value)
})

test_that("$global_pvalue et return=FALSE donnent le test global, pas une comparaison par paires", {
  ref <- oneway.test(Sepal.Length ~ Species, data = iris)$p.value
  expect_equal(bilan(Sepal.Length ~ Species, data = iris, boot = FALSE)$value$global_pvalue, ref)
  set.seed(5)
  d <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"), rep = 1:12, stringsAsFactors = FALSE)
  d$y <- rpois(nrow(d), 3)
  d <- d[-c(1:5, 20:21), ]
  r <- bilan(y ~ A * B, data = d, boot = FALSE)$value
  expect_equal(r$global_pvalue, r$effects$p.value[r$effects$term == "A"])
  expect_equal(r$global_pvalue, p_only(y ~ A * B, data = d))
})

test_that("un groupe constant ne fait pas planter l'analyse", {
  set.seed(7)
  d <- data.frame(y = rnorm(60, 10, 2), g = rep(c("a", "b", "c"), 20))
  d$y[d$g == "a"] <- 0
  expect_false(is.null(bilan(y ~ g, data = d, boot = FALSE)$value))
  two <- d[d$g != "c", ]
  expect_false(is.null(bilan(y ~ g, data = two, boot = FALSE)$value))
})

test_that("un test a un facteur renvoie la p-value du test annonce", {
  expect_equal(p_only(Sepal.Length ~ Species, data = iris),
               oneway.test(Sepal.Length ~ Species, data = iris)$p.value)
  b <- bilan(Sepal.Length ~ Species, data = iris, boot = FALSE)
  expect_identical(b$value$method, "Welch ANOVA")
  expect_true(any(grepl("F = .*ddl = ", b$output)))
})

test_that("$method nomme le test reel a deux groupes", {
  set.seed(5)
  d <- data.frame(y = c(rexp(40, 1), rexp(40, 0.2)), g = rep(c("A", "B"), each = 40))
  r <- silence(m.test(d$y, d$g, plot = FALSE, boot = FALSE))$value
  expect_match(r$method, "Wilcoxon")
  expect_equal(r$global_pvalue, wilcox.test(y ~ g, data = d, exact = FALSE)$p.value, tolerance = 1e-3)
})

test_that("factor() et as.factor() sont acceptes dans la formule", {
  m <- mtcars
  m$am <- factor(m$am)
  ref <- bilan(mpg ~ am + wt, data = m)$value$global_pvalue
  expect_equal(bilan(mpg ~ factor(am) + wt, data = mtcars)$value$global_pvalue, ref)
  expect_equal(bilan(mpg ~ as.factor(am) + wt, data = mtcars)$value$global_pvalue, ref)
})

test_that("une variable absente est nommee", {
  expect_error(silence(m.test(y ~ zzz, data = data.frame(y = 1:10))), "zzz")
})

test_that("la voie rlm() conclut sur le modele robuste", {
  m <- mtcars
  m$am <- factor(m$am)
  b <- bilan(mpg ~ am + wt, data = m)
  expect_true(any(grepl("tests de Wald robustes", b$output)))
  expect_true(any(grepl("ROBUSTE", b$output)))
  expect_equal(p_only(mpg ~ am + wt, data = m), b$value$global_pvalue)
  expect_lt(min(b$value$p.value[upper.tri(b$value$p.value)]), 1)
})

test_that("sans methode robuste applicable, le repli classique est annonce", {
  m <- mtcars
  m$am <- factor(m$am)
  m$cyl <- factor(m$cyl)
  b <- bilan(mpg ~ cyl * am + wt, data = m)
  expect_false(is.null(b$value))
  expect_true(any(grepl("Retour vers le mod.le classique", b$output)))
})

test_that("les post-hocs ANCOVA couvrent tous les facteurs", {
  set.seed(4)
  d <- expand.grid(F = c("a", "b", "c"), G = c("u", "v"), rep = 1:12, stringsAsFactors = FALSE)
  d$x <- runif(nrow(d), 10, 20)
  d$y <- 2 * (d$F == "c") + 3 * (d$G == "v") + 0.3 * d$x + rnorm(nrow(d), 0, 0.5)
  r <- bilan(y ~ F + G + x, data = d)$value
  expect_setequal(names(r$groups_by_factor), c("F", "G"))
  expect_identical(r$groups, r$groups_by_factor$F)
})

test_that("MANOVA : les tests multivaries tournent et les assomptions sont numerotees", {
  set.seed(1)
  d <- data.frame(g = rep(c("A", "B"), each = 10), a = rnorm(20))
  d$b <- d$a + rnorm(20, 0, 0.6)
  b <- bilan(cbind(a, b) ~ g, data = d)
  expect_true(any(grepl("Henze-Zirkler p-value", b$output)))
  for (i in 1:6) expect_true(any(grepl(paste0("ASSOMPTION ", i, "/6"), b$output)))
  noname <- bilan(cbind(a, d$b) ~ g, data = d)
  expect_false(any(grepl("argument manquant", noname$output)))
})

test_that("MANOVA : la non-normalite est detectee et oriente vers la voie robuste", {
  set.seed(2)
  d <- data.frame(g = rep(c("A", "B", "C"), each = 30), a = rexp(90), b = rexp(90)^2)
  d$a[1:3] <- 40
  expect_true(any(grepl("Wilks robuste", bilan(cbind(a, b) ~ g, data = d)$output)))
})

test_that("les colonnes de $groups gardent leur nom", {
  r <- silence(m.test(iris[, 1], iris[, 5], plot = FALSE, iter = 30))$value
  expect_false(any(grepl("..", colnames(r$groups), fixed = TRUE)))
})

test_that("medpb2 : deux groupes proches partagent une lettre", {
  set.seed(6)
  d <- data.frame(x = c(rnorm(60, 0), rnorm(60, 20), rnorm(60, 0.01), rnorm(60, 40)),
                  g = rep(c("P", "Q", "R", "S"), each = 60))
  l <- KefiR:::.pairwise_medpb2(d$x, d$g, nboot = 300)$groups
  expect_identical(l$medpb2_Holm[l$categories == "P"], l$medpb2_Holm[l$categories == "R"])
  expect_length(unique(l$medpb2_Holm), 3)
})

test_that("dunnTest reste silencieux en bootstrap", {
  messages <- character(0)
  withCallingHandlers(
    utils::capture.output(pairwise(iris[, 1], iris[, 5], type = "dunn")),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("Dunn|Kruskal", messages)))
})

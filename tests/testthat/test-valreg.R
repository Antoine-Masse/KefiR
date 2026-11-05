# Le verdict de valreg() est confronte a des modeles dont le defaut est connu.

verdict <- function(reg, ...) {
  value <- NULL
  utils::capture.output(value <- suppressMessages(suppressWarnings(valreg(reg, boot = FALSE, ...))))
  value
}

make_data <- function() {
  set.seed(11)
  n <- 120
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), f8 = factor(rep(letters[1:8], 15)),
                  f3 = factor(rep(c("u", "v", "w"), 40)))
  d$y <- 2 + 1.5 * d$x1 - d$x2 + rnorm(n, 0, 0.8)
  d
}

test_that("un modele sain est valide", {
  expect_true(verdict(lm(y ~ x1 + x2, make_data())))
})

test_that("un facteur a 8 modalites n'est ni bloque ni pris pour de la colinearite", {
  d <- make_data()
  d$y8 <- d$y + as.integer(d$f8) * 0.8
  expect_true(verdict(lm(y8 ~ x1 + f8, d)))
})

test_that("une constante non significative n'invalide pas le modele", {
  d <- make_data()
  d$yc <- d$y - mean(d$y)
  expect_true(verdict(lm(yc ~ x1 + x2, d)))
})

test_that("les defauts connus invalident le modele", {
  d <- make_data()
  expect_false(verdict(lm(y ~ x1 + x2 + f3, d)))              # terme sans effet
  d$x3 <- d$x1 + rnorm(nrow(d), 0, 0.01)
  expect_false(verdict(lm(y ~ x1 + x3 + x2, d)))              # colinearite
  d$yb <- rnorm(nrow(d))
  expect_false(verdict(lm(yb ~ x1, d)))                       # modele non significatif
  d$yh <- d$x1 + rnorm(nrow(d), 0, abs(d$x1) * 2 + 0.1)
  expect_false(verdict(lm(yh ~ x1, d)))                       # heteroscedasticite
  d$yq <- d$x1^2 * 3 + rnorm(nrow(d), 0, 0.3)
  expect_false(verdict(lm(yq ~ x1, d)))                       # non-linearite
  d$ye <- d$x1 + rexp(nrow(d), 0.3)^2
  expect_false(verdict(lm(ye ~ x1, d)))                       # residus non normaux
})

test_that("une limite explicite reste respectee et explique comment la depasser", {
  d <- make_data()
  output <- utils::capture.output(v <- suppressWarnings(valreg(lm(y ~ x1 + f8, d), nvar = 5, boot = FALSE)))
  expect_false(v)
  expect_true(any(grepl("nvar = 9", output)))
})

test_that("un modele sans predicteur est refuse sans erreur", {
  expect_false(verdict(lm(y ~ 1, make_data())))
})

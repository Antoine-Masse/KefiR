# m.test ne doit pas reinitialiser le generateur aleatoire de l'utilisateur.
# WRS2::med1way() appelle set.seed(2) : sans protection, les jeux simules apres
# un appel a m.test() se repetent en boucle.

test_that("la voie med1way laisse le generateur aleatoire intact", {
  set.seed(20261115)
  gen <- function() data.frame(y = c(rnorm(8, 13, 5), rnorm(20, 10, 1), rnorm(20, 10, 1)),
                               g = rep(c("a", "b", "c"), c(8, 20, 20)))
  d <- gen(); d <- gen(); d <- gen()
  r <- suppressWarnings(suppressMessages(
    m.test(y ~ g, data = d, return = TRUE, verbose = FALSE, plot = FALSE)))
  skip_if_not(any(grepl("med1way", r$method)), "la voie med1way n'est pas prise sur ce jeu")
  set.seed(1); attendu <- runif(2)
  set.seed(1); runif(1)
  suppressWarnings(suppressMessages(
    m.test(y ~ g, data = d, return = FALSE, verbose = FALSE, plot = FALSE)))
  expect_identical(runif(1), attendu[2])
})

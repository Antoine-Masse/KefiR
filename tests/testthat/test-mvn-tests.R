# Normalite multivariee : formes closes internes, comparees au paquet MVN.

test_that("Mardia et Henze-Zirkler retrouvent les valeurs de reference sur iris", {
  x <- as.matrix(iris[, 1:4])
  ma <- KefiR:::.mvn_mardia(x)
  # References : MVN 6 sur iris (asymetrie < 0.001, aplatissement 0.818, HZ < 0.001)
  expect_lt(ma[["skew"]], 0.001)
  expect_equal(round(ma[["kurt"]], 3), 0.818)
  expect_lt(KefiR:::.mvn_hz(x), 0.001)
})

test_that("Mardia et Henze-Zirkler donnent les p-values de MVN", {
  skip_if_not_installed("MVN")
  skip_if_not("mvn_test" %in% names(formals(MVN::mvn)))
  set.seed(11)
  jeux <- list(matrix(rnorm(30), 15), matrix(rnorm(38), 19), matrix(rnorm(120), 40),
               matrix(rnorm(800), 200), matrix(rexp(150), 50), matrix(rt(160, 3), 40))
  num <- function(p) suppressWarnings(as.numeric(as.character(p)))
  for (x in jeux) {
    x <- as.data.frame(x)
    ref_ma <- num(MVN::mvn(data = x, mvn_test = "mardia")$multivariate_normality$p.value)
    ref_hz <- num(MVN::mvn(data = x, mvn_test = "hz")$multivariate_normality$p.value)
    moi_ma <- unname(KefiR:::.mvn_mardia(x)); moi_hz <- KefiR:::.mvn_hz(x)
    for (i in 1:2) if (is.na(ref_ma[i])) expect_lt(moi_ma[i], 0.001) else expect_equal(round(moi_ma[i], 3), ref_ma[i])
    if (is.na(ref_hz[1])) expect_lt(moi_hz, 0.001) else expect_equal(round(moi_hz, 3), ref_hz[1])
  }
})

test_that("les cas incalculables rendent NULL au lieu d'une erreur", {
  expect_null(KefiR:::.mvn_mardia(matrix(rnorm(6), 3)))       # trop peu de lignes
  z <- rnorm(30)
  expect_null(KefiR:::.mvn_hz(cbind(z, 2 * z)))               # covariance singuliere
  x <- matrix(rnorm(60), 30); x[1, 1] <- NA
  expect_length(KefiR:::.mvn_mardia(x), 2)                    # manquants ecartes
})

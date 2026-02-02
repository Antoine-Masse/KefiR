test_that("KefiR package loads successfully", {
  # Test que le package se charge sans erreur
  expect_true(require(KefiR, quietly = TRUE))
})

test_that("main function m.test exists", {
  # Test que la fonction principale existe
  expect_true(exists("m.test"))
})

# Audit de m.test() (novembre 2026) : un test par defaut corrige.
# Chaque test echoue avant la correction et passe apres ; voir dev/mtest_audit/RAPPORT.md.

audit_quiet <- function(expr) {
  value <- NULL
  output <- utils::capture.output(value <- suppressMessages(suppressWarnings(expr)))
  list(value = value, output = output)
}

make_split_plot <- function() {
  set.seed(11)
  n <- 12
  ids <- factor(1:n)
  G <- factor(rep(c("A", "B"), each = n / 2))
  d <- expand.grid(id = ids, T = factor(c("t1", "t2", "t3")))
  d$G <- G[as.integer(d$id)]
  u <- rnorm(n, 0, 1)
  d$y <- 10 + u[as.integer(d$id)] + ifelse(d$G == "B" & d$T == "t3", 1.2, 0) +
    rnorm(nrow(d), 0, 0.8)
  d
}

test_that("b. modele mixte force (plan desequilibre) avec verbose = FALSE : pas de plantage", {
  d <- make_split_plot()
  d <- d[-c(1, 20), ]
  expect_no_error(audit_quiet(m.test(y ~ G * T + Error(id / T), data = d,
                                     plot = FALSE, verbose = FALSE)))
})

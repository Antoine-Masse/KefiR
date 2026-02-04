# Tests for wilcoxon.cut.test() and evolreg()

set.seed(1)

# wilcoxon.cut.test()
X <- data.frame(
  y = rnorm(40),
  g = rnorm(40)
)

res_w <- wilcoxon.cut.test(
  x = "y",
  group = "g",
  data = X,
  prop = 0.1,
  verbose = FALSE,
  bootstrap = 0
)

print(res_w)

# evolreg() - keep iterations low for a quick smoke test
res_e <- evolreg(
  data = mtcars,
  Y = "mpg",
  iter = 20,
  nbind = 5,
  plot = FALSE,
  verbose = FALSE,
  fast = TRUE
)

print(is.list(res_e))

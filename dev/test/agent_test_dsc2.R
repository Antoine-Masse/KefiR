# Test minimal for dsc2()

set.seed(1)

data(mtcars)

reg <- lm(mpg ~ wt, data = mtcars)

# Minimal dsc row compatible with the model terms
my_dsc <- data.frame(wt = mean(mtcars$wt))

out <- dsc2(
  data = mtcars,
  reg = reg,
  dsc = my_dsc,
  iter = 20,
  plot = FALSE,
  return = TRUE
)

print(head(out))

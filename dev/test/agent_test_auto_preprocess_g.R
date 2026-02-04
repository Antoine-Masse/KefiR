# Test minimal for auto_preprocess_g()

set.seed(1)

x <- rnorm(20)

g <- data.frame(
  bin = sample(c(0, 1), 20, replace = TRUE),
  num1 = rnorm(20),
  num2 = rnorm(20),
  fac = sample(letters[1:3], 20, replace = TRUE)
)

res <- auto_preprocess_g(
  x = x,
  g = g,
  gamma = 3,
  start_bins = 5,
  cor_threshold = 0.8,
  debug = TRUE,
  verbose = TRUE
)

print(str(res))

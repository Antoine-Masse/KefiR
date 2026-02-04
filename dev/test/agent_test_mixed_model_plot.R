# Test minimal for .mixed_model_analysis() and .plot_with_letters()
# NOTE: Requires {lmerTest} and {lme4} installed.

set.seed(1)

df <- data.frame(
  y = rnorm(30),
  id = factor(rep(1:10, each = 3)),
  time = factor(rep(1:3, times = 10)),
  group = factor(rep(c("A", "B"), each = 15))
)

# Mixed model with random intercept by subject
mixed_res <- .mixed_model_analysis(
  formula = y ~ time * group + (1 | id),
  data = df,
  id = "id",
  within = "time",
  between = "group",
  verbose = FALSE,
  code = FALSE,
  debug = FALSE
)

print(is.list(mixed_res))

# Minimal posthoc_result structure to test plotting with letters
posthoc <- list(
  groups = data.frame(
    categories = levels(df$group),
    letters = c("a", "b")
  )
)

.plot_with_letters(
  x = df$y,
  g = df$group,
  posthoc_result = posthoc,
  main = "Test plot with letters",
  ylab = "y",
  xlab = "group",
  verbose = FALSE,
  debug = FALSE
)

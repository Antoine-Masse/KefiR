.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths())); suppressMessages(library(BayesFactor))
set.seed(2026)
for (i in 1:20000) { n <- sample(c(10, 20), 1); x <- sample(1:5, n, TRUE, c(.1, .2, .4, .2, .1)); y <- sample(1:5, n, TRUE, c(.1, .2, .4, .2, .1))
	r <- tryCatch(correlationBF(y, x, rscale = 1/3, nullInterval = c(0, 1)), error = function(e) e)
	if (inherits(r, "error")) { cat("fail at", i, "n =", n, "r =", cor(x, y), "\n"); dput(list(x = x, y = y)); print(conditionMessage(r)); break } }

# Rate of "point away from the cloud" on normal data, per data set, for several robust distances.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
cat("robustbase:", requireNamespace("robustbase", quietly = TRUE), "\n")
flag <- function(x, y, how) {
	X <- cbind(x, y); p <- 2
	cr <- withr::with_preserve_seed({ set.seed(1); MASS::cov.rob(X, method = "mcd") })
	d0 <- mahalanobis(X, cr$center, cr$cov)
	if (how == "raw") return(any(d0 > qchisq(0.999, p)))
	w <- d0 <= qchisq(0.975, p)
	m <- colMeans(X[w, , drop = FALSE]); S <- cov(X[w, , drop = FALSE]) * 0.975 / pchisq(qchisq(0.975, p), p + 2)
	d1 <- mahalanobis(X, m, S)
	if (how == "rw999") return(any(d1 > qchisq(0.999, p)))
	if (how == "rw_bonf") return(any(d1 > qchisq(1 - 0.01 / nrow(X), p)))
	if (how == "rob_bonf") return(any(robustbase::covMcd(X)$mah > qchisq(1 - 0.01 / nrow(X), p)))
}
set.seed(5)
for (n in c(10, 20, 30, 50, 100, 300)) {
	r <- sapply(c("raw", "rw999", "rw_bonf"), function(h) mean(replicate(300, { x <- rnorm(n); y <- 0.4 * x + rnorm(n); flag(x, y, h) })))
	cat(n, round(100 * r, 1), "\n")
}
# Power to catch one planted point at (3, -3) in standard units (against a link of 0.6)
for (n in c(10, 20, 50, 100)) {
	r <- sapply(c("raw", "rw999", "rw_bonf"), function(h) mean(replicate(300, { x <- rnorm(n); y <- 0.6 * x + 0.8 * rnorm(n); x[n] <- 3; y[n] <- -3; flag(x, y, h) })))
	cat("planted", n, round(100 * r, 1), "\n")
}

.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
flag <- function(x, y, how) {
	X <- cbind(x, y); n <- nrow(X)
	if (how == "covMcd_bonf") return(any(withr::with_preserve_seed({ set.seed(1); robustbase::covMcd(X)$mah }) > qchisq(1 - 0.01 / n, 2)))
	if (how == "covMcd_999") return(any(withr::with_preserve_seed({ set.seed(1); robustbase::covMcd(X)$mah }) > qchisq(0.999, 2)))
	if (how == "cook1") { cd <- cooks.distance(lm(y ~ x)); cd2 <- cooks.distance(lm(x ~ y)); return(any(cd > 1 | cd2 > 1)) }
	if (how == "cook4n") { cd <- cooks.distance(lm(y ~ x)); return(any(cd > 4 / n)) }
}
H <- c("covMcd_bonf", "covMcd_999", "cook1")
set.seed(5)
for (n in c(10, 20, 30, 50, 100, 300)) {
	r <- sapply(H, function(h) mean(replicate(300, { x <- rnorm(n); y <- 0.4 * x + rnorm(n); flag(x, y, h) })))
	cat(n, round(100 * r, 1), "\n")
}
for (n in c(10, 20, 50, 100)) {
	r <- sapply(H, function(h) mean(replicate(300, { x <- rnorm(n); y <- 0.6 * x + 0.8 * rnorm(n); x[n] <- 3; y[n] <- -3; flag(x, y, h) })))
	cat("planted", n, round(100 * r, 1), "\n")
}
for (n in c(20, 50, 100)) {
	r <- sapply(H, function(h) mean(replicate(300, { x <- rnorm(n); y <- rnorm(n); x[n] <- 6; y[n] <- 6; flag(x, y, h) })))
	cat("planted66", n, round(100 * r, 1), "\n")
}

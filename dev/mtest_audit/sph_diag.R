.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
x <- readRDS("sph_after.rds"); x <- x[x$scen == "ar1", ]
gen <- function(i) { set.seed(5000 + i); n <- 12; k <- 6
  sig <- 0.9^abs(outer(1:k, 1:k, "-")); matrix(rnorm(n * k), n, k) %*% chol(sig) }
ref <- t(sapply(x$i, function(i) {
  y <- gen(i)
  fit <- lm(y ~ 1)
  sm <- summary(car::Anova(fit, idata = data.frame(T = factor(paste0("t", 1:6))), idesign = ~T, type = 3),
                multivariate = FALSE)
  adj <- sm$pval.adjustments; uni <- sm$univariate.tests
  gg <- adj["T", "GG eps"]; hf <- min(1, adj["T", "HF eps"])
  f <- uni["T", "F value"]; d1 <- uni["T", "num Df"]; d2 <- uni["T", "den Df"]
  c(gg = gg, p_gg = adj["T", "Pr(>F[GG])"],
    p_rule = if (gg <= 0.75) adj["T", "Pr(>F[GG])"] else pf(f, hf * d1, hf * d2, lower.tail = FALSE))
}))
cat("GG eps moyen", round(mean(ref[, "gg"]), 3), " part GG<=0.75", mean(ref[, "gg"] <= 0.75), "\n")
cat("Rejets : m.test", mean(x$p_mtest < .05), " regle GG/HF directe", mean(ref[, "p_rule"] < .05),
    " GG seul", mean(ref[, "p_gg"] < .05), "\n")
d <- abs(x$p_mtest - ref[, "p_rule"])
cat("ecart max m.test vs regle", max(d), " jeux avec ecart > 1e-6 :", sum(d > 1e-6), "\n")
if (any(d > 1e-6)) print(head(cbind(i = x$i, mtest = x$p_mtest, rule = ref[, "p_rule"], aov = x$p_aov)[d > 1e-6, ], 15))

## Calcul de reference independant : ANCOVA bayesienne (agent_57)
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressPackageStartupMessages({library(BayesFactor); library(bayestestR); library(car); library(emmeans)})
options(width = 120)
ITER_BF <- 50000; ITER_POST <- 20000
med_hdi <- function(v) { h <- bayestestR::hdi(v, ci = 0.95); c(med = median(v), lo = h$CI_low, hi = h$CI_high) }
zones <- function(grid, ok) { # intervalles contigus ou ok est vrai
  if (!any(ok)) return("aucune")
  r <- rle(ok); fin <- cumsum(r$lengths); deb <- fin - r$lengths + 1
  paste(sprintf("[%.3f ; %.3f] (%d pts)", grid[deb[r$values]], grid[fin[r$values]], r$lengths[r$values]), collapse = " U ")
}
cat(R.version.string, "| BayesFactor", as.character(packageVersion("BayesFactor")),
    "| bayestestR", as.character(packageVersion("bayestestR")), "| car", as.character(packageVersion("car")),
    "| emmeans", as.character(packageVersion("emmeans")), "\n")

## ================= Jeu G =================
cat("\n================= JEU G =================\n")
set.seed(31)
d <- data.frame(engrais = factor(rep(c("E1","E2","E3"), each = 15)), surface = runif(45, 5, 15))
d$rendement <- 20 + 1.5 * d$surface + c(0, 0, 3)[d$engrais] + rnorm(45, 0, 2)
d$xc <- d$surface - mean(d$surface)
set.seed(1001)
bf_gx <- lmBF(rendement ~ engrais + xc, data = d, rscaleFixed = 0.5, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
bf_x  <- lmBF(rendement ~ xc, data = d, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
bf <- bf_gx / bf_x
lbf <- extractBF(bf, logbf = TRUE)
cat(sprintf("log BF(engrais) = %.4f (BF = %.4g, erreur rel. = %.3g)\n", lbf$bf, exp(lbf$bf), extractBF(bf)$error))
cat(sprintf("P(engrais), a priori 1/2 = %.6f\n", plogis(lbf$bf)))
set.seed(1002)
po <- posterior(bf_gx, iterations = ITER_POST)
cat("colonnes posterior :", paste(colnames(po), collapse = ", "), "\n")
e1 <- po[, "engrais-E1"]; e2 <- po[, "engrais-E2"]; e3 <- po[, "engrais-E3"]
res <- rbind("E2-E1" = med_hdi(e2 - e1), "E3-E1" = med_hdi(e3 - e1), "E3-E2" = med_hdi(e3 - e2))
cat("Ecarts ajustes bayesiens (mediane, HDI 95 %) :\n"); print(round(res, 4))
fit <- lm(rendement ~ engrais + xc, data = d)
cat("car::Anova type II :\n"); print(Anova(fit, type = 2))
em <- emmeans(fit, ~ engrais)
ct <- pairs(em, reverse = TRUE, adjust = "none")
cat("emmeans (paires, adjust = none, reverse) :\n")
print(cbind(summary(ct)[, c("contrast","estimate","SE","p.value")], confint(ct)[, c("lower.CL","upper.CL")]), digits = 5)

## ================= Jeu H =================
cat("\n================= JEU H =================\n")
set.seed(32)
d <- data.frame(methode = factor(rep(c("classique","projet"), each = 25)), note = runif(50, 6, 18))
d$score <- 5 + 0.6 * d$note + (d$methode == "projet") * (0.5 * (d$note - 12)) + rnorm(50, 0, 1.5)
mx <- mean(d$note); d$xc <- d$note - mx
set.seed(2001)
b_x   <- lmBF(score ~ xc, data = d, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
b_gx  <- lmBF(score ~ methode + xc, data = d, rscaleFixed = 0.5, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
b_ful <- lmBF(score ~ methode + xc + methode:xc, data = d, rscaleFixed = 0.5, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
l_gx <- extractBF(b_gx / b_x, logbf = TRUE)$bf; l_ful <- extractBF(b_ful / b_x, logbf = TRUE)$bf
cat(sprintf("log BF[g+x vs x] = %.4f ; log BF[g+x+g:x vs x] = %.4f\n", l_gx, l_ful))
lw <- c(x = log(1/2), gx = log(1/4) + l_gx, full = log(1/4) + l_ful)
pp <- exp(lw - max(lw)); pp <- pp / sum(pp)
cat(sprintf("Posterieurs : x %.6f | g+x %.6f | g+x+g:x %.6f\n", pp[1], pp[2], pp[3]))
cat(sprintf("P(methode) = %.6f ; P(methode:note) = %.6f\n", pp[2] + pp[3], pp[3]))
set.seed(2002)
po <- posterior(b_ful, iterations = ITER_POST)
cn <- colnames(po)
cat("colonnes posterior :", paste(cn, collapse = ", "), "\n")
ic <- grep(":", cn, value = TRUE)
i_cl <- ic[grepl("classique", ic)]; i_pr <- ic[grepl("projet", ic)]
cat("colonnes d'interaction utilisees :", i_cl, "/", i_pr, "\n")
dint <- po[, "methode-projet"] - po[, "methode-classique"]
dslp <- po[, i_pr] - po[, i_cl]
q <- quantile(d$note, c(.25, .5, .75))
cat("Quartiles de note :", paste(sprintf("%.4f", q), collapse = " / "), " ; moyenne :", round(mx, 4), "\n")
resq <- t(sapply(q, function(x0) med_hdi(dint + dslp * (x0 - mx)))); rownames(resq) <- names(q)
cat("Ecart projet - classique aux quartiles (bayesien, mediane, HDI 95 %) :\n"); print(round(resq, 4))
grid <- seq(min(d$note), max(d$note), length.out = 101)
hh <- t(sapply(grid, function(x0) { h <- bayestestR::hdi(dint + dslp * (x0 - mx), ci = .95); c(h$CI_low, h$CI_high) }))
okb <- hh[, 1] > 0 | hh[, 2] < 0
cat("Etendue de note : [", round(min(d$note), 4), ";", round(max(d$note), 4), "] ; pas de grille", round(diff(grid[1:2]), 4), "\n")
cat("Zone bayesienne (HDI exclut 0) :", zones(grid, okb), "\n")
cat("  dont HDI > 0 :", zones(grid, hh[,1] > 0), "| HDI < 0 :", zones(grid, hh[,2] < 0), "\n")
fit <- lm(score ~ methode * xc, data = d)
cat("car::Anova type II :\n"); print(Anova(fit, type = 2))
fit2 <- lm(score ~ methode * note, data = d)
emq <- emmeans(fit2, ~ methode | note, at = list(note = as.numeric(q)))
ctq <- contrast(emq, "revpairwise", adjust = "none")
cat("emmeans ecart projet - classique aux quartiles :\n")
print(cbind(summary(ctq)[, c("note","estimate","SE","p.value")], confint(ctq)[, c("lower.CL","upper.CL")]), digits = 5)
b <- coef(fit); V <- vcov(fit); tc <- qt(.975, df.residual(fit))
a0 <- b[["methodeprojet"]]; a1 <- b[["methodeprojet:xc"]]
v00 <- V["methodeprojet","methodeprojet"]; v11 <- V["methodeprojet:xc","methodeprojet:xc"]; v01 <- V["methodeprojet","methodeprojet:xc"]
est <- a0 + a1 * (grid - mx); se <- sqrt(v00 + (grid - mx)^2 * v11 + 2 * (grid - mx) * v01)
okf <- abs(est / se) > tc
cat("Johnson-Neyman (grille 101 pts) :", zones(grid, okf), "\n")
cat("  dont ecart > 0 :", zones(grid, okf & est > 0), "| ecart < 0 :", zones(grid, okf & est < 0), "\n")
A <- a1^2 - tc^2 * v11; B <- 2 * (a0 * a1 - tc^2 * v01); C <- a0^2 - tc^2 * v00
disc <- B^2 - 4 * A * C
if (disc >= 0) { z <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
  cat(sprintf("Bornes JN exactes (note) : %.4f et %.4f ; A = %.4g (A > 0 : significatif hors de l'intervalle)\n", z[1] + mx, z[2] + mx, A))
} else cat("Pas de borne JN reelle\n")

## ================= Jeu I =================
cat("\n================= JEU I =================\n")
set.seed(33)
d <- data.frame(regime = factor(rep(c("R1","R2"), each = 15)), age = c(runif(15, 20, 35), runif(15, 32, 50)))
d$poids <- 60 + 0.4 * d$age + 2 * (d$regime == "R2") + rnorm(30, 0, 3)
d$xc <- d$age - mean(d$age)
set.seed(3001)
bf_gx <- lmBF(poids ~ regime + xc, data = d, rscaleFixed = 0.5, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
bf_x  <- lmBF(poids ~ xc, data = d, rscaleCont = sqrt(2)/4, iterations = ITER_BF)
lbf <- extractBF(bf_gx / bf_x, logbf = TRUE)$bf
cat(sprintf("log BF(regime) = %.4f (BF = %.4g) ; P(regime), a priori 1/2 = %.6f\n", lbf, exp(lbf), plogis(lbf)))
set.seed(3002)
po <- posterior(bf_gx, iterations = ITER_POST)
cat("Ecart ajuste R2 - R1 (bayesien) :\n"); print(round(med_hdi(po[, "regime-R2"] - po[, "regime-R1"]), 4))
fit <- lm(poids ~ regime + xc, data = d)
cat("car::Anova type II :\n"); print(Anova(fit, type = 2))
em <- emmeans(fit, ~ regime); ct <- pairs(em, reverse = TRUE, adjust = "none")
cat("emmeans R2 - R1 :\n")
print(cbind(summary(ct)[, c("contrast","estimate","SE","p.value")], confint(ct)[, c("lower.CL","upper.CL")]), digits = 5)
rg <- tapply(d$age, d$regime, range)
lo <- max(sapply(rg, `[`, 1)); hi <- min(sapply(rg, `[`, 2))
cat(sprintf("Etendues : R1 [%.4f ; %.4f], R2 [%.4f ; %.4f] ; commune [%.4f ; %.4f]\n", rg$R1[1], rg$R1[2], rg$R2[1], rg$R2[2], lo, hi))
inside <- d$age >= lo & d$age <= hi
for (g in levels(d$regime)) cat(sprintf("  %s : %d / %d = %.4f dans l'etendue commune\n", g, sum(inside[d$regime == g]), sum(d$regime == g), mean(inside[d$regime == g])))

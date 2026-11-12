## agent_55 : calcul de reference independant (JZS, Rouder 2009 / 2012)
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressPackageStartupMessages(library(BayesFactor))
cat("BayesFactor", as.character(packageVersion("BayesFactor")), "\n\n")

## ---------- A. t apparie ----------
set.seed(21); n <- 14; patient <- 1:14
s <- rnorm(14, 0, 12); e <- rnorm(28, 0, 5)
avant <- 140 + s + e[1:14]; apres <- 140 + s - 6 + e[15:28]
avant <- avant[-3]; apres <- apres[-3]
d <- apres - avant; N <- length(d)
tt <- mean(d) / (sd(d) / sqrt(N)); nu <- N - 1
jzs_bf <- function(t, n, r) {
  nu <- n - 1
  # g ~ IG(1/2, r^2/2) ; integrande sur g
  f <- function(g) exp(-0.5 * log1p(n * g) - (nu + 1) / 2 * log1p(t^2 / ((1 + n * g) * nu)) +
                       log(r) - 0.5 * log(2 * pi) - 1.5 * log(g) - r^2 / (2 * g) +
                       (nu + 1) / 2 * log1p(t^2 / nu))
  integrate(f, 0, Inf, rel.tol = 1e-10)$value
}
bfA <- jzs_bf(tt, N, 0.5)
bfA_pkg <- ttestBF(x = apres, y = avant, paired = TRUE, rscale = 0.5)
cat("=== A ===\nn paires =", N, " moyenne(apres-avant) =", mean(d), " dz =", mean(d) / sd(d), " t =", tt, "\n")
cat("log BF10 (formule) =", log(bfA), "\nlog BF10 (ttestBF) =", extractBF(bfA_pkg, logbf = TRUE)$bf,
    "\necart =", log(bfA) - extractBF(bfA_pkg, logbf = TRUE)$bf, "\n\n")

## ---------- B ----------
set.seed(22)
sujet <- factor(rep(1:12, 4)); session <- factor(rep(paste0("S", 1:4), each = 12))
rt <- 450 + rnorm(12, 0, 40)[sujet] + c(0, -5, -25, -30)[session] + rnorm(48, 0, 15)
dB <- data.frame(rt, sujet, session)

bfB_pkg <- anovaBF(rt ~ session + sujet, data = dB, whichRandom = "sujet", rscaleFixed = 0.5,
                   iterations = 50000, progress = FALSE)
cat("=== B1 ===\nanovaBF :\n"); print(bfB_pkg)
lbfB_pkg <- extractBF(bfB_pkg, logbf = TRUE)

## propre integrale 2D (Rouder 2012) : y = mu + X_f a + X_s b + e,
## a ~ N(0, g_f s2) sur base orthonormee des contrastes somme nulle, b ~ N(0, g_s s2) indicatrices,
## mu plat, p(s2) ~ 1/s2, g ~ IG(1/2, r^2/2). Integrale analytique de mu et s2 :
## p(y|g) prop |Sigma|^{-1/2} (1'Sigma^-1 1)^{-1/2} (y'Q y)^{-(N-1)/2}
Nb <- length(rt)
Xs <- model.matrix(~ sujet - 1)
a <- nlevels(session); C <- diag(a) - 1 / a
Qb <- eigen(C, symmetric = TRUE)$vectors[, 1:(a - 1)]   # base orthonormee des contrastes
Xf <- model.matrix(~ session - 1) %*% Qb
logml <- function(X, gvec, y) {
  Sig <- diag(length(y)) + X %*% (gvec * t(X))
  R <- chol(Sig); ld <- 2 * sum(log(diag(R)))
  one <- rep(1, length(y))
  Si1 <- backsolve(R, forwardsolve(t(R), one)); Siy <- backsolve(R, forwardsolve(t(R), y))
  s11 <- sum(Si1); q <- sum(y * Siy) - sum(one * Siy)^2 / s11
  -0.5 * ld - 0.5 * log(s11) - (length(y) - 1) / 2 * log(q)
}
ldig <- function(g, r) log(r) - 0.5 * log(2 * pi) - 1.5 * log(g) - r^2 / (2 * g)
run_B1 <- function(Xs_use, lab) {
  # integration sur u = log g ; densite de u = g * pi(g)
  ref <- logml(Xs_use, rep(1, ncol(Xs_use)), rt)
  f0 <- function(u) sapply(u, function(ui) { g <- exp(ui)
    exp(logml(Xs_use, rep(g, ncol(Xs_use)), rt) - ref + ldig(g, 1) + ui) })
  m0 <- integrate(f0, -30, 30, rel.tol = 1e-8)$value
  X1 <- cbind(Xf, Xs_use)
  inner <- function(uf) sapply(uf, function(ufi) {
    h <- function(us) sapply(us, function(usi) {
      gv <- c(rep(exp(ufi), ncol(Xf)), rep(exp(usi), ncol(Xs_use)))
      exp(logml(X1, gv, rt) - ref + ldig(exp(ufi), 0.5) + ufi + ldig(exp(usi), 1) + usi) })
    integrate(h, -30, 30, rel.tol = 1e-8)$value })
  m1 <- integrate(inner, -30, 30, rel.tol = 1e-7)$value
  cat(lab, ": log BF (session+sujet vs sujet) =", log(m1 / m0), "\n")
  log(m1 / m0)
}
lb_unc <- run_B1(Xs, "Integrale 2D, indicatrices sujet brutes ")
lb_cen <- run_B1(scale(Xs, scale = FALSE), "Integrale 2D, indicatrices sujet centrees")
cat("anovaBF log BF =", lbfB_pkg$bf, " (+/- erreur rel.", lbfB_pkg$error, ")\n\n")

## ---------- B2 : partitions ----------
parts <- list(); lv <- paste0("S", 1:4)
rec <- function(i, lab) { if (i > 4) { parts[[length(parts) + 1]] <<- lab; return() }
  for (k in 1:(max(c(0, lab)) + 1)) rec(i + 1, c(lab, k)) }
rec(1, integer(0))
cat("=== B2 === nombre de partitions :", length(parts), "\n")
set.seed(2201)
res <- do.call(rbind, lapply(parts, function(p) {
  nom <- paste(sapply(split(lv, p), paste, collapse = "="), collapse = " | ")
  sizes <- table(p); w <- prod(factorial(sizes - 1))
  if (length(sizes) == 1) lbf <- 0 else {
    dB$bloc <- factor(p[as.integer(dB$session)])
    lbf <- extractBF(lmBF(rt ~ bloc + sujet, data = dB, whichRandom = "sujet", rscaleFixed = 0.5,
                          iterations = 50000, progress = FALSE), logbf = TRUE)$bf }
  data.frame(partition = nom, k = length(sizes), w = w, logBF = lbf, stringsAsFactors = FALSE) }))
res$prior <- ifelse(res$k == 1, 0.5, 0.5 * res$w / sum(res$w[res$k > 1]))
lp <- log(res$prior) + res$logBF; res$post <- exp(lp - max(lp)); res$post <- res$post / sum(res$post)
res <- res[order(-res$post), ]
print(res, row.names = FALSE, digits = 4)
cat("\nTrois meilleures :\n"); print(head(res[, c("partition", "post")], 3), row.names = FALSE)
cat("\nP(les deux sessions different) :\n")
labs <- lapply(parts, identity); postv <- setNames(res$post, res$partition)
for (i in 1:3) for (j in (i + 1):4) {
  pr <- sum(sapply(seq_along(parts), function(q) {
    p <- parts[[q]]; nom <- paste(sapply(split(lv, p), paste, collapse = "="), collapse = " | ")
    if (p[i] != p[j]) postv[[nom]] else 0 }))
  cat(sprintf("  %s vs %s : %.4f\n", lv[i], lv[j], pr)) }

## ---------- C ----------
set.seed(23)
dC <- expand.grid(id = factor(1:20), seance = factor(c("s1", "s2", "s3")))
dC$groupe <- factor(ifelse(as.integer(dC$id) <= 10, "temoin", "entraine"))
dC$score <- 50 + rnorm(20, 0, 6)[dC$id] + 4 * (dC$seance != "s1") * (dC$groupe == "entraine") + rnorm(nrow(dC), 0, 3)
set.seed(2301)
bfC <- anovaBF(score ~ groupe * seance + id, data = dC, whichRandom = "id", rscaleFixed = 0.5,
               iterations = 50000, progress = FALSE)
cat("\n=== C ===\n"); print(bfC)
print(extractBF(bfC, logbf = TRUE)[, 1:2])

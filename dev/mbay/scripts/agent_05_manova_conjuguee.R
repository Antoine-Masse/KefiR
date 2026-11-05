# 2f. Vraisemblance marginale exacte, regression multivariee, a priori conjugue MN-IW propre.
# Sigma ~ IW(nu0, S0) ; B | Sigma ~ MN(B0, Lambda0^-1, Sigma)
# Formule : log m(Y) = -(n p / 2) log(pi) + (p/2) (log|L0| - log|Ln|)
#                     + (nu0/2) log|S0| - (nun/2) log|Sn| + lmvgamma(nun/2) - lmvgamma(nu0/2)
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
lmvgamma <- function(a, p) p*(p-1)/4*log(pi) + sum(lgamma(a + (1 - seq_len(p))/2))
ldet <- function(M) as.numeric(determinant(M, logarithm=TRUE)$modulus)

logm <- function(Y, X, B0, L0, nu0, S0) {
  Y <- as.matrix(Y); n <- nrow(Y); p <- ncol(Y)
  Ln <- crossprod(X) + L0
  Bn <- solve(Ln, crossprod(X, Y) + L0 %*% B0)
  nun <- nu0 + n
  E <- Y - X %*% Bn; D <- Bn - B0
  Sn <- S0 + crossprod(E) + t(D) %*% L0 %*% D
  -(n*p/2)*log(pi) + (p/2)*(ldet(L0) - ldet(Ln)) + (nu0/2)*ldet(S0) - (nun/2)*ldet(Sn) +
    lmvgamma(nun/2, p) - lmvgamma(nu0/2, p)
}

# contrastes orthonormes : a priori echangeable entre niveaux
Qstar <- function(a) { e <- eigen(diag(a) - 1/a); e$vectors[, seq_len(a-1), drop=FALSE] }
design <- function(f) { f <- factor(f); a <- nlevels(f); cbind(1, model.matrix(~0+f) %*% Qstar(a)) }

# BF10 facteur vs intercept seul. Y est centre ; intercept : a priori faible commun (lam_int)
bf_manova <- function(Y, f, lam=1, nu0=NULL, S0="diagvar", lam_int=1e-3) {
  Y <- as.matrix(Y); p <- ncol(Y); n <- nrow(Y)
  if (is.null(nu0)) nu0 <- p + 2
  if (is.character(S0)) {
    S0 <- switch(S0,
      identity = diag(p),
      diagvar  = diag(diag(cov(Y)), p) * (nu0 - p - 1),   # E[Sigma] = diag(var(Y)), depend des donnees
      fullcov  = cov(Y) * (nu0 - p - 1))                  # E[Sigma] = cov(Y), equivariant affine
  }
  Yc <- scale(Y, center=TRUE, scale=FALSE)
  X1 <- design(f); X0 <- X1[, 1, drop=FALSE]; k <- ncol(X1)
  L1 <- diag(c(lam_int, rep(lam, k-1)), k); L0 <- matrix(lam_int)
  lm1 <- logm(Yc, X1, matrix(0, k, p), L1, nu0, S0)
  lm0 <- logm(Yc, X0, matrix(0, 1, p), L0, nu0, S0)
  lm1 - lm0
}

rmvn <- function(n, Sig) matrix(rnorm(n*ncol(Sig)), n) %*% chol(Sig)
fmt <- function(x) formatC(x, digits=6, format="g")

cat("=== A. Verification du code p=1 contre calculs independants ===\n")
set.seed(101)
n <- 20; X <- cbind(1, rnorm(n)); y <- 1 + 0.7*X[,2] + rnorm(n, sd=1.3)
B0 <- matrix(c(0.5,0),2); L0 <- diag(c(0.5,2)); nu0 <- 5; S0 <- matrix(3)
a <- logm(matrix(y), X, B0, L0, nu0, S0)
# (i) integration numerique sur sigma2 : y | s2 ~ N(X b0, s2 V), s2 ~ IG(nu0/2, S0/2)
V <- diag(n) + X %*% solve(L0) %*% t(X)
r <- y - X %*% B0; Q <- drop(t(r) %*% solve(V, r)); ldV <- ldet(V)
integrand <- function(s2) exp(-n/2*log(2*pi*s2) - ldV/2 - Q/(2*s2) +
   (nu0/2)*log(S0[1]/2) - lgamma(nu0/2) - (nu0/2+1)*log(s2) - S0[1]/(2*s2) - a)
I <- integrate(Vectorize(integrand), 0, Inf, rel.tol=1e-10)
cat("  log m formule         =", format(a, digits=12), "\n")
cat("  log m integrate(s2)   =", format(a + log(I$value), digits=12), " (ecart ", format(log(I$value), digits=3), ")\n")
# (ii) densite t multivariee fermee
nun <- nu0; d <- lgamma((nun+n)/2) - lgamma(nun/2) - n/2*log(nun*pi) - 0.5*ldet(V*S0[1]/nun) -
  (nun+n)/2*log(1 + Q/S0[1])
cat("  log m t multivariee   =", format(d, digits=12), " (ecart ", format(a-d, digits=3), ")\n")

cat("\n=== B. Verification p=2 par Monte-Carlo direct (Sigma ~ IW, 2e5 tirages) ===\n")
set.seed(102)
n <- 12; f <- gl(3,4); X <- design(f); p <- 2; Y <- rmvn(n, matrix(c(1,.4,.4,1),2))
B0 <- matrix(0,3,2); L0 <- diag(c(.1,1,1)); nu0 <- 6; S0 <- diag(c(2,3))
a <- logm(Y, X, B0, L0, nu0, S0)
V <- diag(n) + X %*% solve(L0) %*% t(X); Vi <- solve(V); ldV <- ldet(V)
M <- 2e5; W <- rWishart(M, nu0, solve(S0))
ll <- numeric(M); YtViY <- t(Y) %*% Vi %*% Y   # B0=0
for (i in seq_len(M)) {
  Prec <- W[,,i]   # Sigma^-1 ~ Wishart(nu0, S0^-1)
  ll[i] <- -(n*p/2)*log(2*pi) - (p/2)*ldV + (n/2)*ldet(Prec) - 0.5*sum(Prec*YtViY)
}
mx <- max(ll); est <- mx + log(mean(exp(ll-mx)))
w <- exp(ll-mx); se <- sd(w)/sqrt(M)/mean(w)
cat("  log m formule =", format(a, digits=10), "  Monte-Carlo =", format(est, digits=10), " +/-", format(se, digits=3), "\n")
cat("  Rq : Rscript R CMD -> ", M, "tirages\n")

cat("\n=== C. p=1 : BF conjugue vs test F, sens de l'effet ===\n")
set.seed(103)
for (eff in c(0, 0.8)) {
  f <- gl(3,15); y <- rnorm(45) + eff*(as.integer(f)-2)
  cat(sprintf("  eff=%.1f logBF10(lam=1)=%.3f  p aov=%.4g\n", eff, bf_manova(cbind(y), f), summary(aov(y~f))[[1]][1,5]))
}

cat("\n=== D. Donnees simulees : p=3, 3 groupes, n=20/groupe ; H0 et effet connu ===\n")
set.seed(104)
Sig <- matrix(.5,3,3); diag(Sig) <- 1
f <- gl(3,20)
Y0 <- rmvn(60, Sig)
mu <- rbind(c(-.6,-.3,0), c(0,0,0), c(.6,.3,0))
Y1 <- rmvn(60, Sig) + mu[as.integer(f),]
pil <- function(Y) { s <- summary(manova(Y~f), test="Pillai")$stats; c(Pillai=s[1,2], p=s[1,6]) }
cat("  H0 : logBF10 =", fmt(bf_manova(Y0,f)), " ; manova:", fmt(pil(Y0)), "\n")
cat("  H1 : logBF10 =", fmt(bf_manova(Y1,f)), " ; manova:", fmt(pil(Y1)), "\n")

cat("\n=== E. Sensibilite a lam, nu0, S0 (logBF10) ===\n")
g <- expand.grid(lam=c(.1,1/3,1,10), nu0=c(5,8,20), S0=c("identity","diagvar","fullcov"), stringsAsFactors=FALSE)
g$H0 <- mapply(function(l,n,s) bf_manova(Y0,f,lam=l,nu0=n,S0=s), g$lam,g$nu0,g$S0)
g$H1 <- mapply(function(l,n,s) bf_manova(Y1,f,lam=l,nu0=n,S0=s), g$lam,g$nu0,g$S0)
print(format(g, digits=4), row.names=FALSE)

cat("\n=== F. Invariance aux unites : Y[,1] en mg, Y[,2] en kg (Y[,3] inchangee) ===\n")
A <- diag(c(1e6, 1, 1e-3))  # on multiplie : 1 -> 1e6, 3 -> 1e-3
Y0u <- Y0 %*% A; Y1u <- Y1 %*% A
for (s in c("identity","diagvar","fullcov")) {
  cat(sprintf("  S0=%-9s H0 : %10.4f -> %12.4f | H1 : %10.4f -> %12.4f\n", s,
     bf_manova(Y0,f,S0=s), bf_manova(Y0u,f,S0=s), bf_manova(Y1,f,S0=s), bf_manova(Y1u,f,S0=s)))
}
cat("  Transformation non diagonale (rotation+echelle), S0=fullcov et diagvar :\n")
set.seed(105); R <- matrix(rnorm(9),3)
cat(sprintf("  fullcov : %.6f -> %.6f ; diagvar : %.6f -> %.6f\n",
   bf_manova(Y1,f,S0="fullcov"), bf_manova(Y1 %*% R,f,S0="fullcov"),
   bf_manova(Y1,f,S0="diagvar"), bf_manova(Y1 %*% R,f,S0="diagvar")))

cat("\n=== G. Etalonnage : 200 jeux, p=3, 3 groupes, n=20/gr, S0=fullcov, nu0=p+2, lam=1/3 et 1 ===\n")
set.seed(106)
sim <- function(effet, lam) {
  t(replicate(200, {
    Y <- rmvn(60,Sig) + effet*mu[as.integer(f),]
    lb <- bf_manova(Y,f,lam=lam,S0="fullcov")
    p <- summary(manova(Y~f),test="Pillai")$stats[1,6]
    c(P=plogis(lb), p=p)
  }))
}
for (lam in c(1/3,1)) for (eff in c(0,1)) {
  r <- sim(eff,lam)
  cat(sprintf("  lam=%.2f effet=%d : P(H1|D) mediane=%.3f ; >0.5 : %.3f ; >0.95 : %.3f ; manova p<0.05 : %.3f\n",
      lam, eff, median(r[,1]), mean(r[,1]>.5), mean(r[,1]>.95), mean(r[,2]<.05)))
}
t0 <- Sys.time(); for (i in 1:1000) bf_manova(Y1,f,S0="fullcov"); cat("\n  Temps : ", format(as.numeric(difftime(Sys.time(),t0,units="secs")), digits=3), "ms par BF\n")
cat("\nFIN\n"); print(sessionInfo())

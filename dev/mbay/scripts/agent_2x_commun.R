# agent_2x_commun.R : fonctions communes aux simulations de calibration de m.bay
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressPackageStartupMessages({library(BayesFactor)})
options(BFprogress = FALSE)
RES <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats"

# ---------------------------------------------------------------------------
# FORME EXACTE DU CODAGE (lue dans la source de BayesFactor 0.9.12-4.8) :
#   BayesFactor:::fixedFromRandomProjection(a) :
#       centering <- diag(a) - 1/a
#       P <- eigen(centering)$vectors[, 1:(a-1)]      # matrice a x (a-1)
#   BayesFactor:::oneDesignMatrix : pour un facteur fixe,
#       X <- model.Matrix(~ A - 1) %*% P              # indicatrices x P
#   pour une interaction, X_AB <- rowMultiply(X_A, X_B) (produit de Kronecker ligne a ligne
#   des deux matrices DEJA reduites). Donc P'P = I_(a-1) et PP' = I_a - 1/a ;
#   effets de groupe alpha = P b, de somme nulle, ||alpha|| = ||b||.
#   L'a priori est invariant par rotation de b : le signe/choix des vecteurs propres est sans effet.
#   b | s2, g ~ N(0, s2 g I), g ~ InvGamma(1/2, r^2/2), g COMMUN aux a-1 coordonnees d'un terme,
#   un g par terme. Le BF est invariant par y -> c y + d : generer avec mu = 0, s = 1 suffit.
#   Pour k = 2 : P = (1/sqrt2, -1/sqrt2)', alpha = +-b/sqrt2, delta = sqrt2 b,
#   d'ou r(delta) = sqrt2 * 1/2 = sqrt2/2 (ttestBF, "medium").
# ---------------------------------------------------------------------------
proj <- function(a) BayesFactor:::fixedFromRandomProjection(a)   # exactement celle du paquet
rg <- function(r) 1 / rgamma(1, shape = 1/2, rate = r^2 / 2)     # g ~ InvGamma(1/2, r^2/2)

# Tirage des effets de groupe sous l'a priori (k groupes)
ralpha <- function(k, r = 1/2) { g <- rg(r); b <- rnorm(k - 1, 0, sqrt(g)); as.numeric(as.matrix(proj(k)) %*% b) }

# log BF10 d'un facteur : k = 2 -> ttestBF ; k >= 3 -> anovaBF (quadrature pour un seul g)
lbf_anova <- function(y, f, r = 1/2) {
  d <- data.frame(y = y, A = f)
  res <- tryCatch(anovaBF(y ~ A, data = d, rscaleFixed = r, progress = FALSE)@bayesFactor$bf,
                  error = function(e) NA_real_)
  res
}
lbf_ttest <- function(y, f, r = sqrt(2)/2) {
  l <- levels(f)
  tryCatch(ttestBF(x = y[f == l[1]], y = y[f == l[2]], rscale = r, progress = FALSE)@bayesFactor$bf,
           error = function(e) NA_real_)
}
lbf <- function(y, f) if (nlevels(f) == 2) lbf_ttest(y, f) else lbf_anova(y, f)
prob <- function(l, pi) plogis(l + qlogis(pi))       # P(H1|D), calcul sur l'echelle log
paov <- function(y, f) summary(aov(y ~ f))[[1]][1, 5]

wilson <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA)); z <- qnorm(1 - (1 - conf)/2); p <- x/n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}

# Repli numerique (quand integrate() de BayesFactor echoue, F tres grand) : meme integrale, en log g, centree sur le pic.
# BF(g) = (1+c g)^((N-1-q)/2) / (1 + c g (1-R2))^((N-1)/2), avec X'X centree = c I (plan equilibre).
lbf_repli <- function(R2, N, q, c, r = 1/2) {
  f <- function(lg) { g <- exp(lg); (N-1-q)/2*log1p(c*g) - (N-1)/2*log1p(c*g*(1-R2)) +
      .5*log(r^2/2) - lgamma(.5) - 1.5*lg - r^2/(2*g) + lg }
  o <- optimize(f, c(-40, 80), maximum = TRUE)
  v <- integrate(function(x) exp(f(x) - o$objective), o$maximum - 30, o$maximum + 30, subdivisions = 2000L, stop.on.error = FALSE)$value
  log(v) + o$objective
}
# Voie rapide : memes fonctions internes que anovaBF/ttestBF (plan equilibre). Compteur de replis dans .rep$n
.rep <- new.env(); .rep$n <- 0L
Fstat <- function(y, n, k) {            # y range par groupe, n par groupe
  m <- colMeans(matrix(y, n, k)); gm <- mean(m)
  ssb <- n*sum((m-gm)^2); sse <- sum((y - rep(m, each = n))^2)
  (ssb/(k-1)) / (sse/(k*(n-1)))
}
lbf_F <- function(F, n, k, r = 1/2) {   # r sur l'echelle reduite ; k=2 -> ttest.tstat avec r*sqrt2
  if (k == 2) l <- tryCatch(ttest.tstat(sqrt(F), n, n, rscale = r*sqrt(2))$bf, error = function(e) NA_real_)
  else l <- tryCatch(suppressWarnings(oneWayAOV.Fstat(F, n, k, rscale = r)$bf), error = function(e) NA_real_)
  if (is.na(l) || !is.finite(l)) { .rep$n <- .rep$n + 1L
    N <- n*k; q <- k-1; l <- lbf_repli(q*F/(q*F+N-1-q), N, q, n, r) }
  l
}
pF <- function(F, n, k) pf(F, k-1, k*(n-1), lower.tail = FALSE)

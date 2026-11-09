# agent_51 : reference independante des facteurs de Bayes de BayesFactor::lmBF
# pour des plans 2x2 (A, B, A + B, A + B + A:B), effets fixes, rscaleFixed = 0.5.
#
# Modele (Rouder, Morey, Speckman & Province 2012, J. Math. Psych. 56:356-374) :
#   y = mu 1 + sum_t X_t theta_t + eps, eps ~ N(0, sigma^2 I), p(mu, sigma^2) ~ 1/sigma^2
#   theta_t | g_t ~ N(0, g_t sigma^2 I), g_t ~ InvGamma(1/2, r^2/2), r = 0.5.
# Codage : facteur a 2 niveaux -> colonne +-1/sqrt(2) (vecteur propre unitaire de
# I - J/2, base orthonormee du complement de 1) ; interaction -> produit des colonnes.
# Le codage est construit ici a la main, sans appeler BayesFactor.
#
# A g fixes (mu, sigma^2, theta integres), avec y_c et X_c centres, n obs :
#   log BF(g) = -1/2 log|G| - 1/2 log|M| - (n-1)/2 log(1 - b' M^-1 b / y_c'y_c),
#   M = X_c'X_c + G^-1, b = X_c'y_c.
# Puis integration sur u = log g par trapezes (1D, 2D, 3D), grille fine contre grossiere,
# deux bornes d'integration ; second controle par Monte-Carlo naif dans l'a priori des g.

.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
sortie <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats/agent_51.txt"

jeu <- function(j) {
  set.seed(5000 + j)
  n <- list(c(5,5,5,5), c(3,7,4,9), c(10,10,10,10))[[j]]   # cases (a1b1, a2b1, a1b2, a2b2)
  A <- factor(rep(c("a1","a2","a1","a2"), n)); B <- factor(rep(c("b1","b1","b2","b2"), n))
  y <- rnorm(sum(n), 0.6 * (A == "a2") + 0.4 * (B == "b2") + c(0, 0, 0, 0.8)[as.integer(interaction(A, B))])
  data.frame(y, A, B)
}

r <- 0.5
s <- r^2 / 2   # echelle de l'InvGamma(1/2, s)

codage <- function(d) {
  xa <- ifelse(d$A == "a1", 1, -1) / sqrt(2)
  xb <- ifelse(d$B == "b1", 1, -1) / sqrt(2)
  list(A = xa, B = xb, "A:B" = xa * xb)
}

lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

# log BF(g) vectorise ; G est une liste de vecteurs (un par terme, meme longueur)
logbf_g <- function(XtX, b, yy, n, G) {
  p <- length(G)
  if (p == 1) {
    M <- XtX[1,1] + 1/G[[1]]
    det <- M; q <- b[1]^2 / M
  } else if (p == 2) {
    a <- XtX[1,1] + 1/G[[1]]; c2 <- XtX[2,2] + 1/G[[2]]; d <- XtX[1,2]
    det <- a*c2 - d^2
    q <- (b[1]^2*c2 + b[2]^2*a - 2*b[1]*b[2]*d) / det
  } else {
    a <- XtX[1,1] + 1/G[[1]]; bb <- XtX[2,2] + 1/G[[2]]; cc <- XtX[3,3] + 1/G[[3]]
    d <- XtX[1,2]; e <- XtX[1,3]; f <- XtX[2,3]
    A11 <- bb*cc - f^2; A22 <- a*cc - e^2; A33 <- a*bb - d^2
    A12 <- e*f - d*cc; A13 <- d*f - bb*e; A23 <- d*e - a*f
    det <- a*A11 + d*A12 + e*A13
    q <- (b[1]^2*A11 + b[2]^2*A22 + b[3]^2*A33 + 2*b[1]*b[2]*A12 + 2*b[1]*b[3]*A13 + 2*b[2]*b[3]*A23) / det
  }
  -0.5 * Reduce(`+`, lapply(G, log)) - 0.5 * log(det) - (n - 1)/2 * log(1 - q / yy)
}

# log densite a priori de u = log g
lprior_u <- function(u) 0.5*log(s) - 0.5*log(pi) - u/2 - s*exp(-u)

prep <- function(d, termes) {
  X <- do.call(cbind, codage(d)[termes])
  Xc <- sweep(X, 2, colMeans(X)); yc <- d$y - mean(d$y)
  list(XtX = crossprod(Xc), b = drop(crossprod(Xc, yc)), yy = sum(yc^2), n = nrow(d))
}

# quadrature par trapezes sur [lo, hi]^p avec N points par dimension
quad <- function(P, p, N, lo, hi) {
  u <- seq(lo, hi, length.out = N); h <- u[2] - u[1]
  w <- log(c(0.5, rep(1, N - 2), 0.5) * h) + lprior_u(u)
  g <- exp(u)
  if (p == 1) return(lse(w + logbf_g(P$XtX, P$b, P$yy, P$n, list(g))))
  if (p == 2) {
    g1 <- rep(g, N); g2 <- rep(g, each = N); W <- rep(w, N) + rep(w, each = N)
    return(lse(W + logbf_g(P$XtX, P$b, P$yy, P$n, list(g1, g2))))
  }
  g1 <- rep(g, N); g2 <- rep(g, each = N); W2 <- rep(w, N) + rep(w, each = N)
  parts <- vapply(seq_len(N), function(k)
    lse(W2 + w[k] + logbf_g(P$XtX, P$b, P$yy, P$n, list(g1, g2, rep(g[k], N*N)))), 0)
  lse(parts)
}

mc <- function(P, p, M = 1e6) {
  G <- lapply(seq_len(p), function(i) 1 / rgamma(M, shape = 0.5, rate = s))  # InvGamma(1/2, s)
  l <- logbf_g(P$XtX, P$b, P$yy, P$n, G)
  est <- lse(l) - log(M)
  w <- exp(l - max(l)); rel <- sd(w) / mean(w) / sqrt(M)
  c(est = est, relerr = rel)
}

modeles <- list("A" = "A", "B" = "B", "A + B" = c("A","B"), "A + B + A:B" = c("A","B","A:B"))
fml <- list("A" = y ~ A, "B" = y ~ B, "A + B" = y ~ A + B, "A + B + A:B" = y ~ A + B + A:B)
Ns <- c(`1` = 4000, `2` = 600, `3` = 200)   # points fins par dimension

suppressMessages(library(BayesFactor))
res <- list()
for (j in 1:3) {
  d <- jeu(j)
  for (m in names(modeles)) {
    tt <- modeles[[m]]; p <- length(tt); P <- prep(d, tt); N <- Ns[[as.character(p)]]
    fin   <- quad(P, p, N, -12, 30)
    gros  <- quad(P, p, N %/% 2, -12, 30)
    large <- quad(P, p, N, -18, 40)
    set.seed(100 + j)
    mcv <- mc(P, p, if (p == 3) 4e5 else 1e6)
    err <- max(abs(fin - gros), abs(fin - large))
    set.seed(1)
    bfo <- extractBF(lmBF(fml[[m]], data = d, rscaleFixed = 0.5, iterations = 100000, progress = FALSE), logbf = TRUE)
    res[[length(res) + 1]] <- data.frame(jeu = j, modele = m, ref = fin, err_quad = err,
      ecart_gros = fin - gros, ecart_bornes = fin - large,
      mc = mcv[["est"]], mc_relerr = mcv[["relerr"]],
      bf_BayesFactor = bfo$bf, error_BayesFactor = bfo$error, ecart = bfo$bf - fin)
    cat(j, m, "fait\n")
  }
}
tab <- do.call(rbind, res)

# controle 1D supplementaire : integrate() en g, sur le modele A
ctl <- sapply(1:3, function(j) {
  P <- prep(jeu(j), "A")
  f <- function(g) exp(logbf_g(P$XtX, P$b, P$yy, P$n, list(g)) + 0.5*log(s) - 0.5*log(pi) - 1.5*log(g) - s/g)
  log(integrate(f, 0, Inf, rel.tol = 1e-10)$value)
})

sink(sortie)
cat("agent_51 - reference independante (quadrature en log g) contre BayesFactor::lmBF\n")
cat("BayesFactor", as.character(packageVersion("BayesFactor")), "-", R.version.string, "-", format(Sys.time()), "\n")
cat("r = 0.5, g_t ~ InvGamma(1/2, r^2/2), codage +-1/sqrt(2), interaction = produit (+-1/2)\n")
cat("Points par dimension (fin) : 1D 4000, 2D 600, 3D 200 ; grossier = moitie ; bornes u in [-12,30] contre [-18,40]\n")
cat("log BF contre le modele constant ; ecart = BayesFactor - reference\n\n")
op <- options(width = 200)
print(within(tab, { ref <- round(ref, 6); mc <- round(mc, 4); bf_BayesFactor <- round(bf_BayesFactor, 6); ecart <- round(ecart, 6) }), row.names = FALSE, digits = 6)
cat("\nControle 1D par integrate() (modele A) :", sprintf("%.8f", ctl), "\n")
cat("Ecart integrate - trapezes (modele A) :", sprintf("%.2e", ctl - tab$ref[tab$modele == "A"]), "\n")
options(op)
sink()
cat(readLines(sortie), sep = "\n")

# Controle complementaire : dispersion de lmBF sur 10 graines (set.seed(1..10)),
# pour les modeles a 2 et 3 termes, afin de juger si les ecarts relevent du bruit
# Monte-Carlo de BayesFactor.
disp <- list()
for (j in 1:3) {
  d <- jeu(j)
  for (m in c("A + B", "A + B + A:B")) {
    v <- sapply(1:10, function(sd) {
      set.seed(sd)
      extractBF(lmBF(fml[[m]], data = d, rscaleFixed = 0.5, iterations = 100000, progress = FALSE), logbf = TRUE)$bf
    })
    ref <- tab$ref[tab$jeu == j & tab$modele == m]
    disp[[length(disp) + 1]] <- data.frame(jeu = j, modele = m, ref = ref, moy_BF10graines = mean(v),
      sd_BF10graines = sd(v), ecart_moy = mean(v) - ref, ecart_moy_sur_se = (mean(v) - ref) / (sd(v) / sqrt(10)),
      min = min(v), max = max(v))
  }
}
sink(sortie, append = TRUE)
cat("\nDispersion de lmBF (iterations = 1e5) sur set.seed(1..10) :\n")
print(do.call(rbind, disp), row.names = FALSE, digits = 5)
sink()
cat(readLines(sortie), sep = "\n")

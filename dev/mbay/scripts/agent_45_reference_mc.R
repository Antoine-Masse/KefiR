# agent_45 : valeurs de reference de vraisemblances marginales, code independant.
# Ecrit a partir de la seule definition mathematique du modele.
#
# Usage : Rscript agent_45_reference_mc.R <mode> [arguments]
#   controle                 controles a solution connue
#   quad                     quadratures deterministes (normal_het 3D, t_het H0)
#   mc <jeu> <loi> <Mtotal>  Monte-Carlo par lots de 10^6, reprenable
#   naifn <jeu> <M>          controle : normal, a tire, tau par quadrature
#   naift <jeu> <M>          controle : t, a tire, (mu, log tau) par quadrature
#   table                    agrege les lots dans resultats/agent_45.txt
#
# ESTIMATEUR PRINCIPAL. Tirages dans l'a priori : theta~Dirichlet(u), nu-2~Exp(1),
# g~InvGamma, et pour la loi t lambda_ij~Gamma(nu/2,nu/2), representation de la
# Student en melange d'echelles. Sachant ces tirages, y est gaussien de variance
# (c_i/tau)*(kappa_ij), kappa_ij=(nu-2)/(nu*lambda_ij), et l'effet a_i s'integre
# analytiquement (variance (c_i/tau)*g, meme facteur tau), puis mu et tau aussi :
#   V_i = c_i (D_i + g_i 1 1'),   Cov(y)=V/tau
#   m(y | tirage) = (2 pi)^{-(N-1)/2} |V|^{-1/2} A^{-1/2} Gamma((N-1)/2) (Q/2)^{-(N-1)/2}
#   A = 1'V^{-1}1,  Q = y'V^{-1}y - (1'V^{-1}y)^2/A
# Aucun ajustement sur la position du posterieur.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "controle"

DIR_RES <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats"
DIR_CAS <- file.path(DIR_RES, "agent_45_cas")
if (!dir.exists(DIR_CAS)) dir.create(DIR_CAS, recursive = TRUE)

U <- 1.5; R0 <- 0.5; K <- 3L
FORMES <- c("H0", "tous", "s1", "s2", "s3")

donnees <- function(j) {
  if (j == 1) { set.seed(4501); y <- c(rnorm(4,0,1), rnorm(5,0.3,1), rnorm(4,0,2)); g <- rep(1:3,c(4,5,4)) }
  if (j == 2) { set.seed(4502); y <- c(rnorm(6,0,1), rnorm(6,2,1), rnorm(6,0,1)); g <- rep(1:3,each=6) }
  if (j == 3) { set.seed(4503); y <- c(rnorm(5,0,0.5), rnorm(5,0,1), rnorm(5,1,3)); g <- rep(1:3,each=5) }
  list(y = y, g = g, N = length(y), idx = split(seq_along(y), g), n = as.integer(table(g)))
}

lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Noyau : W,S,TT,LD matrices M x k (sommes par groupe ponderees par p=1/kappa, somme log kappa),
# cc matrice M x k des c_i, gm matrice M x k (ou scalaire 0) des g_i. Rend log m(y | tirage).
logm_noyau <- function(W, S, TT, LD, cc, gm, N, nv) {
  M <- nrow(W)
  if (length(gm) == 1L) gm <- matrix(gm, M, ncol(W))
  den <- 1 + gm * W
  A <- rowSums(W / (cc * den))
  B <- rowSums(S / (cc * den))
  C <- rowSums((TT - gm * S^2 / den) / cc)
  ldV <- rowSums(sweep(log(cc), 2, nv, `*`) + LD + log(den))
  Q <- C - B^2 / A
  if (any(Q <= 0)) Q[Q <= 0] <- NA
  -(N - 1) / 2 * log(2 * pi) - 0.5 * ldV - 0.5 * log(A) +
    lgamma((N - 1) / 2) - (N - 1) / 2 * log(Q / 2)
}

gm_forme <- function(forme, g, M) {
  gm <- matrix(0, M, K)
  if (forme == "tous") gm[] <- g
  if (forme %in% c("s1", "s2", "s3")) gm[, as.integer(substr(forme, 2, 2))] <- 2 * g
  gm   # 2*g ~ InvGamma(1/2, (r*sqrt2)^2/2) si g ~ InvGamma(1/2, r^2/2)
}

tire_c <- function(M) {
  G <- matrix(rgamma(M * K, U), M, K)
  1 / (K * G / rowSums(G))
}

stats_groupes <- function(d, p) {
  M <- nrow(p); W <- S <- TT <- LD <- matrix(0, M, K)
  for (i in seq_len(K)) {
    id <- d$idx[[i]]; pi_ <- p[, id, drop = FALSE]
    W[, i] <- rowSums(pi_)
    S[, i] <- drop(pi_ %*% d$y[id])
    TT[, i] <- drop(pi_ %*% d$y[id]^2)
    LD[, i] <- -rowSums(log(pi_))
  }
  list(W = W, S = S, TT = TT, LD = LD)
}

stats_normal <- function(d, M) {
  W <- S <- TT <- matrix(0, M, K)
  for (i in seq_len(K)) {
    id <- d$idx[[i]]
    W[, i] <- length(id); S[, i] <- sum(d$y[id]); TT[, i] <- sum(d$y[id]^2)
  }
  list(W = W, S = S, TT = TT, LD = matrix(0, M, K))
}

# Un lot de M tirages : rend data.frame forme, M, logsum, logsum2, maxlog, nb_na
un_lot <- function(d, loi, M, taille = 100000L) {
  acc <- lapply(FORMES, function(f) list(ls = numeric(0), ls2 = numeric(0), mx = -Inf, na = 0))
  names(acc) <- FORMES
  reste <- M
  while (reste > 0) {
    m <- min(reste, taille); reste <- reste - m
    cc <- tire_c(m)
    g <- 1 / rgamma(m, 0.5, rate = R0^2 / 2)
    if (loi == "t_het") {
      e <- rexp(m); nu <- 2 + e
      lam <- matrix(rgamma(m * d$N, shape = rep(nu / 2, d$N), rate = rep(nu / 2, d$N)), m, d$N)
      st <- stats_groupes(d, lam * (nu / e))
      rm(lam)
    } else st <- stats_normal(d, m)
    for (f in FORMES) {
      lw <- logm_noyau(st$W, st$S, st$TT, st$LD, cc, gm_forme(f, g, m), d$N, d$n)
      acc[[f]]$na <- acc[[f]]$na + sum(is.na(lw))
      lw <- lw[!is.na(lw)]
      acc[[f]]$ls <- c(acc[[f]]$ls, lse(lw))
      acc[[f]]$ls2 <- c(acc[[f]]$ls2, lse(2 * lw))
      acc[[f]]$mx <- max(acc[[f]]$mx, max(lw))
    }
  }
  do.call(rbind, lapply(FORMES, function(f) data.frame(
    forme = f, M = M, logsum = lse(acc[[f]]$ls), logsum2 = lse(acc[[f]]$ls2),
    maxlog = acc[[f]]$mx, nb_na = acc[[f]]$na)))
}

fichier_lot <- function(j, loi, b) file.path(DIR_CAS, sprintf("j%d_%s_b%05d.csv", j, loi, b))

ecrit_atomique <- function(df, f) {
  tmp <- paste0(f, ".tmp")
  write.csv(df, tmp, row.names = FALSE)
  file.rename(tmp, f)
}

lance_mc <- function(j, loi, Mtotal, Mlot = 1000000L) {
  d <- donnees(j)
  nb <- ceiling(Mtotal / Mlot)
  for (b in seq_len(nb)) {
    f <- fichier_lot(j, loi, b)
    if (file.exists(f)) next
    set.seed(45000000L + 1000000L * j + 100000L * (loi == "t_het") + b)
    t0 <- Sys.time()
    r <- un_lot(d, loi, Mlot)
    ecrit_atomique(r, f)
    cat(sprintf("j%d %s lot %d : %.1f s\n", j, loi, b, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()
  }
}

agrege <- function(j, loi) {
  fs <- list.files(DIR_CAS, pattern = sprintf("^j%d_%s_b[0-9]+\\.csv$", j, loi), full.names = TRUE)
  if (!length(fs)) return(NULL)
  x <- do.call(rbind, lapply(fs, read.csv))
  do.call(rbind, lapply(FORMES, function(f) {
    z <- x[x$forme == f, ]
    M <- sum(z$M); ls <- lse(z$logsum); ls2 <- lse(z$logsum2)
    ess <- exp(2 * ls - ls2)
    # se delta sur le log, et se par moyennes de lots
    se <- sqrt(max(1 / ess - 1 / M, 0))
    mb <- exp(z$logsum - log(z$M) - (ls - log(M)))
    data.frame(jeu = j, loi = loi, forme = f, log_m = ls - log(M), se = se, ess = ess, M = M,
               se_lots = if (length(mb) > 1) sd(mb) / sqrt(length(mb)) else NA,
               part_max = exp(max(z$maxlog) - ls), nb_na = sum(z$nb_na))
  }))
}

# ---- Quadrature deterministe normal_het : (v1, v2, log g), a, mu, tau analytiques ----
# theta1=v1, theta2=(1-v1)v2, theta3=(1-v1)(1-v2), v1~Beta(u,2u), v2~Beta(u,u).
quad_normal <- function(d, h = 0.2, hg = 0.1, xmax = 16) {
  x <- seq(-xmax, xmax, by = h)
  v <- plogis(x)
  lp1 <- dbeta(v, U, 2 * U, log = TRUE) + log(v) + log1p(-v) + log(h)
  lp2 <- dbeta(v, U, U, log = TRUE) + log(v) + log1p(-v) + log(h)
  ii <- expand.grid(a = seq_along(x), b = seq_along(x))
  v1 <- v[ii$a]; v2 <- v[ii$b]; lpv <- lp1[ii$a] + lp2[ii$b]
  # 1-v pour v proche de 1 : plogis(-x)
  w <- plogis(-x); w1 <- w[ii$a]; w2 <- w[ii$b]
  cc <- 1 / (K * cbind(v1, w1 * v2, w1 * w2))
  M <- nrow(cc); st <- stats_normal(d, M)
  lg <- seq(-14, 60, by = hg)
  b <- R0^2 / 2
  lpg <- 0.5 * log(b) - lgamma(0.5) - 1.5 * lg - b / exp(lg) + lg + log(hg)
  res <- numeric(length(FORMES)); names(res) <- FORMES
  res["H0"] <- lse(lpv + logm_noyau(st$W, st$S, st$TT, st$LD, cc, 0, d$N, d$n))
  for (f in FORMES[-1]) {
    tmp <- numeric(length(lg))
    for (q in seq_along(lg)) {
      gm <- gm_forme(f, exp(lg[q]), M)
      tmp[q] <- lse(lpv + logm_noyau(st$W, st$S, st$TT, st$LD, cc, gm, d$N, d$n))
    }
    res[f] <- lse(tmp + lpg)
  }
  res
}

# ---- Quadrature deterministe H0 par factorisation ----
# Avec omega_i = k tau theta_i : dtau/tau Dir(theta) dtheta =
#   (1/Gamma(u)^k) int dx x^{ku-1} prod_i omega_i^{u-1} e^{-x omega_i} domega_i.
# Pour t, phi_i=omega_i nu/(nu-2) (s_i^2=1/phi_i), le facteur nu/(nu-2) s'elimine en
# changeant x. D'ou m = int p(nu) dnu int dmu int dx x^{ku-1}/Gamma(u)^k
#   prod_i int dphi phi^{u-1} e^{-x phi} prod_j t_nu-echelle(y_ij; mu, 1/phi).
quad_h0 <- function(d, loi, hl = 0.2, nmu = 301, he = 0.25) {
  sdy <- sd(d$y); yb <- mean(d$y)
  tt <- seq(-asinh(400), asinh(400), length.out = nmu); dt <- tt[2] - tt[1]
  mu <- yb + sdy * sinh(tt); ldmu <- log(sdy) + log(cosh(tt)) + log(dt)
  ell <- seq(-16, 16, by = hl) - 2 * log(sdy); phi <- exp(ell)
  lx <- seq(-16, 18, by = hl) - 2 * log(sdy) * (-1); x <- exp(lx)
  E <- exp(U * ell %o% rep(1, length(x)) - phi %o% x) * hl
  val <- function(nu) {
    tot <- matrix(0, length(mu), length(x))
    for (i in seq_len(K)) {
      LL <- matrix(0, length(mu), length(phi))
      for (y1 in d$y[d$idx[[i]]]) {
        dd <- outer((y1 - mu)^2, phi)
        if (is.null(nu)) LL <- LL + 0.5 * rep(ell, each = length(mu)) - 0.5 * log(2 * pi) - 0.5 * dd
        else LL <- LL + lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi) +
          0.5 * rep(ell, each = length(mu)) - (nu + 1) / 2 * log1p(dd / nu)
      }
      mx <- max(LL)
      I <- exp(LL - mx) %*% E
      tot <- tot + mx + log(I)
    }
    tot <- tot + K * U * rep(lx, each = length(mu)) - K * lgamma(U) + log(hl) + ldmu
    lse(as.vector(tot))
  }
  if (loi == "normal_het") return(val(NULL))
  le <- seq(-30, 4.6, by = he)
  v <- vapply(le, function(l) val(2 + exp(l)), 0)
  lse(v + le - exp(le) + log(he))
}

# ---- Integrale directe sur (mu, log phi) de la vraisemblance t, chocs standardises donnes ----
# s_i^2=c_i/phi, phi=tau nu/(nu-2) ; z=(y-mu)/s_i - delta_i sqrt(g nu/(nu-2)) ; dtau/tau=dphi/phi.
grille_t <- function(d, nmu = 121, hl = 0.2) {
  sdy <- sd(d$y)
  tt <- seq(-asinh(60), asinh(60), length.out = nmu)
  ell <- seq(-8, 18, by = hl) - 2 * log(sdy)
  list(mu = mean(d$y) + sdy * sinh(tt), lw = log(sdy) + log(cosh(tt)) + log(tt[2] - tt[1]) ,
       ell = ell, hl = hl)
}
direct_t <- function(d, gr, cc, nu, b) {   # b : effet standardise a_i/s_i, vecteur k
  lc <- lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi)
  nm <- length(gr$mu); np <- length(gr$ell)
  LL <- matrix(0, nm, np)
  for (i in seq_len(K)) {
    inv_s <- exp(gr$ell / 2) / sqrt(cc[i])
    for (y1 in d$y[d$idx[[i]]]) {
      z <- outer(y1 - gr$mu, inv_s) - b[i]
      LL <- LL - (nu + 1) / 2 * log1p(z * z / nu)
    }
    LL <- LL + rep(d$n[i] * (lc + log(inv_s)), each = nm)
  }
  LL <- LL + gr$lw + log(gr$hl)
  lse(as.vector(LL))
}

# ---- Normal, a tire : mu analytique, quadrature sur log tau ----
naif_normal <- function(d, M, taille = 5000L, hl = 0.1) {
  ell <- seq(-14, 26, by = hl) ; tau <- exp(ell); np <- length(ell)
  yb <- vapply(d$idx, function(id) mean(d$y[id]), 0)
  y2 <- vapply(d$idx, function(id) sum(d$y[id]^2), 0)
  out <- setNames(vector("list", 5), FORMES)
  for (f in FORMES) out[[f]] <- list(ls = numeric(0), ls2 = numeric(0))
  reste <- M
  while (reste > 0) {
    m <- min(reste, taille); reste <- reste - m
    cc <- tire_c(m); g <- 1 / rgamma(m, 0.5, rate = R0^2 / 2)
    delta <- matrix(rnorm(m * K), m, K)
    for (f in FORMES) {
      bb <- delta * sqrt(gm_forme(f, g, m))
      A <- B <- C <- P <- matrix(0, m, np)
      for (i in seq_len(K)) {
        wi <- outer(1 / cc[, i], tau); rw <- sqrt(wi); ni <- d$n[i]
        A <- A + ni * wi
        B <- B + ni * (wi * yb[i] - bb[, i] * rw)
        C <- C + wi * y2[i] - 2 * bb[, i] * rw * ni * yb[i] + ni * bb[, i]^2
        P <- P + ni / 2 * (log(wi) - log(2 * pi))
      }
      L <- P + 0.5 * (log(2 * pi) - log(A)) - 0.5 * (C - B^2 / A)
      mx <- do.call(pmax, as.data.frame(L))
      lw <- mx + log(rowSums(exp(L - mx))) + log(hl)
      out[[f]]$ls <- c(out[[f]]$ls, lse(lw)); out[[f]]$ls2 <- c(out[[f]]$ls2, lse(2 * lw))
    }
  }
  do.call(rbind, lapply(FORMES, function(f) {
    ls <- lse(out[[f]]$ls); ess <- exp(2 * ls - lse(out[[f]]$ls2))
    data.frame(forme = f, log_m = ls - log(M), se = sqrt(1 / ess - 1 / M), ess = ess, M = M)
  }))
}

# ---- t, a tire : quadrature (mu, log phi) ----
naif_t <- function(d, M, formes = FORMES, ...) {
  gr <- grille_t(d, ...)
  lw <- matrix(NA_real_, M, length(formes), dimnames = list(NULL, formes))
  for (q in seq_len(M)) {
    cc <- as.vector(tire_c(1L)); g <- 1 / rgamma(1, 0.5, rate = R0^2 / 2)
    e <- rexp(1); nu <- 2 + e; delta <- rnorm(K)
    for (f in formes) {
      b <- delta * sqrt(as.vector(gm_forme(f, g, 1L)) * nu / e)
      lw[q, f] <- direct_t(d, gr, cc, nu, b)
    }
  }
  do.call(rbind, lapply(formes, function(f) {
    ls <- lse(lw[, f]); ess <- exp(2 * ls - lse(2 * lw[, f]))
    data.frame(forme = f, log_m = ls - log(M), se = sqrt(max(1 / ess - 1 / M, 0)), ess = ess, M = M)
  }))
}

# ---- Controles ----
controle <- function() {
  cat("== 1. Normal a variance commune, H0 : formule fermee contre quadrature 1D et 2D ==\n")
  for (j in 1:3) {
    d <- donnees(j); N <- d$N; sc <- sum((d$y - mean(d$y))^2)
    # formule fermee : mu et tau integres
    ferme <- -(N - 1) / 2 * log(2 * pi) - 0.5 * log(N) + lgamma((N - 1) / 2) - (N - 1) / 2 * log(sc / 2)
    # une integrale sur log tau, mu integre a la main
    f1 <- function(l) { t <- exp(l); (t / (2 * pi))^(N / 2) * sqrt(2 * pi / (N * t)) * exp(-t * sc / 2) }
    # bornes centrees sur le mode : sur [-30,30], integrate() manque le pic (jeux 2 et 3)
    lc <- log((N - 1) / sc)
    i1 <- log(integrate(f1, lc - 6, lc + 6, rel.tol = 1e-10, subdivisions = 2000L)$value)
    # noyau du script avec c=1
    st <- stats_normal(d, 1L)
    nq <- logm_noyau(st$W, st$S, st$TT, st$LD, matrix(1, 1, K), 0, N, d$n)
    # 2D brute force sur (mu, log tau)
    mu <- seq(mean(d$y) - 30, mean(d$y) + 30, length.out = 6001); l <- seq(-10, 10, by = 0.01)
    S <- vapply(mu, function(m) sum((d$y - m)^2), 0)
    LL <- outer(S, exp(l)) * (-0.5) + rep(N / 2 * (l - log(2 * pi)), each = length(mu))
    b2 <- lse(as.vector(LL)) + log(mu[2] - mu[1]) + log(0.01)
    cat(sprintf(" jeu %d : ferme %.8f | integrate %.8f | noyau %.8f | 2D %.8f\n", j, ferme, i1, nq, b2))
  }

  cat("\n== 2. Normal variance commune, H1 tous, g fixe : noyau contre integration directe de a ==\n")
  # k=3, c=1, g=0.7 : on integre a sur grille de Gauss produit via quadrature sur delta
  d <- donnees(1); g <- 0.7
  st <- stats_normal(d, 1L)
  nq <- logm_noyau(st$W, st$S, st$TT, st$LD, matrix(1, 1, K), matrix(g, 1, K), d$N, d$n)
  set.seed(1)
  dd <- naif_cond_normal(d, c(1, 1, 1), rep(g, 3), 400000L)
  cat(sprintf(" noyau %.5f | a tire + tau quadrature %.5f (se %.5f)\n", nq, dd[1], dd[2]))

  cat("\n== 3. Tirage de g et de la Dirichlet : moments empiriques ==\n")
  set.seed(2); M <- 2e6; cc <- tire_c(M); th <- 1 / (K * cc)
  cat(sprintf(" E[theta1] %.5f (1/3) ; Var[theta1] %.5f (exacte %.5f) ; somme theta max dev %.1e\n",
              mean(th[, 1]), var(th[, 1]), (1/3)*(2/3)/(3*U+1), max(abs(rowSums(th) - 1))))
  g <- 1 / rgamma(M, 0.5, rate = R0^2 / 2)
  # sqrt(g)*Z de loi de Cauchy d'echelle r : mediane de |.| = r
  cat(sprintf(" mediane |sqrt(g) Z| %.5f (r = %.3f, Cauchy d'echelle r) ; 2*g : %.5f (r*sqrt2 = %.5f)\n",
              median(abs(sqrt(g) * rnorm(M))), R0, median(abs(sqrt(2 * g) * rnorm(M))), R0 * sqrt(2)))

  cat("\n== 4. t, point fixe (c, nu, g) : lambda-MC + analytique contre delta-tire + quadrature 2D ==\n")
  pts <- list(list(cc = c(1, 1, 1), nu = 5, g = 0.4), list(cc = c(0.5, 1.2, 4), nu = 2.5, g = 1),
              list(cc = c(2, 0.6, 1.5), nu = 12, g = 0.1))
  for (j in 1:3) for (p in pts) {
    d <- donnees(j); gr <- grille_t(d)
    set.seed(10 + j)
    M <- 1000000L; e <- p$nu - 2
    lam <- matrix(rgamma(M * d$N, p$nu / 2, rate = p$nu / 2), M, d$N)
    st <- stats_groupes(d, lam * (p$nu / e)); rm(lam)
    ccm <- matrix(p$cc, M, K, byrow = TRUE)
    for (f in c("H0", "tous", "s2")) {
      gm <- matrix(as.vector(gm_forme(f, p$g, 1L)), M, K, byrow = TRUE)
      lw <- logm_noyau(st$W, st$S, st$TT, st$LD, ccm, gm, d$N, d$n)
      ls <- lse(lw); ess <- exp(2 * ls - lse(2 * lw))
      nd <- if (f == "H0") 1L else 1500L
      del <- matrix(rnorm(nd * K), nd, K); if (f == "H0") del[] <- 0
      ld <- vapply(seq_len(nd), function(q)
        direct_t(d, gr, p$cc, p$nu, del[q, ] * sqrt(gm[1, ] * p$nu / e)), 0)
      lsd <- lse(ld); essd <- exp(2 * lsd - lse(2 * ld))
      sed <- if (f == "H0") 0 else sqrt(1 / essd)
      cat(sprintf(" jeu %d nu=%.1f %-4s : lambda %.4f (se %.4f) | direct %.4f (se %.4f) | ecart/se %.2f\n",
                  j, p$nu, f, ls - log(M), sqrt(1 / ess), lsd - log(nd), sed,
                  (ls - log(M) - lsd + log(nd)) / sqrt(1 / ess + sed^2)))
    }
  }
  cat("\n== 5. Factorisation H0 (quad_h0, loi normale) contre quadrature 3D simplex ==\n")
  for (j in 1:3) {
    d <- donnees(j)
    cat(sprintf(" jeu %d : factorisation %.6f | simplex %.6f\n", j, quad_h0(d, "normal_het"),
                quad_normal_h0(d)))
  }
}

quad_normal_h0 <- function(d, h = 0.1, xmax = 16) {
  x <- seq(-xmax, xmax, by = h); v <- plogis(x); w <- plogis(-x)
  lp1 <- dbeta(v, U, 2 * U, log = TRUE) + log(v) + log(w) + log(h)
  lp2 <- dbeta(v, U, U, log = TRUE) + log(v) + log(w) + log(h)
  ii <- expand.grid(a = seq_along(x), b = seq_along(x))
  cc <- 1 / (K * cbind(v[ii$a], w[ii$a] * v[ii$b], w[ii$a] * w[ii$b]))
  st <- stats_normal(d, nrow(cc))
  lse(lp1[ii$a] + lp2[ii$b] + logm_noyau(st$W, st$S, st$TT, st$LD, cc, 0, d$N, d$n))
}

# normal, (c, g) fixes : a tire standardise, mu analytique, quadrature log tau
naif_cond_normal <- function(d, cc, gv, M, hl = 0.05) {
  ell <- seq(-12, 20, by = hl); tau <- exp(ell); np <- length(ell)
  yb <- vapply(d$idx, function(id) mean(d$y[id]), 0)
  y2 <- vapply(d$idx, function(id) sum(d$y[id]^2), 0)
  ls <- numeric(0); ls2 <- numeric(0); reste <- M
  while (reste > 0) {
    m <- min(reste, 5000L); reste <- reste - m
    A <- B <- C <- P <- matrix(0, m, np)
    for (i in seq_len(K)) {
      wi <- tau / cc[i]; rw <- sqrt(wi); ni <- d$n[i]; b <- rnorm(m) * sqrt(gv[i])
      A <- A + rep(ni * wi, each = m)
      B <- B + rep(ni * wi * yb[i], each = m) - outer(b, ni * rw)
      C <- C + rep(wi * y2[i], each = m) - 2 * outer(b, rw) * ni * yb[i] + ni * b^2
      P <- P + rep(ni / 2 * (log(wi) - log(2 * pi)), each = m)
    }
    L <- P + 0.5 * (log(2 * pi) - log(A)) - 0.5 * (C - B^2 / A)
    mx <- do.call(pmax, as.data.frame(L))
    lw <- mx + log(rowSums(exp(L - mx))) + log(hl)
    ls <- c(ls, lse(lw)); ls2 <- c(ls2, lse(2 * lw))
  }
  l <- lse(ls); ess <- exp(2 * l - lse(ls2))
  c(l - log(M), sqrt(1 / ess))
}

# ---- Sorties ----
ecrit_table <- function() {
  x <- list()
  for (j in 1:3) for (loi in c("t_het", "normal_het")) x[[length(x) + 1]] <- agrege(j, loi)
  x <- do.call(rbind, x)
  f <- file.path(DIR_RES, "agent_45.txt")
  cat("jeu loi forme log_m se ess\n", file = f)
  for (q in seq_len(nrow(x)))
    cat(sprintf("%d %s %s %.5f %.5f %.0f\n", x$jeu[q], x$loi[q], x$forme[q], x$log_m[q], x$se[q], x$ess[q]),
        file = f, append = TRUE)
  # controles en commentaires (lus comme tels par read.table)
  cat("# Quadrature deterministe (reference exacte a 1e-6 pres ; grille fine = grossiere) :\n", file = f, append = TRUE)
  for (j in 1:3) {
    fq <- file.path(DIR_CAS, sprintf("quad_j%d.csv", j))
    if (file.exists(fq)) { q <- read.csv(fq); for (r in seq_len(nrow(q)))
      cat(sprintf("# %d %s %s %.6f quad\n", q$jeu[r], q$loi[r], q$forme[r], q$fin[r]), file = f, append = TRUE) }
  }
  cat("# y pris tel que genere. Pour y standardise (y-moy)/sd : log_m + (N-1)*log(sd(y)) :\n", file = f, append = TRUE)
  for (j in 1:3) { d <- donnees(j)
    cat(sprintf("# jeu %d : N=%d sd=%.6f (N-1)*log(sd)=%.6f\n", j, d$N, sd(d$y), (d$N - 1) * log(sd(d$y))),
        file = f, append = TRUE) }
  print(x, digits = 6, row.names = FALSE)
}

# ---- Aiguillage ----
if (mode == "controle") controle()
if (mode == "mc") lance_mc(as.integer(args[2]), args[3], as.numeric(args[4]))
if (mode == "pilote") {
  for (j in 1:3) for (loi in c("t_het", "normal_het")) {
    d <- donnees(j); set.seed(99 + j); t0 <- Sys.time()
    r <- un_lot(d, loi, as.integer(args[2]))
    r$ess <- exp(2 * r$logsum - r$logsum2); r$log_m <- r$logsum - log(r$M)
    cat(sprintf("jeu %d %s (%.1f s)\n", j, loi, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(r[, c("forme", "log_m", "ess", "nb_na")], row.names = FALSE, digits = 6)
  }
}
if (mode == "quad") {
  out <- list()
  for (j in 1:3) {
    d <- donnees(j)
    fn <- file.path(DIR_CAS, sprintf("quad_j%d.csv", j))
    if (file.exists(fn)) next
    qn <- quad_normal(d); qn2 <- quad_normal(d, h = 0.1, hg = 0.05)
    q0 <- quad_h0(d, "t_het"); q02 <- quad_h0(d, "t_het", hl = 0.1, nmu = 601, he = 0.125)
    df <- data.frame(jeu = j, loi = c(rep("normal_het", 5), "t_het"), forme = c(FORMES, "H0"),
                     grossier = c(qn, q0), fin = c(qn2, q02))
    ecrit_atomique(df, fn); print(df, digits = 9, row.names = FALSE); flush.console()
  }
}
if (mode == "naifn") {
  j <- as.integer(args[2]); fn <- file.path(DIR_CAS, sprintf("naifn_j%d.csv", j))
  if (!file.exists(fn)) {
    d <- donnees(j); set.seed(4600 + j); r <- naif_normal(d, as.integer(args[3]))
    r$jeu <- j; ecrit_atomique(r, fn); print(r, digits = 6, row.names = FALSE)
  }
}
if (mode == "naift") {
  j <- as.integer(args[2]); fn <- file.path(DIR_CAS, sprintf("naift_j%d.csv", j))
  if (!file.exists(fn)) {
    d <- donnees(j); set.seed(4700 + j); r <- naif_t(d, as.integer(args[3]))
    r$jeu <- j; ecrit_atomique(r, fn); print(r, digits = 6, row.names = FALSE)
  }
}
if (mode == "table") ecrit_table()

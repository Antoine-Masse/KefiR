# Agent 24 : V1 calcul independant de hetero, V2 zs a F grand.
# Derivation propre : xbar_i | mu, a, s ~ N(mu + a_i, D_i), D_i = s_i^2/n_i.
# a ~ N(0, c (I - J/k)), c = g tau^2 ; mu plat. Alors l'integrale sur mu et a vaut
# (2 pi)^{-(k-1)/2} [W prod(D_i + c)]^{-1/2} exp(-(V - U^2/W)/2),
# w_i = 1/(D_i+c), W = sum w, U = sum w x, V = sum w x^2 ; H0 : c = 0.
# Facteur intra : (2 pi s^2)^{-(n-1)/2} n^{-1/2} exp(-SS/(2 s^2)). Mesure : d log s_i.
# NB p(s^2) prop 1/s^2 <=> p(log s) prop 1, a une constante 2 pres commune a H0 et H1 ;
# on prend d(log s^2) = 2 d(log s), soit 2^k en facteur.
source("agent_24_commun.R")
lse <- function(l) { m <- max(l); m + log(sum(exp(l - m))) }
stats_of <- function(y, g) {
  g <- factor(g); x <- as.numeric(tapply(y, g, mean)); n <- as.numeric(table(g))
  ss <- as.numeric(tapply(y, g, function(v) sum((v - mean(v))^2)))
  list(k = length(n), n = n, x = x, ss = ss, N = sum(n))
}
grille_t <- function(st, i, hm = 1, bornes = NULL) {
  n <- st$n[i]; ss <- st$ss[i]; sdl <- 1/sqrt(2*(n-1)); h <- min(.25, .5*sdl)*hm
  t0 <- .5*log(ss/(n-1))
  dl <- uniroot(function(d) (n-1)/2*(exp(2*d)-1-2*d) - 45, c(0, 10))$root + .3
  dist <- log(max(abs(st$x - mean(st$x)))*sqrt(n) + 1e-300)
  tr <- max(t0, dist) + 32/(n-1) + .5
  if (!is.null(bornes)) { # a priori tronque : trapeze avec demi-poids aux bornes
    m <- max(5, ceiling((bornes[2]-bornes[1])/h)); t <- seq(bornes[1], bornes[2], length.out = m+1)
    w <- rep(diff(t)[1], m+1); w[c(1,m+1)] <- w[1]/2; return(list(t=t, w=w)) }
  t <- seq(t0 - dl, tr, by = h); list(t = t, w = rep(h, length(t)))
}
# Route A : grille deterministe sur (log s_i, log g), reduction analytique de mu et a.
route_A <- function(y, g, r = .5, hm = 1, bornes = NULL) {
  st <- stats_of(y, g); k <- st$k
  G <- lapply(1:k, function(i) grille_t(st, i, hm, bornes))
  dims <- vapply(G, function(z) length(z$t), 1L)
  idx <- as.matrix(expand.grid(lapply(dims, seq_len)))
  s2 <- vapply(1:k, function(i) exp(2*G[[i]]$t[idx[,i]]), numeric(nrow(idx)))
  lw <- Reduce(`+`, lapply(1:k, function(i) log(G[[i]]$w[idx[,i]])))
  D <- sweep(s2, 2, st$n, "/"); tau2 <- rowMeans(s2)
  li <- lw + k*log(2) - (k-1)/2*log(2*pi)
  for (i in 1:k) li <- li - (st$n[i]-1)/2*log(2*pi*s2[,i]) - .5*log(st$n[i]) - st$ss[i]/(2*s2[,i])
  kern <- function(c) { w <- 1/(D + c); W <- rowSums(w); U <- rowSums(sweep(w,2,st$x,"*"))
    V <- rowSums(sweep(w,2,st$x^2,"*"))
    -.5*(log(W) - rowSums(log(w))) - .5*(V - U^2/W) }
  lh0 <- lse(li + kern(0))
  # Elagage (k>=3 seulement) : garde les points a moins de 45 unites log du maximum, a 7 valeurs de g balayees
  if (k >= 3 && is.null(bornes)) { gs <- exp(seq(-4, max(log(max((st$x-mean(st$x))^2)/min(st$ss/st$n/(st$n-1)))+log(st$N), 3) + 8, length.out = 7))
    best <- li + kern(0); for (gi in gs) best <- pmax(best, li + kern(gi*tau2) - .5*(k-1)*0)
    keep <- best > max(best) - 45
    s2 <- s2[keep,,drop=FALSE]; D <- D[keep,,drop=FALSE]; tau2 <- tau2[keep]; li <- li[keep] }
  # grille en log g
  ming <- r^2/2/90; dd <- max((st$x-mean(st$x))^2)
  gp <- dd/min(st$ss/(st$n-1)/st$n, (st$ss/(st$n-1)))*st$N + r^2
  gmax <- max(log(gp), 3) + 80/k
  hg <- .4*hm; tg <- seq(log(ming), gmax, by = hg)
  lg <- numeric(length(tg))
  for (j in seq_along(tg)) { gv <- exp(tg[j])
    lpg <- log(r^2/2)/2 - lgamma(.5) - 1.5*tg[j] - r^2/(2*gv) + tg[j] # densite IG * jacobien
    lg[j] <- lse(li + kern(gv*tau2)) + lpg + log(hg) }
  list(lh1 = lse(lg), lh0 = lh0, lbf = lse(lg) - lh0, npts = nrow(s2)*length(tg))
}
# Route B (k=2) : sans g ni matrice : a priori Cauchy(0, r tau) sur z = (a1-a2)/sqrt2, grille z = r tau sinh(u)
route_B <- function(y, g, r = .5, hm = 1) {
  st <- stats_of(y, g); G <- lapply(1:2, function(i) grille_t(st, i, hm))
  dx <- st$x[1]-st$x[2]; hu <- .02*hm; u <- seq(-40, 40, by = hu); lf <- numeric(length(G[[1]]$t))
  base2 <- -(st$n[2]-1)*G[[2]]$t - st$ss[2]/(2*exp(2*G[[2]]$t)) - (st$n[2]-1)/2*log(2*pi) - .5*log(st$n[2])
  su <- sinh(u); lcu <- -log(pi) - log(cosh(u))
  for (a in seq_along(lf)) { s1 <- exp(2*G[[1]]$t[a])
    s2 <- exp(2*G[[2]]$t); tau <- sqrt((s1+s2)/2); v <- s1/st$n[1]+s2/st$n[2]
    Zm <- outer(r*tau, su) # z
    lz <- -.5*log(2*pi*v) - (matrix(dx, length(s2), length(u)) - sqrt(2)*Zm)^2/(2*v)
    # masse : Cauchy en z, dz = r tau cosh(u) du => (1/pi)/cosh(u) du, independant de tau
    m <- apply(sweep(lz, 2, lcu, "+"), 1, lse) + log(hu)
    lf[a] <- lse(m + base2 + log(G[[2]]$w)) - (st$n[1]-1)*G[[1]]$t[a] - st$ss[1]/(2*s1) - (st$n[1]-1)/2*log(2*pi) - .5*log(st$n[1]) + log(G[[1]]$w[a]) }
  lh1 <- lse(lf) + 2*log(2)
  list(lh1 = lh1)
}
# Controle direct de la reduction analytique, k=3 : integrale numerique sur (mu, z1, z2) a (s, g) fixes
reduc_k3 <- function(x, n, s2, c) {
  D <- s2/n; C <- cbind(c(1,-1,0)/sqrt(2), c(1,1,-2)/sqrt(6)); sdz <- sqrt(c)
  sdm <- sqrt(max(D)+c); nod <- seq(-9, 9, by = .12)
  gr <- as.matrix(expand.grid(mu = mean(x)+nod*sdm, z1 = nod*min(sdz, 3*sqrt(max(D))), z2 = nod*min(sdz, 3*sqrt(max(D)))))
  # centrage grossier : z autour de 0 -> on elargit si donnees loin
  vol <- (.12*sdm)*(.12*min(sdz,3*sqrt(max(D))))^2
  zp <- gr[,2:3]; a <- zp %*% t(C)
  l <- -rowSums(zp^2)/(2*c) - log(2*pi*c)
  for (i in 1:3) l <- l - .5*log(2*pi*D[i]) - (x[i]-gr[,1]-a[,i])^2/(2*D[i])
  num <- lse(l) + log(vol)
  w <- 1/(D+c); W <- sum(w); U <- sum(w*x); V <- sum(w*x^2)
  ana <- -log(2*pi) - .5*(log(W) - sum(log(w))) - .5*(V - U^2/W)
  c(num = num, ana = ana)
}
faire <- function(n, sd, mu, seed) { set.seed(seed); list(y = rnorm(sum(n), rep(mu, n), rep(sd, n)), g = factor(rep(seq_along(n), n)), n=n, sd=sd, mu=mu, seed=seed) }
JEUX <- list(
  petit_3_3      = faire(c(3,3), c(1,1), c(0,1), 2411),
  deseq_5_30     = faire(c(5,30), c(1,1), c(0,.6), 2412),
  var_diff_10_10 = faire(c(10,10), c(1,10), c(0,3), 2413),
  effet_nul      = faire(c(10,10), c(1,1), c(0,0), 2414),
  effet_fort     = faire(c(10,10), c(1,1), c(0,5), 2415),
  deseq_var_3_20 = faire(c(3,20), c(5,1), c(0,2), 2416),
  k3_5_5_40_A    = faire(c(5,5,40), c(4,4,1), c(0,0,0), 2417),
  k3_5_5_40_B    = faire(c(5,5,40), c(1,1,4), c(0,1,3), 2418),
  k3_4_4_4_fort  = faire(c(4,4,4), c(1,1,1), c(0,0,6), 2419))
OUT <- list()
tab <- list()
for (nm in names(JEUX)) { J <- JEUX[[nm]]; t0 <- Sys.time()
  A1 <- route_A(J$y, J$g, .5, 1); A2 <- route_A(J$y, J$g, .5, .7)
  F <- HET(J$y, J$g, .5, 20000L, 1L)
  B <- if (length(J$n) == 2) route_B(J$y, J$g, .5, 1)$lh1 else NA
  st <- stats_of(J$y, J$g); tauh <- sqrt(mean(st$ss/(st$n-1)))
  tab[[nm]] <- data.frame(jeu = nm, lBF_grille = A2$lbf, err_quad_pct = 100*abs(A1$lbf-A2$lbf),
    lBF_fn = F$log_bf, err_pct_fn = F$error_pct, ecart_pct = 100*(F$log_bf-A2$lbf),
    z = 100*(F$log_bf-A2$lbf)/F$error_pct, lh1_grille = A2$lh1, lh0_grille = A2$lh0,
    lh1_fn = F$log_ml_h1, lh0_fn = F$log_ml_h0,
    lh1_routeB = B, sec = as.numeric(difftime(Sys.time(), t0, units="secs")), ess = F$ess, npts=A2$npts,
    d_h1_brut = F$log_ml_h1 - A2$lh1, d_h0_brut = F$log_ml_h0 - A2$lh0, N1logtau = (st$N-1)*log(tauh))
  print(tab[[nm]][,c(1:7,12,13)], digits=6); saveRDS(list(tab=tab), paste0(SORTIES,"agent_24_v1.rds")); flush.console() }
OUT$tab <- do.call(rbind, tab)
# Reduction analytique k=3 (3 points)
OUT$reduc <- rbind(reduc_k3(c(0,1,3), c(5,5,40), c(1,2,.5), .7), reduc_k3(c(0,0,0), c(4,4,4), c(1,1,1), .3), reduc_k3(c(0,-1,.5), c(5,5,20), c(9,1,4), 2))
print(OUT$reduc)
# Couverture : 30 graines
COUV <- list()
for (nm in c("petit_3_3","k3_5_5_40_A","var_diff_10_10")) { J <- JEUX[[nm]]; ref <- OUT$tab$lBF_grille[OUT$tab$jeu==nm]
  for (nd in c(4000L, 20000L)) { if (nm!="k3_5_5_40_A" && nd==20000L) next
    v <- t(vapply(1:30, function(s) { r <- HET(J$y,J$g,.5,nd,s); c(r$log_bf, r$error_pct, r$ess) }, numeric(3)))
    COUV[[paste(nm,nd)]] <- data.frame(jeu=nm, n_draws=nd, sd_logBF_pct=100*sd(v[,1]), err_pct_moy=mean(v[,2]), rapport=100*sd(v[,1])/mean(v[,2]),
      biais_pct=100*(mean(v[,1])-ref), se_biais=100*sd(v[,1])/sqrt(30), part_dans_1err=mean(abs(100*(v[,1]-ref))<=v[,2]), part_dans_2err=mean(abs(100*(v[,1]-ref))<=2*v[,2]), ess_min=min(v[,3]))
    print(COUV[[paste(nm,nd)]], digits=4) } }
OUT$couv <- do.call(rbind, COUV)
OUT$jeux <- lapply(JEUX, function(J) J[c("n","sd","mu","seed")]); saveRDS(OUT, paste0(SORTIES,"agent_24_v1.rds"))
cat("V1 fini", format(Sys.time()), "\n")

# ---- V2 ----
# Marginale JZS propre : contrastes orthonormes C, A = C'(diag n - n n'/N)C = V L V'.
# BF(g) = prod(1+g l)^{-1/2} [(SSE + sum c^2/(1+g l))/SST]^{-(N-1)/2}
zs_comp <- function(n, F) { # jeu de F exact, residus fixes
  k <- length(n); set.seed(2420); g <- factor(rep(1:k, n)); e <- rnorm(sum(n)); e <- e - ave(e, g)
  m <- seq(-1, 1, length.out = k); mg <- rep(m, n)
  dfb <- k-1; dfe <- sum(n)-k; sb <- sum((mg-mean(mg))^2); se <- sum(e^2)
  lam <- sqrt(F*(dfb/dfe)*se/sb); y <- lam*mg + e
  list(y=y, g=g, n=n, k=k, sse=sum((y-ave(y,g))^2), group_mean=as.numeric(tapply(y,g,mean)), lam=lam) }
log_bf_int <- function(Z, r = .5) { k <- Z$k; n <- Z$n; N <- sum(n)
  C <- qr.Q(qr(cbind(1, diag(k)[, -k, drop=FALSE])))[, -1, drop=FALSE]
  A <- t(C) %*% (diag(n) - tcrossprod(n)/N) %*% C; e <- eigen(A, symmetric=TRUE); l <- e$values
  xs <- n*(Z$group_mean - sum(n*Z$group_mean)/N) # X'y en coordonnees groupes
  cc <- as.numeric((t(e$vectors) %*% t(C) %*% xs))/sqrt(l); sse <- Z$sse; ssb <- sum(cc^2)
  lfn <- function(t) { g <- exp(t); vapply(g, function(gi) sum(-.5*log1p(gi*l)) - (N-1)/2*(log(sse+sum(cc^2/(1+gi*l))) - log(sse+ssb)), 0) + (.5*log(r^2/2) - lgamma(.5) - .5*t - r^2/(2*g)) }
  o <- optimize(lfn, c(-20, 60 + log(1+ssb/sse)*1.2), maximum=TRUE, tol=1e-6); pk <- o$objective
  # fenetre autour du pic et integrate sur echelle centree (fonction exponentiellement bornee)
  ii <- integrate(function(t) exp(lfn(t)-pk), o$maximum-40, o$maximum+80/k, rel.tol=1e-10, subdivisions=2000L, stop.on.error=FALSE)
  list(lbf = pk + log(ii$value), err = ii$abs.error/ii$value, t_pic = o$maximum, lfn=lfn, pk=pk) }
mp <- requireNamespace("Rmpfr", quietly=TRUE)
log_bf_mpfr <- function(Z, I, r=.5) { # meme integrande, trapezes, tout en mpfr a 300 bits, equilibre uniquement
  if (!mp || diff(range(Z$n))!=0) return(NA)
  suppressPackageStartupMessages(library(Rmpfr)); k<-Z$k; n<-Z$n[1]; N<-k*n; dfb<-k-1
  ssb <- sum(n*(Z$group_mean-mean(Z$group_mean))^2); P<-300
  sse<-mpfr(Z$sse,P); ssb<-mpfr(ssb,P); tt<-seq(I$t_pic-30, I$t_pic+60, by=.05); s<-mpfr(0,P); rr<-mpfr(r,P)
  gv <- exp(mpfr(tt,P)); li <- -(k-1)/2*log(1+gv*n) - (N-1)/2*(log(sse+ssb/(1+gv*n)) - log(sse+ssb)) + .5*log(rr^2/2) - lgamma(mpfr(.5,P)) - .5*mpfr(tt,P) - rr^2/(2*gv)
  as.numeric(log(sum(exp(li))*mpfr(.05,P))) }
library(BayesFactor)
V2 <- list()
for (cas in list(list(n=c(10,10,10)), list(n=c(5,5,40)), list(n=c(20,20)))) for (Fv in c(1e1,1e3,1e4,1e5,1e6,1e8,1e12,1e16,1e20,1e24,1e28)) {
  Z <- zs_comp(cas$n, Fv); fg <- if (is.finite(Fv)) NULL
  fn <- tryCatch(as.numeric(ZS(Z$y, Z$g, .5)), error=function(e) NA)
  I <- tryCatch(log_bf_int(Z), error=function(e) list(lbf=NA, err=NA))
  bf <- if (Fv<=1e6) tryCatch(as.numeric(BayesFactor::extractBF(anovaBF(y~g, data.frame(y=Z$y,g=Z$g), rscaleFixed=.5, progress=FALSE))[,1:2]), error=function(e) c(NA,NA)) else c(NA,NA)
  m <- tryCatch(log_bf_mpfr(Z, I), error=function(e) NA)
  obs <- anova(lm(Z$y~Z$g))[1,4]
  V2[[length(V2)+1]] <- data.frame(n=paste(cas$n,collapse=","), F_vise=Fv, F_obs=obs, zs=fn, integrate=I$lbf, err_int=I$err, mpfr=m,
     anovaBF=if(is.na(bf[1])) NA else log(bf[1]), err_anovaBF=bf[2], ecart_zs_int=fn-I$lbf, ecart_zs_BF=if(is.na(bf[1])) NA else fn-log(bf[1]), ecart_zs_mpfr=fn-m)
  print(V2[[length(V2)]], digits=10); flush.console() }
OUT$V2 <- do.call(rbind, V2); saveRDS(OUT, paste0(SORTIES,"agent_24_v1.rds"))
cat("V2 fini", format(Sys.time()), "\n")

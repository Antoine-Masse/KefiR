# Estimateur des lois normales : grille en g, pilote, troncature, erreur annoncee.
# Reference k=2 et k=3 : quadrature deterministe sur la proportion de precision, avec
# integrale en g sur une grille large ecrite ici.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
ns <- asNamespace("KefiR")
orig_gint <- get(".bay_law_gint", ns)
lse <- function(l) { t <- max(l); t + log(sum(exp(l - t))) }
# g sur grille tres large, sans pilote, pas h/hf, sans passer par le code teste
gint_wide <- function(st, cc, rs, lo=-60, hi=130, h=0.02) {
  tg <- seq(lo, hi, by=h)
  L <- vapply(tg, function(t) .bay_law_kernel(st, cc, exp(t)), numeric(nrow(cc)))
  if (nrow(cc) == 1L) L <- matrix(L, nrow=1L)
  res <- matrix(NA_real_, nrow(cc), length(rs))
  for (i in seq_along(rs)) {
    r <- rs[i]
    lp <- 0.5*log(r^2/2) - 0.5*log(pi) - 0.5*tg - r^2/(2*exp(tg))
    A <- sweep(L, 2L, lp, "+")
    top <- apply(A, 1L, max)
    res[, i] <- top + log(rowSums(exp(A - top))) + log(h)
  }
  res
}
# reference : quadrature sur le simplexe via alr, trapezes, k = 2 ou 3
ref_het <- function(st, u=1.5, r=.5, lim=60, hz=.05) {
  k <- st$k; const <- .bay_law_const(st)
  z <- seq(-lim, lim, by=hz)
  grid <- if (k==2) matrix(z,ncol=1) else as.matrix(expand.grid(z,z))
  th <- cbind(exp(grid),1)/(1+rowSums(exp(grid)))
  th <- pmax(th,1e-300)
  lprior <- lgamma(k*u)-k*lgamma(u)+u*rowSums(log(th))
  cc <- 1/(k*th)
  h0 <- const + lse(.bay_law_kernel(st,cc,0)+lprior) + (k-1)*log(hz)
  # H1 : grille large par blocs
  b <- split(seq_len(nrow(cc)), ceiling(seq_len(nrow(cc))/2000))
  v <- unlist(lapply(b,function(i) gint_wide(st,cc[i,,drop=FALSE],r,h=.1)[,1]))
  h1 <- const + lse(v+lprior) + (k-1)*log(hz)
  c(h0=h0,h1=h1)
}
mk <- function(n,m,s,seed=1) { set.seed(seed); y<-unlist(Map(function(a,b,c) { v<-rnorm(a); if(a>1) v<-(v-mean(v))/sd(v); b+c*v },n,m,s)); list(y=y,g=rep(seq_along(n),n)) }

cat("=== A. Reference contre IS, normal_het (S=4000, 3 graines) : ecart et erreur annoncee\n")
cas <- list(
 "(2,30) d=1e4"=list(c(2,30),c(1e4,0),c(1,1)),
 "(1,30) d=100"=list(c(1,30),c(100,0),c(1,1)),
 "(2,2) d=1e6"=list(c(2,2),c(1e6,0),c(1,1)),
 "(5,5) sd 1e-3,1 d=.5"=list(c(5,5),c(.5,0),c(1e-3,1)),
 "(4,4) sd 1e-6,1 d=.5"=list(c(4,4),c(.5,0),c(1e-6,1)),
 "(30,30) sd 1,1 d=.5"=list(c(30,30),c(.5,0),c(1,1)),
 "(5,5,40) sd4,4,1"=list(c(5,5,40),c(0,0,0),c(4,4,1)),
 "(2,2,40) d=1e3"=list(c(2,2,40),c(1e3,-1e3,0),c(1,1,1)),
 "(1,1,30) d=50"=list(c(1,1,30),c(50,-50,0),c(1,1,1)),
 "(3,3,3)"=list(c(3,3,3),c(0,1,2),c(1,1,1)),
 "(10,10,10) sd1,1,1e4 d=1"=list(c(10,10,10),c(0,1,0),c(1,1,1e4))
)
for (nm in names(cas)) {
  d <- do.call(mk,cas[[nm]]); st <- .bay_laws_stats(d$y,d$g)
  rf <- ref_het(st)
  cat(sprintf("%-26s ref H0 %.3f H1 %.3f logBF %.3f\n",nm,rf[1],rf[2],rf[2]-rf[1]))
  for (sd in 1:3) {
    tb <- tab_of(d$y,d$g,draws=4000L,seed=sd,laws="normal_het")
    cat(sprintf("   seed %d  H0 ecart %+.3f (annonce %.3f, z=%+.1f, ess %.0f) H1 ecart %+.3f (annonce %.3f, z=%+.1f, ess %.0f)\n",sd,
      tb$log_ml[1]-rf[1],tb$rel_se[1],(tb$log_ml[1]-rf[1])/tb$rel_se[1],tb$ess[1],
      tb$log_ml[2]-rf[2],tb$rel_se[2],(tb$log_ml[2]-rf[2])/tb$rel_se[2],tb$ess[2]))
  }
  flush.console()
}

cat("\n=== B. Troncature droite : deux groupes a sd infime, un troisieme disperse (normal_het H1)\n")
for (ratio in c(1e2,1e4,1e6,1e8,1e10)) {
  d <- mk(c(10,10,10),c(0,1,.5),c(1/ratio,1/ratio,100)); st <- .bay_laws_stats(d$y,d$g)
  assign(".bay_law_gint",orig_gint,ns)
  a <- tab_of(d$y,d$g,draws=2000L,laws="normal_het")
  assign(".bay_law_gint",function(st,cc,rs,h=NULL) gint_wide(st,cc,rs,h=.05),ns)
  b <- tab_of(d$y,d$g,draws=2000L,laws="normal_het")
  assign(".bay_law_gint",orig_gint,ns)
  cat(sprintf("d/s=%g : code logBF %.3f (H1 %.3f err %.3f) ; grille large logBF %.3f (H1 %.3f) ; diff H1 %.3f\n",ratio,a$log_ml[2]-a$log_ml[1],a$log_ml[2],a$rel_se[2],b$log_ml[2]-b$log_ml[1],b$log_ml[2],b$log_ml[2]-a$log_ml[2]))
  flush.console()
}

cat("\n=== C. Pilote et pas : memes tirages, grille du code contre grille large (normal_het H1)\n")
plans <- list(
 "k=8 (10x7,3) d=1e3"=list(c(rep(10,7),3),c(rep(0,7),1e3),rep(1,8)),
 "k=20 n=3, d=30"=list(rep(3,20),c(rep(0,19),30),rep(1,20)),
 "k=30 n=2, sd 1"=list(rep(2,30),seq(0,1,length.out=30),rep(1,30)),
 "k=20 n=10 sd varies"=list(rep(10,20),seq(0,3,length.out=20),exp(seq(-3,3,length.out=20))),
 "k=50 n=5 d=1"=list(rep(5,50),c(rep(0,49),5),rep(1,50)),
 "k=3 N=3e5, d=.01"=list(rep(1e5,3),c(0,.01,0),c(1,1,1)),
 "k=2 N=2e6, d=.01"=list(rep(1e6,2),c(0,.01),c(1,2)),
 "k=5 n=1e4, d=3 sd var"=list(rep(1e4,5),c(0,0,0,3,0),c(1,3,.1,1,1))
)
for (nm in names(plans)) {
  d <- do.call(mk,plans[[nm]])
  t0 <- Sys.time(); a <- tab_of(d$y,d$g,draws=2000L,laws="normal_het"); t1 <- Sys.time()
  assign(".bay_law_gint",function(st,cc,rs,h=NULL) gint_wide(st,cc,rs,lo=log(min(rs)^2/(1000*st$k)),hi=60,h=min(.25,.5/sqrt(st$k))/4),ns)
  b <- tab_of(d$y,d$g,draws=2000L,laws="normal_het")
  assign(".bay_law_gint",orig_gint,ns)
  cat(sprintf("%-26s H1 code %.4f (ess %.0f err %.3f, %.1fs) large %.4f diff %+.2e\n",nm,a$log_ml[2],a$ess[2],a$rel_se[2],as.numeric(difftime(t1,t0,units="secs")),b$log_ml[2],b$log_ml[2]-a$log_ml[2]))
  flush.console()
}

cat("\n=== D. Invariances sur normal_het, memes graines\n")
d <- mk(c(5,8,12),c(0,1,.3),c(1,2,.5))
base <- tab_of(d$y,d$g,draws=20000L,laws="normal_het")
cat("base", round(base$log_ml,4), "err", round(base$rel_se,4), "\n")
rf <- tab_of(-d$y,d$g,draws=20000L,laws="normal_het"); cat("y -> -y", round(rf$log_ml,4), "\n")
pm <- tab_of(d$y,c(3,1,2)[d$g],draws=20000L,laws="normal_het"); cat("etiquettes", round(pm$log_ml,4), "\n")
ord <- sample(length(d$y)); oo <- tab_of(d$y[ord],d$g[ord],draws=20000L,laws="normal_het"); cat("ordre lignes", round(oo$log_ml,4), "\n")

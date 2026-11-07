# Lois t : posterieur multimodal (contamination pres de 1/3). Reference : quadrature 3D
# sur (mu, log s^2, log(nu - 2)), sans fonction du paquet, loi t_eq, H0.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
ref_t_eq_h0 <- function(y, nm=4000L, ns2=90L, nv=30L) {
  rng <- range(y); pad <- diff(rng)*.3
  N <- length(y)
  ls2 <- seq(log(1e-4), log(diff(rng)^2*4), length.out=ns2)
  lv <- seq(-8, 3, length.out=nv)
  mug <- seq(rng[1]-pad, rng[2]+pad, length.out=nm)
  hm <- mug[2]-mug[1]; hs <- ls2[2]-ls2[1]; hv <- lv[2]-lv[1]
  out <- matrix(-Inf, ns2, nv); peak <- array(NA_real_, c(ns2,nv))
  D2 <- outer(mug, y, "-")^2
  for (j in seq_len(nv)) { nu <- 2+exp(lv[j])
    norm_t <- lgamma((nu+1)/2)-lgamma(nu/2)-.5*log(nu*pi)
    for (i in seq_len(ns2)) {
      ll <- N*(norm_t-.5*ls2[i]) - (nu+1)/2*rowSums(log1p(D2*exp(-ls2[i])/nu))
      t0 <- max(ll); out[i,j] <- t0+log(sum(exp(ll-t0)))+log(hm)
      peak[i,j] <- mug[which.max(ll)]
    }
    out[,j] <- out[,j] - exp(lv[j]) + lv[j]
  }
  t0 <- max(out); lm <- t0+log(sum(exp(out-t0)))+log(hs)+log(hv)
  # masse par region : log s^2 sous ou sur la mediane des modes
  marg_s <- apply(out,1L,function(v){m<-max(v);m+log(sum(exp(v-m)))})
  list(lm=lm, marg_s=marg_s, ls2=ls2)
}
# Donnees sur l'echelle standardisee de .bay_laws_stats ; le lambda lt=log taubar a un jacobien 1
# par rapport a log s^2 a nu fixe, le signe s'integre sans facteur.
cas <- list()
set.seed(1); nA <- 20; for (dd in c(30,100)) for (f in c(.25,1/3)) {
  m <- round(30*f)
  y1 <- c(rnorm(30-m),rnorm(m,dd)); y2 <- c(rnorm(30-m),rnorm(m,dd))
  cas[[sprintf("d=%g f=%.2f (N=60, 2 groupes)",dd,f)]] <- list(y=c(y1,y2),g=rep(1:2,each=30))
}
for (nm in names(cas)) {
  d <- cas[[nm]]; st <- .bay_laws_stats(d$y,d$g)
  ys <- unlist(st$y,use.names=FALSE)
  rf <- ref_t_eq_h0(ys, nm=800L, ns2=60L, nv=20L)
  ms <- rf$marg_s - max(rf$marg_s); pk <- which(diff(sign(diff(ms)))==-2)+1L
  cat(sprintf("%s\n  ref t_eq H0 %.3f ; modes en log s^2 : %s (log-hauteur rel. %s)\n",nm,rf$lm,paste(round(rf$ls2[pk],1),collapse=","),paste(round(ms[pk],1),collapse=",")))
  for (sd in 1:2) { tb <- tab_of(d$y,d$g,draws=3000L,seed=sd,laws="t_eq")
    cat(sprintf("   seed %d H0 %.3f ecart %+.3f annonce %.3f ess %.0f | H1 %.3f annonce %.3f ess %.0f | logBF %.3f\n",sd,tb$log_ml[1],tb$log_ml[1]-rf$lm,tb$rel_se[1],tb$ess[1],tb$log_ml[2],tb$rel_se[2],tb$ess[2],tb$log_ml[2]-tb$log_ml[1]))
  }
  flush.console()
}

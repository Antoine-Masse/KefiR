# Ecart silencieux de .bay_law_kernel_grid sur des lignes extremes, et correction candidate.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
y2 <- c(rep(5,10), 1,2,2,3,3,3,3,4,4,3); g2 <- rep(1:2,each=10); st <- .bay_laws_stats(y2,g2)
gs <- exp(seq(-6,12,by=1))
# correction candidate : centrer M par ligne sur le groupe de plus grand poids
kfix <- function(st,cc,gs) {
  n<-st$n;k<-st$k;S<-nrow(cc);lc<-log(cc);w0<-1/cc
  M <- matrix(rep(st$m,each=S),S,k); M <- M - st$m[max.col(w0,ties.method="first")]
  a<-1/outer(1/n,gs,"+"); W<-w0%*%a; Wm<-(w0*M)%*%a; S2<-(w0*M^2)%*%a
  Q<-drop(w0%*%st$ss)+S2-Wm^2/W
  base<- -0.5*drop(lc%*%(n-1))-0.5*rowSums(lc)
  sweep(-0.5*log(W)-(st$N-1)/2*log(Q),2L,0.5*colSums(log(a)),"+")+base
}
cat("c2 ; ecart max |grille-ancien| ; NaN grille ; ecart max |corrige-ancien| ; NaN corrige ; ecart max vers le haut grille\n")
for (c2 in 10^c(4,8,10,12,14,15,16,17,18,20)) {
  cc <- matrix(c(.5,c2),1)
  A <- .bay_law_kernel_grid(st,cc,gs); F <- kfix(st,cc,gs); B <- vapply(gs,function(g) .bay_law_kernel(st,cc,g),0)
  cat(sprintf("%.0e  %.3g  %d  %.3g  %d  %.3g\n",c2,max(abs(A-B),na.rm=TRUE),sum(is.nan(A)),max(abs(F-B),na.rm=TRUE),sum(is.nan(F)),max(A-B,na.rm=TRUE)))
}
cat("Plan d/s=1e8, ligne c=(1e31,1/3,1e35)\n")
set.seed(1); mk <- function(n,m,s){y<-unlist(Map(function(a,b,c){v<-rnorm(a);v<-(v-mean(v))/sd(v);b+c*v},n,m,s));list(y=y,g=rep(seq_along(n),n))}
d <- mk(c(10,10,10),c(0,1,.5),c(1e-8,1e-8,100)); st3 <- .bay_laws_stats(d$y,d$g)
cc <- matrix(c(1.14e31,1/3,8.8e35),1); A <- .bay_law_kernel_grid(st3,cc,gs); F <- kfix(st3,cc,gs); B <- vapply(gs,function(g) .bay_law_kernel(st3,cc,g),0)
cat(max(abs(A-B),na.rm=TRUE),sum(is.nan(A)),max(abs(F-B),na.rm=TRUE),sum(is.nan(F)),"\n")

# .bay_law_kernel_grid (nouveau) contre .bay_law_kernel (somme de termes positifs).
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
mk <- function(n,m,s,seed=1) { set.seed(seed); y<-unlist(Map(function(a,b,c) { v<-rnorm(a); if(a>1) v<-(v-mean(v))/sd(v); b+c*v },n,m,s)); list(y=y,g=rep(seq_along(n),n)) }
cmp <- function(lab, st, cc, gs=exp(seq(-8,30,by=2))) {
  A <- .bay_law_kernel_grid(st,cc,gs)
  B <- vapply(gs,function(g) .bay_law_kernel(st,cc,g),numeric(nrow(cc))); if(nrow(cc)==1L) B <- matrix(B,nrow=1L)
  d <- A-B
  cat(sprintf("%-40s NaN grille %d / %d ; NaN ancien %d ; max|ecart| %.3g ; ecart > 1 : %d ; ecart max vers le haut %.3g\n",lab,sum(is.nan(A)),length(A),sum(is.nan(B)),max(abs(d),na.rm=TRUE),sum(abs(d)>1,na.rm=TRUE),max(d,na.rm=TRUE)))
}
th_rand <- function(k,S,a=1.5,z=NULL) { x <- matrix(rgamma(k*S,a),S,k); t <- x/rowSums(x); t }
set.seed(1)
cat("--- ordinaire\n")
d <- mk(c(10,14,12),c(0,1,.3),c(1,2,.5)); st <- .bay_laws_stats(d$y,d$g); cmp("3 groupes ordinaires",st,1/(3*th_rand(3,500)))
cat("--- th extremes (alr jusqu'a +- 40, +- 300)\n")
for (lim in c(10,40,300)) { z <- matrix(runif(2*500,-lim,lim),500,2); th <- .bay_alr_inv(z); cmp(sprintf("3 groupes, alr +-%d",lim),st,1/(3*th)) }
cat("--- groupe a variance nulle, moyenne non nulle, tirages du melange defensif\n")
y2 <- c(rep(5,10), 1,2,2,3,3,3,3,4,4,3); g2 <- rep(1:2,each=10); st2 <- .bay_laws_stats(y2,g2)
pr <- .bay_law_proposal(st2,1.5); th <- .bay_law_draw_def(pr,2000L); cmp("n=(10,10), g1 constant",st2,1/(2*th))
cat("--- rapports d/s croissants, tirages du melange defensif\n")
for (ratio in c(1e2,1e4,1e6,1e7,1e8,1e9,1e10)) {
  d <- mk(c(10,10,10),c(0,1,.5),c(1/ratio,1/ratio,100)); st3 <- .bay_laws_stats(d$y,d$g)
  pr <- .bay_law_proposal(st3,1.5); th <- .bay_law_draw_def(pr,2000L); cmp(sprintf("d/s=%g",ratio),st3,1/(3*th))
}
cat("--- effet sur m.bay : erreurs brutes\n")
tc <- function(lab,expr) { r <- tryCatch(expr,error=function(e) {cat(lab,": ERREUR",conditionMessage(e),"\n");NULL}); if(!is.null(r)) cat(sprintf("%s : P=%.4f err=%.2f%% lois=%s\n",lab,r$prob_effect,r$effects$error_pct,paste(sprintf("%s=%.3f",r$laws$law,r$laws$posterior),collapse=" "))); flush.console() }
q <- function(...) m.bay(...,verbose=FALSE,sensitivity=FALSE)
for (s in 1:3) tc(sprintf("n=(10,10) g1 constant a 5, seed %d",s),q(y2,g2,draws=1000L,seed=s))
y3 <- c(rep(3,10), 1,2,2,3,3,3,3,4,4,5); for (s in 1:3) tc(sprintf("n=(10,10) g1 constant a 3, seed %d",s),q(y3,g2,draws=1000L,seed=s))
tc("y=3x5|1..5 (toutes lois)",q(c(3,3,3,3,3,1,2,3,4,5),rep(1:2,each=5),draws=1000L))
tc("y=1,1,1,1,2,3",q(c(1,1,1,1,2,3),c(1,1,1,1,2,2),draws=300L))
for (ratio in c(1e6,1e8,1e10)) { d <- mk(c(10,10,10),c(0,1,.5),c(1/ratio,1/ratio,100)); tc(sprintf("d/s=%g m.bay",ratio),q(d$y,d$g,draws=500L)) }

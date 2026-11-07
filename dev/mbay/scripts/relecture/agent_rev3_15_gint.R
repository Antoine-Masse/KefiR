# Grille de .bay_law_gint contre grille large et fine, memes lignes : ecart sur la marginale.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
wide <- function(st,cc,rs,hf) { tg <- seq(log(min(rs)^2/(1000*st$k)),70,by=hf); L <- .bay_law_kernel_grid(st,cc,exp(tg)); r<-rs
  lp <- 0.5*log(r^2/2)-0.5*log(pi)-0.5*tg-r^2/(2*exp(tg)); A <- sweep(L,2L,lp,"+"); top <- A[cbind(seq_len(nrow(A)),max.col(A,ties.method="first"))]; top+log(rowSums(exp(A-top)))+log(hf) }
mk <- function(n,m,s,seed=1) { set.seed(seed); y<-unlist(Map(function(a,b,c) { v<-rnorm(a); if(a>1) v<-(v-mean(v))/sd(v); b+c*v },n,m,s)); list(y=y,g=rep(seq_along(n),n)) }
plans <- list(
 "k=8 (10x7,3) d=1e3"=list(c(rep(10,7),3),c(rep(0,7),1e3),rep(1,8)),
 "k=20 n=3, d=30"=list(rep(3,20),c(rep(0,19),30),rep(1,20)),
 "k=30 n=2, sd 1"=list(rep(2,30),seq(0,1,length.out=30),rep(1,30)),
 "k=20 n=10 sd varies"=list(rep(10,20),seq(0,3,length.out=20),exp(seq(-3,3,length.out=20))),
 "k=50 n=5 d=5"=list(rep(5,50),c(rep(0,49),5),rep(1,50)),
 "k=2 N=2e6 d=.01"=list(rep(1e6,2),c(0,.01),c(1,2)),
 "k=5 n=1e4 d=3"=list(rep(1e4,5),c(0,0,0,3,0),c(1,3,.1,1,1)),
 "(2,30) d=1e4"=list(c(2,30),c(1e4,0),c(1,1)),
 "(10,10,10) d/s=1e5"=list(c(10,10,10),c(0,1,.5),c(1e-5,1e-5,100)))
for (nm in names(plans)) { d <- do.call(mk,plans[[nm]]); st <- .bay_laws_stats(d$y,d$g); k<-st$k
  set.seed(2); pr <- .bay_law_proposal(st,1.5); th <- .bay_law_draw_def(pr,1000L); cc <- 1/(k*th)
  a <- .bay_law_gint(st,cc,.5)[,1]; b <- wide(st,cc,.5,min(.25,.5/sqrt(k))/4)
  lw <- b + (1.5-1)*rowSums(log(th)) - .bay_law_lq_def(pr,th); wt <- exp(lw-max(lw))
  cat(sprintf("%-22s ecart marginale %+.2e ; ligne max (poids>1e-6) %.2e ; toutes lignes %.2e\n",nm,log(sum(wt*exp(a-b))/sum(wt)),max(abs(a-b)[wt>1e-6]),max(abs(a-b)))); flush.console()
}

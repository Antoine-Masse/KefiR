source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
set.seed(11)
y <- c(rnorm(5,0,1), rnorm(8,.7,1), rnorm(11,-.3,1.4)); g <- rep(c("a","b","c"),c(5,8,11))
st <- .bay_laws_stats(y,g); k<-3; one<-matrix(1,1,3); const <- .bay_law_const(st); r<-.5
cat("--- A. forme s, c=1, contre ZS a 2 groupes s contre reste (meme r)\n")
h0 <- const+.bay_law_kernel(st,one,0)
for(s in 1:3){ gs <- ifelse(g==st$levels[s],"s","o"); lbf <- .bay_bf_oneway_zs(y,gs,r)
 x <- const+drop(.bay_law_gint(st,one,r*sqrt(2),s=s)); cat(sprintf("s=%d grille %.8f  ZS2 %.8f ecart %.2e\n",s,x-h0,lbf,x-h0-lbf))}
cat("--- A2. forme tous, c=1, contre ZS\n")
x <- const+drop(.bay_law_gint(st,one,r)); cat(sprintf("grille %.8f ZS %.8f ecart %.2e\n",x-h0,.bay_bf_oneway_zs(y,g,r),x-h0-.bay_bf_oneway_zs(y,g,r)))
cat("--- B. k=2 : forme s == forme tous a variances egales (annonce en-tete)\n")
y2<-y[g!="c"]; g2<-g[g!="c"]; st2<-.bay_laws_stats(y2,g2); o2<-matrix(1,1,2)
for(s in 1:2) cat(sprintf("s=%d forme s %.8f, tous %.8f\n",s,.bay_law_gint(st2,o2,r*sqrt(2),s=s),.bay_law_gint(st2,o2,r)))
cat("    k=2 hetero c=(.5,2) : s1 s2 tous (les formes different, le melange n'est pas utilise)\n")
c2<-matrix(c(.5,2),1,2); print(c(.bay_law_gint(st2,c2,r*sqrt(2),s=1),.bay_law_gint(st2,c2,r*sqrt(2),s=2),.bay_law_gint(st2,c2,r)))
cat("--- C. noyau contre formule matricielle independante (c, g fixes)\n")
brute <- function(st,cc,g,s=NULL){ yy<-unlist(st$y,use.names=FALSE); grp<-rep(seq_len(st$k),st$n); N<-st$N
 V<-diag(cc[grp]); for(i in seq_len(st$k)) if(is.null(s)||s==i){ id<-which(grp==i); V[id,id]<-V[id,id]+g*cc[i]}
 Vi<-solve(V); one<-rep(1,N); W<-sum(Vi); mu<-sum(Vi%*%yy)/W; RSS<-drop(t(yy-mu)%*%Vi%*%(yy-mu))
 -(N-1)/2*log(2*pi)-.5*determinant(V)$modulus-.5*log(W)+lgamma((N-1)/2)-(N-1)/2*log(RSS/2)}
cc<-matrix(c(.3,1.7,1),1,3)
for(gg in c(.01,.3,5)) for(s in list(NULL,1,2,3)) {
 a<-const+.bay_law_kernel_grid(st,cc,gg,s); b<-brute(st,cc,gg,s); cat(sprintf("g=%g s=%s kernel %.8f brute %.8f ecart %.1e\n",gg,format(s),a,b,a-b)) }
cat("    H0 : "); a<-const+.bay_law_kernel(st,cc,0); b<-brute(st,cc,0,NULL); cat(a,b,a-b,"\n")
cat("--- D. loi de a_s : dcauchy(r sqrt2 sigma) contre integrale N(0,g sigma^2) x IG(1/2,(r sqrt2)^2/2)\n")
for(a in c(.1,1,5)){sig<-1.3;I<-integrate(function(gg) dnorm(a,0,sqrt(gg)*sig)*(r/sqrt(pi))*gg^(-1.5)*exp(-r^2/gg),0,Inf)$value
 cat(sprintf("a=%g integ %.8f cauchy %.8f\n",a,I,dcauchy(a,0,r*sqrt(2)*sig)))}
cat("--- E. noyau t forme s, nu tres grand, contre dnorm (c=1 puis hetero) ; a_s dans le mu\n")
x<-c(.2,.5,log(1.1),log(1e7)); pa<-.bay_t_parts(matrix(x,1),st,TRUE,FALSE,1.5,s=2)
ll<-sum(dnorm(st$y[[1]],.2,sqrt(1/1.1),log=TRUE))+sum(dnorm(st$y[[2]],.7,sqrt(1/1.1),log=TRUE))+sum(dnorm(st$y[[3]],.2,sqrt(1/1.1),log=TRUE))
cat(sprintf("t base-lp_nu %.5f normale %.5f\n",pa$base-(-exp(x[4])+x[4]),ll))
th<-c(.2,.5,.3); z<-log(th[1:2]/th[3]); x<-c(.2,.5,log(1.1),z,log(1e7)); pa<-.bay_t_parts(matrix(x,1),st,TRUE,TRUE,1.5,s=2)
sg<-sqrt(1/(3*th)/1.1);ll<-sum(dnorm(st$y[[1]],.2,sg[1],log=TRUE))+sum(dnorm(st$y[[2]],.7,sg[2],log=TRUE))+sum(dnorm(st$y[[3]],.2,sg[3],log=TRUE))
ldir<-lgamma(4.5)-3*lgamma(1.5)+1.5*sum(log(th))
cat(sprintf("t het %.5f normale+dir %.5f\n",pa$base-(-exp(x[6])+x[6]),ll+ldir))
cat("--- F. IS t_eq contre grille brute (H0 ; forme s=2) sur les memes donnees, 3 graines\n")
grid_int <- function(st,h1,s,r){ 
  p <- if(!h1) 1 else 2
  fn <- function(x) .bay_t_logpost(x,st,h1,FALSE,1.5,r,s)
  gen <- function(lo,hi,m) {h<-(hi-lo)/(m-1); list(x=seq(lo,hi,by=h),h=h)}
  G<-list(); if(!h1){G$mu<-gen(-2.5,2.5,121)} else {G$mu<-gen(-2.5,2.5,81);G$a<-gen(-5,5,81)}
  G$lt<-gen(-4,5,61);G$lv<-gen(-6,5,31)
  X<-as.matrix(expand.grid(lapply(G,`[[`,"x")))
  out<-numeric(nrow(X)); bl<-split(seq_len(nrow(X)),ceiling(seq_len(nrow(X))/2e5))
  for(b in bl){pa<-.bay_t_parts(X[b,,drop=FALSE],st,h1,FALSE,1.5,s); out[b]<-pa$base+(if(h1).bay_t_prior_theta(pa,r,st$k) else 0)}
  m<-max(out); m+log(sum(exp(out-m)))+sum(log(vapply(G,`[[`,numeric(1),"h")))
}
g0<-grid_int(st,FALSE,NULL,r); g1<-grid_int(st,TRUE,2L,r)
cat(sprintf("grille H0 %.4f  forme s=2 %.4f\n",g0,g1))
for(sd in 1:3){ set.seed(sd); a<-.bay_t_ml(st,FALSE,FALSE,1.5,r,4000L); b<-.bay_t_ml(st,TRUE,FALSE,1.5,r,4000L,s=2L)
 cat(sprintf(" graine %d IS H0 %.4f (ess %.0f) s=2 %.4f (ess %.0f)\n",sd,a$log_ml,a$ess,b$log_ml,b$ess))}
cat("--- G. melange t : log 0,5/k ; recalcul a la main\n")
set.seed(5); m<-.bay_t_h1_mix(st,FALSE,1.5,r,2000L); set.seed(5); z<-.bay_t_ml(st,TRUE,FALSE,1.5,r,2000L); sp<-sapply(1:3,function(s).bay_t_ml(st,TRUE,FALSE,1.5,r,2000L,s=s)$log_ml)
cat(sprintf("mix %.4f ; main %.4f\n",m$log_ml,.bay_log_sum_exp(c(z$log_ml+log(.5),sp+log(.5/3)))))
cat("--- H. .bay_law_normal egal a la main, et continuite k=2 contre k=3 (groupe vide ajoute presque nul)\n")
nm<-.bay_law_normal(st,r,1.5,500L,1L); hand<-.bay_lse2(nm$eq_h1$log_ml_zs,const+.bay_row_lse(rbind(sapply(1:3,function(s) drop(.bay_law_gint(st,one,r*sqrt(2),s=s)))))-log(3))+log(.5)
cat(sprintf("eq_h1 %.6f main %.6f\n",nm$eq_h1$log_ml,hand))
cat("eq_h1 via gint_mix (ligne c=1) ",const+.bay_law_gint_mix(st,one,r)[1],"\n")

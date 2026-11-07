# Controle independant k=3 (agent_40 ne couvre que k=2).
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
gl <- function(n){i<-seq_len(n-1L);b<-i/sqrt(4*i^2-1);J<-matrix(0,n,n);J[cbind(i,i+1L)]<-b;J[cbind(i+1L,i)]<-b;e<-eigen(J,TRUE);o<-order(e$values);list(x=e$values[o],w=2*e$vectors[1L,o]^2)}
nod <- function(lo,hi,n){g<-gl(n);list(x=(hi+lo)/2+(hi-lo)/2*g$x,w=(hi-lo)/2*g$w)}
lse <- function(l){t<-max(l);t+log(sum(exp(l-t)))}
set.seed(11)
y <- c(rnorm(5,0,1),rnorm(6,.8,2),rnorm(4,-.3,.7)); g <- rep(1:3,c(5,6,4))
st <- .bay_laws_stats(y,g)
cat("(1) kernel+const normal, c et g fixes, par integration brute sur a_i, mu, log tau\n")
raw_normal <- function(cc,gg) {
  mus <- nod(-6,6,50); lts <- nod(-5,4,40)
  a <- nod(-8,8,100)
  v <- numeric(0)
  for (im in seq_along(mus$x)) for (it in seq_along(lts$x)) {
    mu<-mus$x[im]; tau<-exp(lts$x[it]); ll<-0
    for (i in 1:3) {
      s2 <- cc[i]/tau
      if (gg==0) { ll <- ll+sum(dnorm(st$y[[i]],mu,sqrt(s2),log=TRUE)) } else {
        lpa <- dnorm(a$x,0,sqrt(gg*s2),log=TRUE)+log(a$w)
        lly <- rowSums(outer(a$x,st$y[[i]],function(z,yy) dnorm(yy,mu+z,sqrt(s2),log=TRUE)))
        ll <- ll+lse(lpa+lly)
      }
    }
    v <- c(v,ll+log(mus$w[im])+log(lts$w[it]))  # p(log tau)=1 (1/tau dtau)
  }
  lse(v)
}
for (cfg in list(list(c(1,1,1),0),list(c(.5,2,1.5),0),list(c(.5,2,1.5),.3),list(c(1,1,1),1.2),list(c(3,.4,.9),.05))) {
  cc <- matrix(cfg[[1]],1,3); gg<-cfg[[2]]
  cat(sprintf("c=(%s) g=%g : formule %.5f  brut %.5f\n",paste(cfg[[1]],collapse=","),gg,.bay_law_const(st)+.bay_law_kernel(st,cc,gg),raw_normal(cfg[[1]],gg)))
  flush.console()
}
cat("(2) a priori de theta sachant sigma, k=3 : formule contre integration numerique sur g et mu\n")
th <- c(.3,-.5,1.1); ds <- c(.5,2,1.2); r<-.7
parts <- list(th_loc=matrix(th,1),lc=matrix(log(ds),1),lt=0)
num <- integrate(function(lg) vapply(lg,function(l){G<-exp(l)
  inner<-integrate(function(m) vapply(m,function(mm)exp(sum(dnorm(th,mm,sqrt(G*ds),log=TRUE))),0),-Inf,Inf,subdivisions=2000L,rel.tol=1e-9)$value
  inner*exp(.5*log(r^2/2)-.5*log(pi)-1.5*l-r^2/(2*G))*G},0),-12,14,subdivisions=500L,rel.tol=1e-8)$value
cat(sprintf("formule %.6f  numerique %.6f\n",.bay_t_prior_theta(parts,r,3),log(num)))
cat("(3) .bay_t_parts contre dt brut, k=3 het H0\n")
X <- matrix(c(.1,.2,.3,-.4,1.5),1)
p <- .bay_t_parts(X,st,FALSE,TRUE,1.5)
z<-c(.3,-.4); vth<-c(exp(z),1)/(1+sum(exp(z))); nu<-2+exp(1.5); sg2<-1/(3*vth)/exp(.2); s<-sqrt(sg2*(nu-2)/nu)
ll<-sum(vapply(1:3,function(i)sum(dt((st$y[[i]]-.1)/s[i],nu,log=TRUE)-log(s[i])),0))
dir<-lgamma(4.5)-3*lgamma(1.5)+.5*sum(log(vth))
jac<-sum(log(vth)); pn<- -(nu-2)+1.5
cat(sprintf("code %.6f brut %.6f\n",p$base,ll+dir+jac+pn))
cat("(4) .bay_law_ldgd : integrale sur le simplexe k=3\n")
f <- function(a,b){ n<-nod(0,1,200); tot<-0; for(i in seq_along(n$x)) { m<-nod(0,1-n$x[i],200); th<-cbind(n$x[i],m$x,1-n$x[i]-m$x); tot<-tot+n$w[i]*sum(m$w*exp(.bay_law_ldgd(th,a,b)))}; tot}
cat(f(c(1.5,2,3),c(1,.5,4)),"\n")
cat("(5) Densite mvt : integrale d=2 et covariance\n")
fit<-.bay_mvt_make(c(1,2),c(2,.5),matrix(c(cos(.3),sin(.3),-sin(.3),cos(.3)),2))
set.seed(1);D<-.bay_mvt_draw(2e5,fit,5,1); cat("cov tirages*(df-2)/df:",round(cov(D)*3/5,3),"; attendu:",round(fit$root%*%t(fit$root),3),"\n")
gx<-nod(-40,40,400);gr<-expand.grid(gx$x,gx$x);ww<-as.vector(outer(gx$w,gx$w));cat("integrale ld:",sum(ww*exp(.bay_mvt_ld(as.matrix(gr),fit,5,1))),"\n")

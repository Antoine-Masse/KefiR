# agent_29 : controle separe par le coordinateur de la voie hetero.
# (1) nouveau tirage H0, autres graines ; (2) BF a k=2 par Monte-Carlo direct depuis l'a priori sur (sigma1/sigma2, g, a),
# ecrit a part : log BF = log E[ L(mu integre, a) ] via grille sur rho et delta ~ Cauchy.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
R <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/R"; for (f in c("sys_msg.R","sys_exit.R","sys_bay_core.R","sys_bay_oneway.R")) source(file.path(R,f))
set.seed(909090)
n <- c(5,5,40); s <- c(4,4,1); g <- factor(rep(1:3,n)); N <- 300
P <- replicate(N, { y <- rnorm(50,0,rep(s,n)); plogis(.bay_bf_oneway_hetero(y,g,n_draws=2000L,seed=1L)$log_bf) })
k <- sum(P>.95); cat("H0 4,4,1 N",N,"P>.95",k/N,"IC",round(binom.test(k,N)$conf.int,3),"P>.5",mean(P>.5),"\n")
# (2) k=2 : marginale de (y1,y2) sous H0 et H1 par integration directe en (log s1, log s2) et delta ; mu integre analytiquement
# H(a) = |.| avec V_i=s_i^2/n_i, moyenne d=m1-m2 ~ N(delta*tau, V1+V2) ; a mu plat seul d compte.
bf2 <- function(y1,y2,r=sqrt(2)/2) {
  n1<-length(y1);n2<-length(y2);d<-mean(y1)-mean(y2);S1<-sum((y1-mean(y1))^2);S2<-sum((y2-mean(y2))^2)
  ls <- seq(log(min(S1/n1,S2/n2))/2-6, log(max(S1/n1,S2/n2))/2+6, length.out=301); h<-ls[2]-ls[1]
  gr <- expand.grid(l1=ls,l2=ls); s1<-exp(gr$l1); s2<-exp(gr$l2); V<-s1^2/n1+s2^2/n2; tau<-sqrt((s1^2+s2^2)/2)
  lw <- -(n1-1)*gr$l1-S1/(2*s1^2)-(n2-1)*gr$l2-S2/(2*s2^2)-log(V)/2
  m0 <- exp(lw-max(lw))*exp(-d^2/(2*V))
  # H1 : delta~Cauchy(r) : integrale 1D par changement t=atan
  u <- (seq_len(400)-.5)/400*pi-pi/2; del <- r*tan(u)
  m1 <- numeric(nrow(gr)); keep<-which(lw-max(lw)> -40)
  for(i in keep) m1[i] <- mean(exp(-(d-del*tau[i])^2/(2*V[i])))
  log(sum(exp(lw[keep]-max(lw))*m1[keep]))-log(sum(m0[keep]))
}
set.seed(11); for (cas in list(c(10,10,1,1,0),c(6,20,3,1,1),c(4,4,2,1,2))) {
  y1<-rnorm(cas[1],cas[5],cas[3]);y2<-rnorm(cas[2],0,cas[4]);G<-factor(rep(1:2,c(cas[1],cas[2])))
  a<-.bay_bf_oneway_hetero(c(y1,y2),G,n_draws=50000L)
  cat("cas",cas," fonction",round(a$log_bf,4),"err%",round(a$error_pct,3)," coord",round(bf2(y1,y2),4),"\n") }

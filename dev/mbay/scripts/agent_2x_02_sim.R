# Etapes 1 a 3 : k=2 et k=3. Voie rapide validee dans agent_2x_01_valid.R (ecart anovaBF/ttestBF < 1e-6 sur le log BF).
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
args <- commandArgs(TRUE); N1 <- as.numeric(args[1]); N2 <- as.numeric(args[2]); N3 <- as.numeric(args[3]); NH <- as.numeric(args[4])
BUDGET <- as.numeric(args[5]); T0 <- Sys.time(); el <- function() as.numeric(Sys.time()-T0, units="secs")
logt <- function(...) cat(format(Sys.time(),"%H:%M:%S"), ..., "\n")
sc <- 0
# ---- Etape 1 : sous l'a priori
res1 <- list()
for (pi0 in c(.5,.2)) for (k in c(3,2)) for (n in c(10,30)) {
  sc <- sc+1; set.seed(20261100+sc); .rep$n <- 0L; t1 <- el()
  H <- rbinom(N1,1,pi0); l <- numeric(N1); F <- numeric(N1); rp <- logical(N1)
  for (i in 1:N1) { al <- if (H[i]) ralpha(k) else rep(0,k)
    y <- rnorm(n*k) + rep(al, each=n); F[i] <- Fstat(y,n,k); b <- .rep$n; l[i] <- lbf_F(F[i],n,k); rp[i] <- .rep$n > b }
  res1[[length(res1)+1]] <- list(pi=pi0,k=k,n=n,H=H,l=l,F=F,rp=rp,sec=el()-t1)
  logt("E1 pi",pi0,"k",k,"n",n,"s",round(el()-t1),"replis",.rep$n)
  saveRDS(res1, file.path(RES,"agent_2x_E1.rds"))
}
# ---- Etape 2 : effets fixes.
# ECHELLE : delta = (moyenne max - moyenne min)/sigma. Motifs alpha = (-d/2, d/2) si k=2 ; (-d/2, 0, d/2) si k=3.
# Alors ||b||^2 = sum(alpha^2) = d^2/2 pour k=2 ET k=3 : meme norme sur l'echelle de l'a priori et meme
# non-centralite n d^2/2. f de Cohen = d/2 (k=2), d/sqrt6 (k=3). Le motif est sans effet sur F en plan equilibre normal.
alph <- function(k,d) if (k==2) c(-d/2,d/2) else c(-d/2,0,d/2)
rerr <- list(normal=function(m) rnorm(m),
             t3=function(m) rt(m,3)/sqrt(3),                                     # variance 1
             lognormale=function(m) (rlnorm(m)-exp(.5))/sqrt((exp(1)-1)*exp(1))) # moyenne 0, variance 1
run2 <- function(k,n,d,err,N) { .rep$n <- 0L; a <- rep(alph(k,d),each=n); e <- rerr[[err]]; F <- numeric(N); l <- numeric(N)
  for (i in 1:N) { F[i] <- Fstat(a+e(n*k),n,k); l[i] <- lbf_F(F[i],n,k) }
  list(k=k,n=n,d=d,err=err,l=l,p=pF(F,n,k),rep=.rep$n) }
res2 <- list()
for (k in c(3,2)) for (n in c(10,30,100)) for (d in c(0,.2,.5,.8)) { if (el()>BUDGET) {logt("BUDGET E2",k,n,d); next}
  sc <- sc+1; set.seed(20261100+sc); res2[[length(res2)+1]] <- run2(k,n,d,"normal",N2) }
logt("E2 normal fini"); saveRDS(res2, file.path(RES,"agent_2x_E2.rds"))
for (err in c("t3","lognormale")) for (k in c(3,2)) for (d in c(0,.5)) { if (el()>BUDGET) {logt("BUDGET E2nn",err,k,d); next}
  sc <- sc+1; set.seed(20261100+sc); res2[[length(res2)+1]] <- run2(k,30,d,err,N2) }
logt("E2 non-normal fini"); saveRDS(res2, file.path(RES,"agent_2x_E2.rds"))
# Heteroscedasticite : k=3, n=10,30,50, sigma de reference = 1, motif (-d/2,0,d/2). anovaBF complet (plan desequilibre).
resH <- list(); ns <- c(10,30,50); f <- factor(rep(1:3,ns))
for (cas in list(list(nom="SD 1,1,3 (SD 3 sur n=50)",sd=c(1,1,3)), list(nom="SD 3,1,1 (SD 3 sur n=10)",sd=c(3,1,1)))) for (d in c(0,.5)) {
  if (el()>BUDGET) {logt("BUDGET H",cas$nom,d); next}
  sc <- sc+1; set.seed(20261100+sc); l <- p <- numeric(NH); t1 <- el()
  for (i in 1:NH) { y <- rep(alph(3,d),ns) + rnorm(90)*rep(cas$sd,ns); l[i] <- lbf_anova(y,f); p[i] <- paov(y,f) }
  resH[[length(resH)+1]] <- list(nom=cas$nom,d=d,l=l,p=p,sec=el()-t1); saveRDS(resH, file.path(RES,"agent_2x_EH.rds"))
}
logt("heterosc fini")
# ---- Etape 3 : r, memes jeux de donnees
res3 <- list()
for (k in c(3,2)) for (n in if (k==3) c(10,30,100) else 30) for (d in c(0,.5)) { if (el()>BUDGET) {logt("BUDGET E3",k,n,d); next}
  sc <- sc+1; set.seed(20261100+sc); a <- rep(alph(k,d),each=n); L <- matrix(NA,N3,3); P <- numeric(N3)
  for (i in 1:N3) { F <- Fstat(a+rnorm(n*k),n,k); P[i] <- pF(F,n,k); L[i,] <- sapply(c(.25,.5,1), function(r) lbf_F(F,n,k,r)) }
  res3[[length(res3)+1]] <- list(k=k,n=n,d=d,L=L,p=P); saveRDS(res3, file.path(RES,"agent_2x_E3.rds"))
}
logt("FINI s=",round(el()))

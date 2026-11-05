# Etape 4 : A*B equilibre 2x3, n=10 par cellule. Simulation STRATIFIEE par modele vrai : chaque generateur
# (melange de modeles) s'analyse ensuite par reponderation exacte des strates (voir agent_2x_04_AxB_ana.R).
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
args <- commandArgs(TRUE); max_s <- as.numeric(args[1]); N_max <- as.numeric(args[2])
set.seed(20261142)
d <- expand.grid(i = 1:10, A = gl(2,1), B = gl(3,1))
PA <- as.matrix(proj(2)); PB <- as.matrix(proj(3))
XA <- as.matrix(model.matrix(~A-1,d)) %*% PA; XB <- as.matrix(model.matrix(~B-1,d)) %*% PB
XAB <- cbind(XA[,1]*XB[,1], XA[,1]*XB[,2])   # = rowMultiply de BayesFactor, controle dans agent_2x_00
M <- list(null=c(0,0,0), A=c(1,0,0), B=c(0,1,0), AplusB=c(1,1,0), AxB=c(1,1,1))
noms <- c("A","B","A + B","A + B + A:B")
out <- list(); t0 <- Sys.time(); i <- 0
sv <- function() saveRDS(do.call(rbind,out), file.path(RES,"agent_2x_04_AxB.rds"))
while (as.numeric(Sys.time()-t0,units="secs") < max_s && i < N_max) {
  i <- i+1; m <- ((i-1) %% 5)+1; z <- M[[m]]; mu <- rep(0,60)
  if (z[1]) mu <- mu + XA  %*% rnorm(1,0,sqrt(rg(.5)))
  if (z[2]) mu <- mu + XB  %*% rnorm(2,0,sqrt(rg(.5)))
  if (z[3]) mu <- mu + XAB %*% rnorm(2,0,sqrt(rg(.5)))
  d$y <- as.numeric(mu) + rnorm(60)
  o <- tryCatch(anovaBF(y~A*B, data=d, whichModels="withmain", rscaleFixed=.5, progress=FALSE)@bayesFactor,
                error=function(e) NULL)
  if (is.null(o)) { l <- rep(NA_real_,4); er <- rep(NA_real_,4) } else { l <- o[noms,"bf"]; er <- o[noms,"error"] }
  nrep <- 0L
  ssT <- sum((d$y-mean(d$y))^2); nb <- 0
  if (is.na(l[1])) { nrep<-nrep+1L; R2 <- summary(lm(y~A,d))$r.squared; l[1] <- lbf_repli(R2,60,1,30) }
  if (is.na(l[2])) { nrep<-nrep+1L; R2 <- summary(lm(y~B,d))$r.squared; l[2] <- lbf_repli(R2,60,2,20) }
  out[[i]] <- c(i=i, model=m, l0=0, lA=l[1], lB=l[2], lAB=l[3], lAxB=l[4], eAB=er[3], eAxB=er[4], nrep=nrep,
                sec=as.numeric(Sys.time()-t0,units="secs"))
  if (i %% 25 == 0) sv()
}
sv(); cat("termine, jeux :",i," s :",as.numeric(Sys.time()-t0,units="secs"),"\n")

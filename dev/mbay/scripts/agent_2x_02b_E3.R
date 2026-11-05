# Complement etape 3 (scenarios coupes par le budget de agent_2x_02_sim.R), avec moins de jeux pour k=2.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/agent_2x_commun.R")
alph <- function(k,d) if (k==2) c(-d/2,d/2) else c(-d/2,0,d/2)
res3 <- readRDS(file.path(RES,"agent_2x_E3.rds")); sc <- 20261300
for (cf in list(c(3,100,600),c(2,30,300))) for (d in c(0,.5)) { k<-cf[1];n<-cf[2];N3<-cf[3]
  sc <- sc+1; set.seed(sc); a <- rep(alph(k,d),each=n); L <- matrix(NA,N3,3); P <- numeric(N3)
  for (i in 1:N3) { F <- Fstat(a+rnorm(n*k),n,k); P[i] <- pF(F,n,k); L[i,] <- sapply(c(.25,.5,1), function(r) lbf_F(F,n,k,r)) }
  res3[[length(res3)+1]] <- list(k=k,n=n,d=d,L=L,p=P); saveRDS(res3, file.path(RES,"agent_2x_E3.rds")) }
cat("E3b fini\n")

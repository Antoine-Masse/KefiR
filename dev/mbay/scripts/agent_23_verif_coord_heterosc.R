# agent_23 : recalcul separe par le coordinateur du taux de P>0,95 (anovaBF one-way, H0, n=(5,5,40), sd=(4,4,1)).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); set.seed(777001)
n <- c(5,5,40); s <- c(4,4,1); g <- factor(rep(1:3,n)); N <- 1500
P <- replicate(N, { y <- rnorm(50,0,rep(s,n)); bf <- extractBF(anovaBF(y~g,data=data.frame(y,g),progress=FALSE),onlybf=TRUE); bf/(1+bf) })
k <- sum(P>.95); cat("N",N,"P>0.95:",k/N,"IC",binom.test(k,N)$conf.int,"P>0.5:",mean(P>.5),"\n")
r <- readRDS("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/resultats/agent_21_heterosc.rds")
cat("Relecture rds agent_21 $A : N",nrow(r$A),"P>0.95",mean(r$A$P>.95),"P>0.5",mean(r$A$P>.5),"aov",mean(r$A$p_aov<.05),"Welch",mean(r$A$p_welch<.05),"\n")

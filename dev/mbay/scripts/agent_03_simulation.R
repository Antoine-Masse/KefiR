# 3. Simulation : 200 jeux y~F (3 groupes, n=15/groupe), H0 et d=0.5, P(H1|D) avec a priori 1/2.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
# "d=0.5" : ecart de 0.5 sd entre groupes extremes voisins, moyennes -0.5, 0, 0.5
set.seed(301); f <- gl(3,15)
sim <- function(delta) t(replicate(200,{
  y <- rnorm(45)+delta*(as.integer(f)-2); dd <- data.frame(y,f)
  bf <- extractBF(anovaBF(y~f,dd,progress=FALSE))$bf
  c(P=bf/(1+bf), p=summary(aov(y~f,dd))[[1]][1,5])
}))
t0 <- Sys.time(); H0 <- sim(0); H1 <- sim(.5); el <- as.numeric(difftime(Sys.time(),t0,units="secs"))
cat(sprintf("Temps 400 jeux : %.1f s (%.0f ms/jeu)\n",el,el/400*1000))
rep_ <- function(nm,r) { cat(sprintf("\n%s : P(H1|D) quantiles 5/25/50/75/95 = %s\n", nm, paste(round(quantile(r[,1],c(.05,.25,.5,.75,.95)),3),collapse=" / ")))
  cat(sprintf("  P>0.5 : %.3f ; P>0.95 : %.3f ; P<0.05 : %.3f ; p<0.05 aov : %.3f\n",mean(r[,1]>.5),mean(r[,1]>.95),mean(r[,1]<.05),mean(r[,2]<.05))) }
rep_("H0",H0); rep_("Effet d=0.5 (moyennes -0.5, 0, +0.5)",H1)
cat("\nIC binomiaux 95 % (P>0.5 sous H0) :", round(binom.test(sum(H0[,1]>.5),200)$conf.int,3), " ; p<0.05 sous H0 :", round(binom.test(sum(H0[,2]<.05),200)$conf.int,3),"\n")
cat("IC 95 % sous effet : P>0.5 :", round(binom.test(sum(H1[,1]>.5),200)$conf.int,3), " ; P>0.95 :", round(binom.test(sum(H1[,1]>.95),200)$conf.int,3), " ; p<0.05 :", round(binom.test(sum(H1[,2]<.05),200)$conf.int,3),"\n")
cat("Correlation rang(P, -p) sous effet :", round(cor(H1[,1],-H1[,2],method="spearman"),3),"\n")
cat("Seuil de p correspondant a P=0.5 (valeurs voisines) :", signif(max(H1[H1[,1]>.5,2]),3),"\n")
cat("Seuil de p correspondant a P=0.95 :", signif(max(H1[H1[,1]>.95,2]),3),"\n")
cat("\nFIN\n")

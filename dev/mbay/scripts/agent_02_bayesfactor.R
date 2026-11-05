# 2a-d. BayesFactor : verifications numeriques independantes, BF d'inclusion, plans, limites, temps.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor)); options(BFprogress=FALSE)
cat("BayesFactor", as.character(packageVersion("BayesFactor")), "\n")
tryit <- function(label, expr) { cat("\n--", label, "\n"); r <- tryCatch({ v <- withCallingHandlers(expr, warning=function(w){cat("  [warning]",conditionMessage(w),"\n"); invokeRestart("muffleWarning")}); v }, error=function(e){ cat("  ERREUR:", conditionMessage(e), "\n"); NULL }); r }
bfx <- function(b) { x <- extractBF(b); x[,c("bf","error")] }

cat("\n=== 2a. ttestBF vs integration JZS (Rouder 2009) ===\n")
jzs10 <- function(t, N, nu, r=sqrt(2)/2) {
  f <- function(g) (1+N*g*r^2)^(-1/2) * (1 + t^2/((1+N*g*r^2)*nu))^(-(nu+1)/2) * (2*pi)^(-1/2) * g^(-3/2) * exp(-1/(2*g))
  integrate(f, 0, Inf, rel.tol=1e-10)$value / (1+t^2/nu)^(-(nu+1)/2)
}
for (cas in list(list(s=201,n1=20,n2=20,d=0), list(s=202,n1=20,n2=25,d=.6), list(s=203,n1=100,n2=100,d=.3))) {
  set.seed(cas$s); x <- rnorm(cas$n1); y <- rnorm(cas$n2)+cas$d; tt <- t.test(y,x,var.equal=TRUE)
  m <- jzs10(tt$statistic, cas$n1*cas$n2/(cas$n1+cas$n2), cas$n1+cas$n2-2); b <- bfx(ttestBF(y,x))
  cat(sprintf("  n=%d,%d d=%.1f t=%.4f : ttestBF=%.8g (propErr %.2g) ; integrate=%.8g ; ecart relatif=%.2e\n",cas$n1,cas$n2,cas$d,tt$statistic,b$bf,b$error,m,abs(b$bf-m)/m))
}

cat("\n=== 2b. anovaBF 3 groupes vs Zellner-Siow programme a la main (MEME a priori : Rouder 2012, g~IG(1/2,r^2/2), contrastes orthonormes, mu plat, sigma Jeffreys) ===\n")
Qs <- function(a) { e <- eigen(diag(a)-1/a); e$vectors[,seq_len(a-1),drop=FALSE] }
Xf <- function(f) { f <- factor(f); model.matrix(~0+f) %*% Qs(nlevels(f)) }
bf_g <- function(y, Xs, gidx, gs) { # BF(g) en log, X centre
  n <- length(y); yc <- y-mean(y); X <- scale(Xs, center=TRUE, scale=FALSE)
  G <- diag(gs[gidx], ncol(X)); M <- crossprod(X) + diag(1/gs[gidx], ncol(X))
  ld <- as.numeric(determinant(diag(ncol(X)) + G %*% crossprod(X))$modulus)
  q <- sum(yc^2) - drop(crossprod(crossprod(X,yc), solve(M, crossprod(X,yc))))
  -ld/2 - (n-1)/2*(log(q)-log(sum(yc^2)))
}
dig <- function(g, r) (r^2/2)^(1/2)/gamma(1/2) * g^(-3/2) * exp(-r^2/(2*g))
zs1 <- function(y, f, r=.5) { X <- Xf(f); integrate(Vectorize(function(g) exp(bf_g(y,X,rep(1,ncol(X)),g))*dig(g,r)),0,Inf,rel.tol=1e-8)$value }
for (cas in list(list(s=211,n=15,eff=0), list(s=212,n=15,eff=.5), list(s=213,n=c(10,15,20),eff=.5))) {
  set.seed(cas$s); f <- factor(rep(1:3, times=cas$n)); if(length(cas$n)==1) f <- gl(3,cas$n); y <- rnorm(length(f)) + cas$eff*(as.integer(f)-2)
  b <- bfx(anovaBF(y~f,data.frame(y,f))); m <- zs1(y,f)
  cat(sprintf("  n=%s eff=%.1f : anovaBF=%.8g (propErr %.2g) ; integrate ZS=%.8g ; ecart relatif=%.2e\n",paste(cas$n,collapse=","),cas$eff,b$bf,b$error,m,abs(b$bf-m)/m))
}
set.seed(214); dd <- expand.grid(A=gl(2,1),B=gl(3,1),r=1:8); dd$y <- rnorm(nrow(dd))+as.integer(dd$A)*.5+as.integer(dd$B)*.3
XA <- Xf(dd$A); XB <- Xf(dd$B)
mab <- integrate(Vectorize(function(g1) integrate(Vectorize(function(g2) exp(bf_g(dd$y,cbind(XA,XB),c(1,2,2),c(g1,g2)))*dig(g2,.5)),0,Inf,rel.tol=1e-5)$value*dig(g1,.5)),0,Inf,rel.tol=1e-5)$value
bab <- bfx(lmBF(y~A+B,dd,iterations=50000))
cat(sprintf("  y~A+B (2 g, 2D) : lmBF=%.6g (propErr %.2g) ; double integrate=%.6g ; ecart relatif=%.2e\n",bab$bf,bab$error,mab,abs(bab$bf-mab)/mab))
set.seed(215); n <- 40; x <- rnorm(n); y <- .4*x+rnorm(n); R2 <- summary(lm(y~x))$r.squared; r <- sqrt(2)/4
zc <- integrate(function(g) exp((n-2)/2*log(1+g)-(n-1)/2*log(1+g*(1-R2)))*(r^2*n/2)^(1/2)/gamma(1/2)*g^(-3/2)*exp(-r^2*n/(2*g)),0,Inf)$value
bc <- bfx(lmBF(y~x,data.frame(x,y)))
cat(sprintf("  y~x continu (Liang, g~IG(1/2,r^2 n/2), r=sqrt(2)/4) : lmBF=%.6g ; integrate=%.6g ; ecart relatif=%.2e\n",bc$bf,zc,abs(bc$bf-zc)/zc))

cat("\n=== 2c. y~A*B : BF d'inclusion a la main ===\n")
set.seed(221); d <- expand.grid(A=gl(2,1,labels=c("a1","a2")),B=gl(3,1,labels=c("b1","b2","b3")),r=1:10); d$y <- rnorm(nrow(d)) + .7*(d$A=="a2") + .6*(d$A=="a2"&d$B=="b3")
t0 <- Sys.time(); ba <- anovaBF(y~A*B,d,whichModels="all"); cat("  Temps whichModels=all :", format(difftime(Sys.time(),t0),digits=3),"\n")
x <- extractBF(ba); x$modele <- rownames(x); print(x[,c("modele","bf","error")],row.names=FALSE)
mods <- c(Intercept=1, setNames(x$bf,gsub(" ","",x$modele)))
has <- function(m,t) t %in% strsplit(m,"[+]")[[1]]
incl_all <- function(t) { w <- sapply(names(mods),has,t); sum(mods[w])/sum(mods[!w]) }
get <- function(...) { k <- c(...); f <- function(m) setequal(strsplit(m,"[+]")[[1]],k); mods[sapply(names(mods),f)] }
i0 <- 1
incA <- (get("A")+get("A","B"))/(i0+get("B")); incB <- (get("B")+get("A","B"))/(i0+get("A")); incI <- get("A","B","A:B")/get("A","B")
ia <- sapply(c("A","B","A:B"),incl_all)
cat(sprintf("  Inclusion 'tous modeles' (8 modeles) : A=%.4g B=%.4g A:B=%.4g\n",ia[1],ia[2],ia[3]))
cat(sprintf("  Inclusion 'apparies'  : A=%.4g B=%.4g A:B=%.4g\n",incA,incB,incI))
if (requireNamespace("bayestestR",quietly=TRUE)) {
  r <- tryCatch(as.data.frame(bayestestR::bayesfactor_inclusion(ba,match_models=FALSE)),error=function(e) conditionMessage(e)); cat("  bayestestR match_models=FALSE:\n"); print(r)
  r <- tryCatch(as.data.frame(bayestestR::bayesfactor_inclusion(ba,match_models=TRUE)),error=function(e) conditionMessage(e)); cat("  bayestestR match_models=TRUE:\n"); print(r)
}
print(anova(lm(y~A*B,d)))

cat("\n-- ANCOVA y~F+x : lmBF(y~F+x)/lmBF(y~x)\n")
set.seed(222); dc <- data.frame(F=gl(3,20)); dc$x <- rnorm(60); dc$y <- .5*(as.integer(dc$F)-2)+.6*dc$x+rnorm(60)
r <- lmBF(y~F+x,dc)/lmBF(y~x,dc); print(bfx(r)); print(anova(lm(y~x+F,dc))[2,])

cat("\n-- Apparie : lmBF whichRandom=id vs ttestBF(paired=TRUE)\n")
set.seed(223); ni <- 20; id <- factor(rep(1:ni,2)); tm <- gl(2,ni); ui <- rnorm(ni); yy <- ui[id]+.5*(tm==2)+rnorm(40,sd=.7); dp <- data.frame(y=yy,tm,id)
tp <- bfx(lmBF(y~tm+id,dp,whichRandom="id")/lmBF(y~id,dp,whichRandom="id")); print(tp)
tt <- bfx(ttestBF(dp$y[tm==2],dp$y[tm==1],paired=TRUE)); print(tt); print(t.test(dp$y[tm==2],dp$y[tm==1],paired=TRUE)$p.value)

cat("\n-- Split-plot : y~G*T+id, id aleatoire\n")
set.seed(224); ns <- 16; sp <- expand.grid(T=gl(3,1),id=factor(1:ns)); sp$G <- factor(ifelse(as.integer(sp$id)<=ns/2,"g1","g2")); u <- rnorm(ns); sp$y <- u[sp$id]+.8*(sp$T==3)+.5*(sp$G=="g2")+rnorm(nrow(sp),sd=.7)
t0 <- Sys.time(); bs <- anovaBF(y~G*T+id,sp,whichRandom="id",whichModels="withmain"); print(bfx(bs)); cat("  temps:",format(difftime(Sys.time(),t0),digits=3),"\n")
r <- lmBF(y~G*T+id,sp,whichRandom="id")/lmBF(y~G+T+id,sp,whichRandom="id"); cat("  A:B apparie :\n"); print(bfx(r))

cat("\n=== 2d. Limites : on essaie ===\n")
set.seed(231); dm <- data.frame(F=gl(3,20)); dm$y1 <- rnorm(60); dm$y2 <- rnorm(60); dm$id <- factor(rep(1:20,3)); dm$x <- rnorm(60); dm$bin <- rbinom(60,1,.5); dm$cnt <- rpois(60,3)
print(tryit("cbind multivarie anovaBF", anovaBF(cbind(y1,y2)~F,dm)))
print(tryit("cbind multivarie lmBF", lmBF(cbind(y1,y2)~F,dm)))
print(tryit("pente aleatoire (1+x|id) lmBF", lmBF(y1~x+(1+x|id),dm)))
print(tryit("pente aleatoire par x:id whichRandom (interaction continue-facteur)", lmBF(y1~x+id+x:id,dm,whichRandom=c("id","x:id"))))
print(tryit("pente aleatoire par T:id (facteur:id, ce que BayesFactor sait faire)", lmBF(y~G+T+id+T:id,sp,whichRandom=c("id","T:id"))))
cat("\n-- variances heterogenes : aucun argument dans lmBF/anovaBF :\n"); print(names(formals(lmBF))); print(names(formals(anovaBF)))
print(tryit("var.equal dans ttestBF", ttestBF(dm$y1[1:20],dm$y1[21:40],var.equal=FALSE)))
print(tryit("family=binomial lmBF", lmBF(bin~F,dm,family=binomial())))
print(tryit("y binaire lmBF traite comme gaussien", lmBF(bin~F,dm)))
print(tryit("y comptage anovaBF traite comme gaussien", anovaBF(cnt~F,dm)))
print(tryit("generalTestBF multivarie", generalTestBF(cbind(y1,y2)~F,dm)))
print(tryit("contingencyTableBF (binaire en table, seul non gaussien)", contingencyTableBF(table(dm$F,dm$bin),sampleType="jointMulti")))
print(tryit("donnees manquantes", {dna <- dm; dna$y1[1] <- NA; lmBF(y1~F,dna)}))
print(tryit("y~F desequilibre avec une cellule a 1 obs", {dz<-dm[c(1,21:60),]; lmBF(y1~F,dz)}))

cat("\n=== 2c bis. Temps anovaBF(whichModels='withmain'), limite 60 s/cas ===\n")
tcell <- function(n, k) {
  set.seed(n+k); L <- rep(list(gl(2,1)),k); names(L) <- LETTERS[seq_len(k)]
  g <- expand.grid(L); dd <- g[rep(seq_len(nrow(g)),length.out=n),,drop=FALSE]; dd$y <- rnorm(n)+.3*(dd$A==2)
  fm <- as.formula(paste("y~",paste(LETTERS[seq_len(k)],collapse="*")))
  t0 <- Sys.time()
  setTimeLimit(elapsed=60,transient=TRUE)
  on.exit(setTimeLimit())
  r <- tryCatch(anovaBF(fm,dd,progress=FALSE),error=function(e) NULL)
  setTimeLimit()
  el <- as.numeric(difftime(Sys.time(),t0,units="secs"))
  if (is.null(r)) return(sprintf("  n=%3d facteurs=%d : ECHEC/limite 60 s depassee (%.1f s)",n,k,el))
  x <- extractBF(r); sprintf("  n=%3d facteurs=%d : %.2f s ; %d modeles ; propErr max=%.3g ; mediane=%.3g",n,k,el,nrow(x),max(x$error,na.rm=TRUE),median(x$error,na.rm=TRUE))
}
for (n in c(30,100,500)) for (k in 1:4) { cat(tcell(n,k),"\n"); flush.console() }
t0 <- Sys.time(); dd <- expand.grid(A=gl(2,1),B=gl(2,1),C=gl(2,1),r=1:13); dd$y<-rnorm(nrow(dd)); r<-tryCatch({setTimeLimit(elapsed=30,transient=TRUE); anovaBF(y~A*B*C,dd,whichModels="all")},error=function(e)NULL);setTimeLimit()
cat(sprintf("  n=%d 3 facteurs whichModels=all : %s %.1f s\n",nrow(dd),if(is.null(r)) "ECHEC/limite 30 s" else paste(length(r),"modeles"),as.numeric(difftime(Sys.time(),t0,units="secs"))))
cat("\nFIN\n"); print(sessionInfo())

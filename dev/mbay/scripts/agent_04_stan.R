# 2e. rstanarm (modeles precompiles) : temps, R-hat, ESS, bridge sampling, multivarie ; brms sans compilateur.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
cat("g++:",Sys.which("g++"),"\n")
Tm <- function(label, expr) { t0 <- Sys.time(); cat("\n--",label,"\n"); r <- tryCatch(withCallingHandlers(expr, warning=function(w){cat("  [warning]",substr(conditionMessage(w),1,200),"\n"); invokeRestart("muffleWarning")}), error=function(e){cat("  ERREUR:",substr(conditionMessage(e),1,800),"\n"); NULL}); cat(sprintf("  temps: %.1f s\n",as.numeric(difftime(Sys.time(),t0,units="secs")))); flush.console(); r }
diag_ <- function(fit) if(!is.null(fit)) { s <- summary(fit); cat(sprintf("  Rhat max=%.3f ; n_eff min=%.0f\n",max(s[,"Rhat"],na.rm=TRUE),min(s[,"n_eff"],na.rm=TRUE))) }
suppressMessages(ok <- require(rstanarm)); cat("rstanarm charge:",ok,"\n")
if (ok) {
  options(mc.cores=1)
  set.seed(401); d <- data.frame(F=gl(3,20)); d$y <- .5*(as.integer(d$F)-2)+rnorm(60); d$x <- rnorm(60); d$id <- factor(rep(1:20,3)); d$y2 <- d$y*.5+rnorm(60); d$tm <- d$F
  d$yl <- d$y + rnorm(20)[d$id]
  a1 <- Tm("stan_aov n=60 (2 chaines 1000 it)", stan_aov(y~F,d,prior=R2(.5),chains=2,iter=1000,refresh=0,seed=1)); diag_(a1)
  g1 <- Tm("stan_glm n=60 (2 chaines 1000 it)", stan_glm(y~F,data=d,chains=2,iter=1000,refresh=0,seed=1)); diag_(g1)
  g2 <- Tm("stan_lmer (1|id) n=60 (2 chaines 1000 it)", stan_lmer(yl~F+(1|id),data=d,chains=2,iter=1000,refresh=0,seed=1)); diag_(g2)
  cat("\n=== bridge sampling sur rstanarm ===\n")
  if (requireNamespace("bridgesampling",quietly=TRUE)) {
    b1 <- Tm("bridge_sampler(stanreg) stan_glm", bridgesampling::bridge_sampler(g1,silent=TRUE)); print(b1)
    b2 <- Tm("bridge_sampler(fit$stanfit)", bridgesampling::bridge_sampler(g1$stanfit,silent=TRUE)); print(b2)
    g1b <- Tm("stan_glm diagnostic_file puis bridge", stan_glm(y~F,data=d,chains=2,iter=2000,refresh=0,seed=1,diagnostic_file=file.path(tempdir(),"df.csv")))
    b3 <- Tm("bridge_sampler(stanreg avec diagnostic_file)", bridgesampling::bridge_sampler(g1b,silent=TRUE)); print(b3)
    g0b <- Tm("modele nul stan_glm", stan_glm(y~1,data=d,chains=2,iter=2000,refresh=0,seed=1,diagnostic_file=file.path(tempdir(),"df0.csv")))
    b4 <- Tm("bridge_sampler nul", bridgesampling::bridge_sampler(g0b,silent=TRUE)); print(b4)
    if (!is.null(b3)&&!is.null(b4)) { cat("  BF10 rstanarm+bridge :\n"); print(bridgesampling::bf(b3,b4)); cat("  erreur:",bridgesampling::error_measures(b3)$percentage,bridgesampling::error_measures(b4)$percentage,"\n") }
    suppressMessages(library(BayesFactor)); cat("  anovaBF meme donnees (A PRIORI DIFFERENTS) :\n"); print(extractBF(anovaBF(y~F,d,progress=FALSE))[,1:2])
  }
  mv <- Tm("stan_mvmer (2 reponses, (1|id))", stan_mvmer(list(y~F+(1|id),y2~F+(1|id)),data=d,chains=1,iter=500,refresh=0,seed=1))
  if(!is.null(mv)) { s<-summary(mv); cat(sprintf("  Rhat max=%.3f ; n_eff min=%.0f\n",max(s[,"Rhat"],na.rm=TRUE),min(s[,"n_eff"],na.rm=TRUE))) }
  mv2 <- Tm("stan_mvmer sans effet aleatoire", stan_mvmer(list(y~F,y2~F),data=d,chains=1,iter=200,refresh=0,seed=1))
}
cat("\n=== brms sans compilateur ===\n")
suppressMessages(okb <- require(brms)); cat("brms charge:",okb,"\n")
if (okb) {
  set.seed(401); d <- data.frame(F=gl(3,20)); d$y <- rnorm(60); d$y2 <- rnorm(60)
  m <- Tm("brm(y~F) sans compilateur", brm(y~F,d,chains=1,iter=200,refresh=0,silent=2))
  cat("  stancode (sans compilation) :\n"); cat(substr(tryCatch(stancode(y~F,d),error=function(e)conditionMessage(e)),1,200),"...\n")
  sc <- Tm("make_stancode mvbind (sans compilation)", brms::make_stancode(bf(mvbind(y,y2)~F)+set_rescor(TRUE),d)); cat("  mvbind code genere:",!is.null(sc),"\n")
}
cat("\nFIN\n"); print(sessionInfo())

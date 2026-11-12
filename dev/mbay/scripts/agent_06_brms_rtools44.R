# brms avec Rtools 4.4 present (C:/RBuildTools/4.4, hors PATH, version ANTERIEURE a R 4.6) : mise au PATH pour ce processus seulement.
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
rt <- "C:/RBuildTools/4.4"
Sys.setenv(PATH=paste(file.path(rt,"x86_64-w64-mingw32.static.posix/bin"),file.path(rt,"usr/bin"),Sys.getenv("PATH"),sep=";"), RTOOLS45_HOME=rt)
cat("g++:",Sys.which("g++"),"\n"); print(system2("g++","--version",stdout=TRUE)[1])
print(tryCatch(pkgbuild::has_build_tools(debug=FALSE),error=function(e)conditionMessage(e)))
suppressMessages(library(brms))
set.seed(401); d <- data.frame(F=gl(3,20)); d$y <- .5*(as.integer(d$F)-2)+rnorm(60); d$y2 <- d$y*.5+rnorm(60)
Tm <- function(label, expr) { t0 <- Sys.time(); cat("\n--",label,"\n"); r <- tryCatch(expr, error=function(e){cat("  ERREUR:",substr(conditionMessage(e),1,800),"\n"); NULL}); cat(sprintf("  temps: %.1f s\n",as.numeric(difftime(Sys.time(),t0,units="secs")))); flush.console(); r }
m <- Tm("brm(y~F) n=60, 2 chaines, 1000 it, save_pars(all)", brm(y~F,d,chains=2,iter=1000,refresh=0,silent=2,seed=1,save_pars=save_pars(all=TRUE)))
if(!is.null(m)) { print(summary(m)$fixed); b <- Tm("bridge_sampler brms",bridgesampling::bridge_sampler(m,silent=TRUE)); print(b) }
mv <- Tm("brm(mvbind(y,y2)~F) n=60, 1 chaine, 500 it", brm(bf(mvbind(y,y2)~F)+set_rescor(TRUE),d,chains=1,iter=500,refresh=0,silent=2,seed=1))
if(!is.null(mv)) print(summary(mv)$fixed)
cat("\nFIN\n")

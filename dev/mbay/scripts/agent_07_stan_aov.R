# Complement : stan_aov (l'appel de agent_04 echouait avec R2(.5) mode, 2 predicteurs).
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(rstanarm))
set.seed(401); d <- data.frame(F=gl(3,20)); d$y <- .5*(as.integer(d$F)-2)+rnorm(60)
t0 <- Sys.time(); a <- stan_aov(y~F,d,prior=R2(.5,what="mean"),chains=2,iter=1000,refresh=0,seed=1)
s <- summary(a); cat(sprintf("stan_aov : %.1f s ; Rhat max=%.3f ; n_eff min=%.0f\n",as.numeric(difftime(Sys.time(),t0,units="secs")),max(s[,"Rhat"],na.rm=TRUE),min(s[,"n_eff"],na.rm=TRUE)))

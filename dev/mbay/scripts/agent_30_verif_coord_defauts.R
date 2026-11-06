# agent_30 : reproduction par le coordinateur des defauts B1, B2, I6 de RELECTURE_ADVERSE_2026-11-06.md
.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
R <- "C:/Users/masse/Desktop/KefiR/KefiR-mbay/R"; for (f in c("sys_msg.R","sys_exit.R","sys_bay_core.R","sys_bay_oneway.R")) source(file.path(R,f))
set.seed(7); y <- c(rnorm(40), rnorm(3,50,1)); g <- factor(rep(1:5,c(10,10,10,10,3)))
cat("B1 hetero", round(.bay_bf_oneway_hetero(y,g)$log_bf,3), "zs", round(.bay_bf_oneway_zs(y,g),3), "Welch p", signif(oneway.test(y~g)$p.value,3), "\n")
for (d in c(5,50,500,5000)) { yy <- y; yy[41:43] <- y[41:43]-50+d; cat(" d",d,"log BF",round(.bay_bf_oneway_hetero(yy,g)$log_bf,2),"\n") }
set.seed(3); a <- rnorm(10); b <- rnorm(10,1.3); z <- rnorm(10)
for (S in c(1,10,100,1e4,1e6)) cat("B2 S",S,"log BF",round(.bay_bf_oneway_hetero(c(a,b,z*S),factor(rep(1:3,each=10)))$log_bf,2),"\n")
cat("2 groupes seuls",round(.bay_bf_oneway_hetero(c(a,b),factor(rep(1:2,each=10)))$log_bf,2),"\n")
r <- suppressWarnings(tryCatch(.bay_bf_oneway_zs(c(.3,.3,.3,.1+.2,.1+.2,.1+.2),factor(rep(1:2,each=3))),error=function(e) conditionMessage(e))); cat("I6 zs",format(r),"\n")

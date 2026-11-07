source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
ns <- asNamespace("KefiR"); unlockBinding(".bay_law_gint", ns); orig_gint <- get(".bay_law_gint", ns)
gint_wide <- function(st, cc, rs, lo=-60, hi=130, h=0.05) {
  tg <- seq(lo, hi, by=h); L <- vapply(tg, function(t) .bay_law_kernel(st, cc, exp(t)), numeric(nrow(cc)))
  if (nrow(cc) == 1L) L <- matrix(L, nrow=1L)
  res <- matrix(NA_real_, nrow(cc), length(rs))
  for (i in seq_along(rs)) { r <- rs[i]; lp <- 0.5*log(r^2/2)-0.5*log(pi)-0.5*tg-r^2/(2*exp(tg)); A <- sweep(L,2L,lp,"+"); top <- apply(A,1L,max); res[,i] <- top+log(rowSums(exp(A-top)))+log(h) }
  res
}
mk <- function(n,m,s,seed=1) { set.seed(seed); y<-unlist(Map(function(a,b,c) { v<-rnorm(a); if(a>1) v<-(v-mean(v))/sd(v); b+c*v },n,m,s)); list(y=y,g=rep(seq_along(n),n)) }
cat("=== B2. d/s 1e8, 1e10, code seul\n")
for (ratio in c(1e7,1e8,1e10)) {
  d <- mk(c(10,10,10),c(0,1,.5),c(1/ratio,1/ratio,100)); st <- .bay_laws_stats(d$y,d$g)
  cat("ratio",ratio,"m",signif(st$m,4),"ss",signif(st$ss,4),"\n")
  o <- tryCatch(tab_of(d$y,d$g,draws=1000L,laws="normal_het"),error=function(e) conditionMessage(e)); if(is.character(o)) cat(" code ERREUR",o,"\n") else cat(" code logBF",o$log_ml[2]-o$log_ml[1],"H1",o$log_ml[2],"err",o$rel_se[2],"\n")
  assign(".bay_law_gint",function(st,cc,rs,h=NULL) gint_wide(st,cc,rs),ns)
  o2 <- tryCatch(tab_of(d$y,d$g,draws=1000L,laws="normal_het"),error=function(e) conditionMessage(e)); if(is.character(o2)) cat(" large ERREUR",o2,"\n") else cat(" large H1",o2$log_ml[2],"\n")
  assign(".bay_law_gint",orig_gint,ns); flush.console()
}
cat("\n=== C. Pilote et pas (normal_het H1)\n")
plans <- list(
 "k=8 (10x7,3) d=1e3"=list(c(rep(10,7),3),c(rep(0,7),1e3),rep(1,8)),
 "k=20 n=3, d=30"=list(rep(3,20),c(rep(0,19),30),rep(1,20)),
 "k=30 n=2, sd 1"=list(rep(2,30),seq(0,1,length.out=30),rep(1,30)),
 "k=20 n=10 sd varies"=list(rep(10,20),seq(0,3,length.out=20),exp(seq(-3,3,length.out=20))),
 "k=2 N=2e5 d=.01"=list(rep(1e5,2),c(0,.01),c(1,2)),
 "k=5 n=1e4 d=3 sd var"=list(rep(1e4,5),c(0,0,0,3,0),c(1,3,.1,1,1))
)
for (nm in names(plans)) {
  d <- do.call(mk,plans[[nm]])
  t0 <- Sys.time(); a <- tab_of(d$y,d$g,draws=500L,laws="normal_het"); t1 <- Sys.time()
  assign(".bay_law_gint",function(st,cc,rs,h=NULL) gint_wide(st,cc,rs,lo=log(min(rs)^2/(1000*st$k)),hi=60,h=min(.25,.5/sqrt(st$k))/4),ns)
  b <- tab_of(d$y,d$g,draws=500L,laws="normal_het")
  assign(".bay_law_gint",orig_gint,ns)
  cat(sprintf("%-26s H1 code %.5f (ess %.0f err %.3f, %.1fs) large %.5f diff %+.2e\n",nm,a$log_ml[2],a$ess[2],a$rel_se[2],as.numeric(difftime(t1,t0,units="secs")),b$log_ml[2],b$log_ml[2]-a$log_ml[2]))
  flush.console()
}
cat("\n=== D. Invariances\n")
d <- mk(c(5,8,12),c(0,1,.3),c(1,2,.5))
base <- tab_of(d$y,d$g,draws=4000L,laws="normal_het"); cat("base",round(base$log_ml,4),"err",round(base$rel_se,4),"\n")
rf <- tab_of(-d$y,d$g,draws=4000L,laws="normal_het"); cat("y->-y",round(rf$log_ml,4),"\n")
pm <- tab_of(d$y,c(3,1,2)[d$g],draws=4000L,laws="normal_het"); cat("etiquettes",round(pm$log_ml,4),"\n")

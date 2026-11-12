# Origine exacte du NaN dans la loi normale a variances par groupe.
source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
ns <- asNamespace("KefiR")
for (f in c(".bay_law_kernel",".bay_law_kernel_grid")) unlockBinding(f, ns)
k0 <- get(".bay_law_kernel",ns); kg <- get(".bay_law_kernel_grid",ns)
seen <- new.env(); seen$n <- 0
assign(".bay_law_kernel_grid",function(st,cc,gs) {
  A <- kg(st,cc,gs)
  bad <- which(rowSums(is.nan(A))>0)
  if (length(bad) && seen$n<3) { seen$n <- seen$n+1; i <- bad[1]
    cat("  grille : lignes NaN",length(bad),"/",nrow(cc),"; c =",format(cc[i,],digits=3)," ; NaN en g =",format(range(gs[is.nan(A[i,])]),digits=3),"\n")
    B <- vapply(gs,function(g) k0(st,cc[i,,drop=FALSE],g),numeric(1))
    cat("  ancien noyau au meme point : NaN",sum(is.nan(B)),"; plage",format(range(B,na.rm=TRUE),digits=5),"\n")
    w0<-1/cc[i,];a<-1/outer(1/st$n,gs,"+");jj<-which(is.nan(A[i,]))[1];W<-sum(w0*a[,jj]);Wm<-sum(w0*st$m*a[,jj]);S2<-sum(w0*st$m^2*a[,jj])
    cat("  W",W,"Wm",Wm,"S2",S2,"S2-Wm^2/W",S2-Wm^2/W," SS/c",sum(w0*st$ss),"\n")
  }
  A
},ns)
tc <- function(lab,y,g,seed=1,laws="normal_het") { cat(lab,"\n")
  r <- tryCatch(tab_of(y,g,draws=1000L,seed=seed,laws=laws),error=function(e) {cat("  ERREUR",conditionMessage(e),"\n");NULL}); if(!is.null(r)) print(r[,c("hypothesis","log_ml","rel_se","ess")],digits=4); flush.console() }
mk <- function(n,m,s,seed=1) { set.seed(seed); y<-unlist(Map(function(a,b,c) { v<-rnorm(a); if(a>1) v<-(v-mean(v))/sd(v); b+c*v },n,m,s)); list(y=y,g=rep(seq_along(n),n)) }
y2 <- c(rep(5,10), 1,2,2,3,3,3,3,4,4,3); g2 <- rep(1:2,each=10)
tc("g1 constant a 5",y2,g2)
seen$n<-0; y3 <- c(rep(3,10), 1,2,2,3,3,3,3,4,4,5); tc("g1 constant a 3",y3,g2)
seen$n<-0; d <- mk(c(10,10,10),c(0,1,.5),c(1e-8,1e-8,100)); tc("d/s=1e8",d$y,d$g)
seen$n<-0; tc("t_het, y=3x5|1..5",c(3,3,3,3,3,1,2,3,4,5),rep(1:2,each=5),laws="t_het")

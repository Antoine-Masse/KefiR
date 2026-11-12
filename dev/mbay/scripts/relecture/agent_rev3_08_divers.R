source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev3_lib.R")
q <- function(...) m.bay(...,verbose=FALSE,draws=300L,sensitivity=FALSE)
tc <- function(lab,expr) { r <- tryCatch(expr,error=function(e) {cat(lab,": ERREUR",conditionMessage(e),"\n");NULL})
 if(!is.null(r)) cat(sprintf("%s : n=%d groupes=%s P=%.3f\n",lab,r$n,paste(r$descriptive_stats$n,collapse=","),r$prob_effect)); flush.console() }
cat("--- g recycle\n")
set.seed(1); y <- rnorm(30)
tc("g longueur 15 pour 30 y", q(y,rep(c("a","b"),length.out=15)))
tc("g longueur 10 pour 30 y", q(y,rep(c("a","b"),length.out=10)))
tc("g longueur 7 pour 30 y", q(y,rep(c("a","b"),length.out=7)))
cat("--- Origine du NaN : groupe ex-aequo, t\n")
y <- c(3,3,3,3,3, 1,2,3,4,5); g <- rep(1:2,each=5); st <- .bay_laws_stats(y,g)
set.seed(1)
for (cfg in list(c(F,F),c(T,F),c(F,T),c(T,T))) {
  a <- tryCatch(.bay_t_ml(st,cfg[1],cfg[2],1.5,.5,1000L),error=function(e) conditionMessage(e))
  cat("h1=",cfg[1]," het=",cfg[2]," : ",if(is.character(a)) a else paste(round(a$log_ml,3),round(a$rel_se,3),round(a$ess)),"\n")
}
cat("--- Crash t sur plans simples, sans ex-aequo : balayage de petits plans\n")
set.seed(100); nerr <- 0; ntot <- 0; msgs <- character()
for (it in 1:25) {
  k <- sample(2:4,1); n <- sample(c(1,2,3,5),k,TRUE); if (sum(n)<=k) next
  y <- rnorm(sum(n))*sample(c(1,10,1000),1) ; g <- rep(seq_len(k),n)
  if (runif(1)<.3) y <- round(y)
  ntot <- ntot+1
  r <- tryCatch(.bay_laws_average(tab_of(y,g,draws=300L,seed=it,laws=c("t_eq","t_het")),.5),error=function(e) e)
  if (inherits(r,"error")) { nerr <- nerr+1; msgs <- c(msgs,sprintf("n=%s y=%s : %s",paste(n,collapse=","),paste(signif(y,3),collapse=","),conditionMessage(r))) }
}
cat(ntot,"plans,",nerr,"erreurs\n"); print(head(msgs,8))

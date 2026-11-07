source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
n<-20; g<-rep(c("a","b","c"),each=n)
base <- function(seed=3){set.seed(seed); round(c(rnorm(n,0,1),rnorm(n,.3,1),rnorm(n,-.2,1)),3)}
tie_group <- function(y,grp,share,val=NULL){ id<-which(g==grp); m<-floor(share*length(id)); v<-if(is.null(val)) round(mean(y[id]),3) else val; y[id[1:m]]<-v; y}
sets <- list()
for(s in c(.45,.49,.50)) sets[[sprintf("groupe a %d%% ex-aequo",100*s)]] <- tie_group(base(),"a",s)
sets[["3 groupes 45% meme valeur"]] <- Reduce(function(y,gr) tie_group(y,gr,.45,val=0),c("a","b","c"),base())
sets[["3 groupes 50% meme valeur"]] <- Reduce(function(y,gr) tie_group(y,gr,.50,val=0),c("a","b","c"),base())
sets[["3 groupes 45% valeurs distinctes"]] <- Reduce(function(y,gr) tie_group(y,gr,.45),c("a","b","c"),base())
y<-base(); for(gr in c("a","b","c")){id<-which(g==gr); y[id[1:9]]<- -1; y[id[10:18]]<- 1}; sets[["bimodal 45+45% dans chaque groupe"]]<-y
y<-base(); y[g=="a"][1:9]<-5; sets[["groupe a 45% loin (val=5)"]]<-y
y<-base(); y[g=="a"]<-3+1e-8*rnorm(n); sets[["groupe a quasi constant (1e-8)"]]<-y
y<-base(); y[g=="a"][1:16]<-3+1e-9*rnorm(16); sets[["groupe a 80% quasi ex-aequo (1e-9)"]]<-y
y<-base(); y[1:40]<-3+1e-9*rnorm(40); sets[["2/3 echantillon quasi ex-aequo"]]<-y
y<-base(); y[g=="a"][1:11]<-3; sets[["groupe a 55% (seuil franchi)"]]<-y
for(lab in names(sets)){ y<-sets[[lab]]; st<-.bay_laws_stats(y,g)
 cat(sprintf("\n## %s : retenues=%s\n",lab,paste(.bay_laws_allowed(st,.bay_law_names()),collapse=",")))
 for(sd in 1:4){ r<-try_it(paste0(" graine ",sd), { v<-m.bay(y,g,seed=sd,verbose=FALSE,sensitivity=FALSE); cat(sprintf("   ess=%s  logBF/loi=%s  notes=%s\n",paste(round(tapply(v$models$ess,v$models$law,min)),collapse="/"),paste(round(v$laws$log_bf,2),collapse="/"),paste(substr(v$notes$fr,1,35),collapse=" | "))); v}) }
 flush.console()
}

source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
set.seed(51); y<-c(rnorm(10,5,2),rnorm(10,8,2),rnorm(10,5,2)); g<-rep(c("a","b","c"),each=10)
cat("##### k=3 defaut\n"); v<-m.bay(y,g); sp<-sqrt(sum(tapply(y,g,var)*9)/27); cat("ecart a la main",(max(tapply(y,g,mean))-min(tapply(y,g,mean)))/sp," objet ",v$largest_gap,"\n")
cat("##### k=2\n"); v2<-m.bay(y[g!="c"],g[g!="c"])
cat("##### ex-aequo : une seule loi\n"); yy<-y; yy[g=="a"][1:6]<-5; v3<-m.bay(yy,g)
cat("##### by_u, grid\n"); print(v3$sensitivity$by_u); print(v3$sensitivity$grid)
cat("##### normales seulement\n"); y4<-y; y4[c(1:5,11:15,21:25)]<-5; print(table(y4==5,g)); v4<-tryCatch(m.bay(y4,g),error=function(e)print(e))
cat("##### prior_h1=.8\n"); v5<-m.bay(y,g,prior_h1=.8)
cat("##### y*1e6+1e9\n"); v6<-m.bay(y*1e6+1e9,g,sensitivity=FALSE,verbose=FALSE); cat("gap",v6$largest_gap," P",v6$prob_effect," vs ",v$prob_effect,"\n")
cat("##### invariances (graine 1)\n")
for(lab in c("-y","labels")){ w<-if(lab=="-y") m.bay(-y,g,sensitivity=FALSE,verbose=FALSE) else m.bay(y,chartr("abc","cab",g),sensitivity=FALSE,verbose=FALSE); cat(lab,w$prob_effect,w$effects$log_bf," err ",w$effects$error_pct,"\n")}
cat("##### 2 groupes eloignes, k=20 (B1 etendu)\n")
set.seed(52); g2<-rep(sprintf("g%02d",1:20),each=3); yb<-rnorm(60); for(d in c(10,1000)){ y5<-yb; y5[g2=="g19"]<-y5[g2=="g19"]+d; y5[g2=="g20"]<-y5[g2=="g20"]-d; try_it(paste("2 eloignes d=",d),m.bay(y5,g2,verbose=FALSE,sensitivity=FALSE,draws=2000L)) ; y6<-yb; y6[g2=="g20"]<-y6[g2=="g20"]+d; try_it(paste("1 eloigne d=",d),m.bay(y6,g2,verbose=FALSE,sensitivity=FALSE,draws=2000L))}

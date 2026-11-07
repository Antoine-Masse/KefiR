source("C:/Users/masse/Desktop/KefiR/KefiR-mbay/dev/mbay/scripts/relecture/agent_rev4_lib.R")
chk <- function(lab,y,g,u=1.5,r=.5){ st<-.bay_laws_stats(y,g); pr<-.bay_law_proposal(st,u); set.seed(1); th<-.bay_law_draw_def(pr,1500L); cc<-1/(st$k*th)
 for(s in c(list(NULL),as.list(seq_len(st$k)))){ rr<-if(is.null(s)) r else r*sqrt(2)
  a<-if(is.null(s)) .bay_law_gint(st,cc,rr) else .bay_law_gint(st,cc,rr,s=s)
  # sans pilote : par paquets de 50 lignes
  b<-unlist(lapply(split(seq_len(nrow(cc)),ceiling(seq_len(nrow(cc))/50)),function(i) if(is.null(s)) .bay_law_gint(st,cc[i,,drop=FALSE],rr) else .bay_law_gint(st,cc[i,,drop=FALSE],rr,s=s)))
  h<-min(.25,.5/sqrt(st$k))/4; f<-unlist(lapply(split(seq_len(nrow(cc)),ceiling(seq_len(nrow(cc))/50)),function(i) if(is.null(s)) .bay_law_gint(st,cc[i,,drop=FALSE],rr,h=h) else .bay_law_gint(st,cc[i,,drop=FALSE],rr,h=h,s=s)))
  lq<-.bay_law_lq_def(pr,th); ld<-lgamma(st$k*u)-st$k*lgamma(u)+(u-1)*rowSums(log(th)); lw<-drop(a)+ld-lq; wt<-exp(lw-max(lw));wt<-wt/sum(wt)
  cat(sprintf("%s s=%s pilote-sans: max %.2e, pondere %.2e | h/4: max %.2e, pondere %.2e\n",lab,format(s),max(abs(drop(a)-b)),sum(wt*abs(drop(a)-b)),max(abs(b-f)),sum(wt*abs(b-f))))
 };flush.console()}
set.seed(61); g<-rep(letters[1:3],each=6); chk("k=3 ordinaire",rnorm(18),g)
set.seed(62); g<-rep(letters[1:5],each=8); y<-rnorm(40); y[g=="e"]<-y[g=="e"]+1000; chk("k=5 d=1000",y,g)
set.seed(63); g<-rep(letters[1:3],c(400,400,5)); y<-rnorm(805); y[g=="c"]<-y[g=="c"]+8; chk("k=3 N=805, petit groupe loin",y,g)
set.seed(64); g<-rep(letters[1:4],each=300); y<-rnorm(1200,c(0,0,0,3)[as.integer(factor(g))],c(.01,1,1,10)[as.integer(factor(g))]); chk("k=4 N=1200 sd 0,01..10",y,g)

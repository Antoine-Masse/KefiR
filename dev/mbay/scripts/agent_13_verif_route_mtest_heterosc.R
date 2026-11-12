suppressMessages(pkgload::load_all("C:/Users/masse/Desktop/KefiR/KefiR-mbay", quiet=TRUE))
set.seed(4242); n<-c(5,5,40); s<-c(4,4,1); g<-factor(rep(1:3,n))
out<-sapply(1:30,function(i){y<-rnorm(sum(n),0,rep(s,n));
 r<-tryCatch(suppressWarnings(suppressMessages(m.test(y,g,verbose=FALSE,plot=FALSE,return=TRUE,boot=FALSE))),error=function(e)NULL)
 if(is.null(r)) return(c(NA,NA)); c(r$global_pvalue, r$method)})
cat("methodes:\n");print(table(out[2,],useNA="ifany"))
p<-as.numeric(out[1,]); cat("N=",sum(!is.na(p))," m.test global p<.05 :",mean(p<.05,na.rm=TRUE),"\n")

.libPaths(c("C:/Users/masse/Desktop/KefiR/mbay-lib", .libPaths()))
suppressMessages(library(BayesFactor))
set.seed(11)
# Croisement pur 2x2 : A inverse l'effet de B, moyenne de A nulle
d<-expand.grid(A=factor(c("a1","a2")),B=factor(c("b1","b2")),r=1:15)
d$y<-rnorm(nrow(d))+ifelse(d$A=="a1",1,-1)*ifelse(d$B=="b1",1,-1)
print(aggregate(y~A+B,d,mean)); print(anova(lm(y~A*B,d)))
m<-anovaBF(y~A*B,data=d,whichModels="withmain",progress=FALSE,iterations=50000)
lb<-c(m@bayesFactor$bf,0); names(lb)<-c(rownames(m@bayesFactor),"nul"); print(round(lb,2))
l<-function(x) exp(lb[[x]]-max(lb))
P<-function(bf) bf/(1+bf)
cat("apparie A   : BF=",l("A + B")/l("B")," P=",P(l("A + B")/l("B")),"\n")
cat("apparie A:B : BF=",l("A + B + A:B")/l("A + B")," P=",P(l("A + B + A:B")/l("A + B")),"\n")
all<-sapply(names(lb),l); inA<-grepl("A",names(lb))
cat("Uniforme marginal, A intervient (A ou A:B): P=",sum(all[inA])/sum(all)," a priori=",mean(inA),"\n")
# Question 'A joue-t-il un role, seul ou en interaction' : A+B+A:B contre B
bf<-l("A + B + A:B")/l("B"); cat("A+B+A:B contre B : BF=",bf," P(pi=1/2)=",P(bf),"\n")
# Identite annoncee : BF apparie = cotes conditionnelles de l'uniforme, sachant le sous-espace {B, A+B}
cat("Identite: uniforme restreint {B,A+B} P(A)=",l("A + B")/(l("A + B")+l("B")),"\n")

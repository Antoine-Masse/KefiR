#' Automatic mean/median comparison
#'
#' `mm.test()` performs all possible mean comparison tests, even if not mathematically correct and even if hypotheses are wrong. It executes tests to verify hypotheses, and the results are summarized in a graph. The user needs to choose the correct path, which can be done automatically by `m.test()`. By default, it finds the path that is more to the left.
#' @author Julien BOUSQUET (2021)
#'
#' @param x The values: a vector of values, with the samples to be compared.
#' @param group Categories: a vector of factors with fewer than maxcat categories.
#' @param pval The usual confidence level, 0.05 by default.
#' @param verbose `FALSE` by default. Increase information.
#' @param return Logical. Should the results be returned? Default is `TRUE`.
#' @param paired `FALSE` by default, but can be set to `TRUE` if each value of x in categories are paired and in the same order.
#' @param pval_ks The p-value for distribution tests like Kolmogorov-Smirnov or Shapiro: 0.01 by default.
#' @param maxcat The maximum number of categories: 50 by default.
#' @param plot `TRUE` to plot the graph with all the results.
#' @param silent `TRUE` to keep the function silent: no warnings.
#' @param boot `FALSE` by default. If `TRUE`, perform bootstrap sampling with `iter` iterations.
#' @param iter Number of iterations in the bootstrap, 100 by default.
#' @param conf Confidence level of the bootstrap, 0.95 by default.
#' @param code If `TRUE`, return the code that performs the different tests. `FALSE` by default.
#' @param debug. `FALSE` by default. If `TRUE`, print debugging information.
#' @return A list containing the results of the tests.
#' @import foreach
#' @import doParallel
#' @importFrom stats var.test t.test wilcox.test fligner.test aov oneway.test kruskal.test bartlett.test
#' @importFrom graphics arrows par plot.new points text
#' @importFrom grDevices dev.new
#' @examples
#' # Example 1
#' data(iris)
#' #mm.test(iris$Sepal.Length, iris$Species)
mm.test <- function(x, 
	group, 
	pval=0.05, 
	verbose=FALSE,
	return=TRUE,
	paired=FALSE,
	pval_ks=0.01,
	maxcat=50,
	plot=TRUE,
	silent=TRUE,
	boot=FALSE,
	iter=100,
	conf=0.95,
	code=FALSE,
	debug.=FALSE){
  arr.len <- 0.2 #length of arrows
  arr.lwd <- 2 #line width of arrows
  arr.angle <- 20 #angles of arrows
  arr.vspace <- 0.5 #vertival space between arrow and aim
  vspace <- 5#vertical spaces between strata representation
  # Validate Number of categories
  C <- length(unique(group)) #number of categories
  if(debug.)
  if(C>maxcat){stop("Too much categories : is 'group' really a  factor?\n See parameter maxcat=\n ")}
  if(C>0.2*length(group)){warning("Too much categories in group")}
  if(C<2){stop("Less than 2 categories. Not enough")}
  CAT <- unique(group) 
  L <- list() 
  L$category.number <- C
  L$category.name <- CAT


  if(debug.)cat("Strate 0 : sample size \n")
  foreach(i = CAT, .combine=rbind)%dopar%{
     sum( group==i )
  } -> SamplesSizes # vector of sizes of samples by category
L$category.size = data.frame(category=CAT,  size=SamplesSizes)

  r <- range(SamplesSizes)
  Nmin <- r[1] # min size of sample 
  Nmax <- r[2] # max size of sample
  if(Nmin ==0) stop("Some category is empty") 
  if(Nmin <4) stop("Some category as less than 3 elements")
  if(Nmin < 10) warning("Some category as less than 10 elements")

 #Strate 1 : repartition
   # Normality inside  categories 
   sidak <- 1-(1-pval)^(1/C) # pvalue adjusted by Sidak method
   foreach(i = CAT, .packages='normtest', .combine=rbind)%dopar%{
	#ind is the vector of indexes of the current category
	ind <-  group==i
	#N.temp is the size of the current category
	N.temp <- L$category.size$size[L$category.size$category==i]
                if(N.temp <100){
                   data.frame(category=i, test='Shapiro', sidak=sidak, p.value=round(shapiro.test(x[ind])$p.value,4))
                }else if(N.temp<1000){
                   data.frame(category=i, test='Jarque Bera', sidak=sidak, p.value=round(jb.norm.test(x[ind])$p.value,4))
                }else{
                    data.frame(category=i, test='Big number', sidak=sidak, p.value=1)
	}
   } -> L$repartition
   if(debug.){cat('Shapiro-Wilk or Jarque Bera p.values : \n');  print(L$repartition)}

#Only two categories!!#######
# Strate 2 : test of variance for 2 categories
   if(C==2){
	#ind is the vector of indexes of the current category
	ind <-  group==CAT[1]
	data.frame(p.value=var.test(x[ind], x[!ind])$p.value) -> L$var.test
	data.frame(p.value=lawstat::levene.test(x, group)$p.value) -> L$levene.test
    	data.frame(p.value=fligner.test(x, group)$p.value) -> L$fligner.test

   if(debug.){cat('Fisher test with 2 samples p.values : \n'); print(L$var.test)}
# Strate 3 : tests of mean for 2 categories
                #var.equal=TRUE
	#ind is the vector of indexes of the current category
	ind <-  group==CAT[1]
	data.frame(p.value=t.test(x[ind], x[!ind], var.equal=TRUE, paired=paired)$p.value) -> L$t.test.varequal
	data.frame(p.value=t.test(x[ind], x[!ind], var.equal=FALSE, paired=paired)$p.value) -> L$t.test.vardiff
	data.frame(p.value=wilcox.test(x[ind], x[!ind])$p.value) -> L$wilcoxon.mann.whitney
	if(paired)data.frame(p.value=wilcox.test(x[ind], x[!ind], paired=TRUE)$p.value) -> L$wilcoxon
   	if(debug.){cat('Student test with 2 samples  -- varequal = T p.values : \n'); print(L$t.test.varequal)}
	# Strate 4 : No POstHoc with 2 ech.
# Strate 5 : representation for 2 categories
	# Start new device with appropriate size
	dev.new(width=10, height=7, unit='in')
	graphical.Height <- 70
	graphical.Width <- 100

	plot(NA, axes=F, xlab='', ylab='',#empty plot
	  xlim=c(0, graphical.Width), ylim=c(graphical.Height,0)
	)

	pos0 <- device.print(0, 0,  title= "----Normality of repartition-------------------------")
	pos1 <- device.print(50, pos0[2],  L$repartition)  #, "Normality of repartition")
	
	pos2 <- device.print(0, pos1[2]+vspace, title="----Homogeneity of variance -------------------------")
	pos3 <- device.print(12, pos2[2], info=L$var.test, title="Fisher var.test()")
	pos4 <- device.print(50, pos2[2], L$levene.test, "Levene (Brown Forsythe)")

	# only for paired data, and do wilcox rank sign test
	if(paired){
	  # Wilcoxon rank sign need paired values
	  pos5 <- device.print(88, pos2[2], L$fligner.test, "Fligner test (most robust)")
	}#else{
	  #Mann Witney don't need hypothesis on equality of variance
	  #Nothing	
	  #}

	pos6 <- device.print(0, pos4[2]+vspace, title="----Equality of means -----------------------------")
	pos7 <- device.print(12, pos6[2], L$t.test.varequal, "Student  with equal var")

	if(paired){x.correction <- -5}else{ x.correction <- 0}
	pos8 <- device.print(50+x.correction, pos6[2], L$t.test.vardiff, "Student with var diff")

 	if(paired){
	  pos9 <- device.print(75, pos6[2], L$wilcoxon, "Wilcoxon rank")
	  pos10 <- device.print(95,pos6[2], L$wilcoxon.mann.whitney, "Mann-Whitney")
	}else{
	  pos10 <- device.print(88,pos6[2], L$wilcoxon.mann.whitney, "Mann-Whitney")
	}

	# arrow from normality to Fisher var.test()
	col0.1 <- ifelse(min(L$repartition$p.value)>0.05, 'green', 'grey')
	arrows(pos1[1], pos1[2], pos3[1], pos2[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.1) #arc0.1

	# arrow from normality to levene
	col0.2 <- ifelse(min(L$repartition$p.value)>L$repartition$sidak[1] & min(L$repartition$p.value)<0.05 , 'green', 'grey')
	arrows(pos1[1], pos1[2], pos4[1],pos2[2]-arr.vspace , code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.2) #arc0.2

	if(paired){# arrow from normality to fligner Killeen
	   col0.3 <- ifelse(min(L$repartition$p.value)<L$repartition$sidak[1], 'green', 'grey')
	   arrows(pos1[1], pos1[2], pos5[1], pos2[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.3) #arc0.3
	}

	#arrow from var.test() to Student equal var
	col1.1 <- ifelse(L$var.test$p.value>pval  , 'green', 'grey')
	arrows(pos3[1],pos3[2],pos7[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.1) #arc1.1

	#arrow from var.test to student with diff var
	col1.2 <- ifelse(L$var.test$p.value<pval  , 'green', 'grey')
	arrows(pos3[1],pos3[2],pos8[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.2) #arc1.2

	#arrow from levene-Brown Forsythe to Student with equal means
	col1.3 <- ifelse(L$levene.test$p.value>pval  , 'green', 'grey')
	arrows(pos4[1],pos4[2],pos7[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.3) #arc1.3

	#arrow from levene-Brown Forsythe to Student with diff means
	col1.4 <- ifelse(L$levene.test$p.value<pval  , 'green', 'grey')
	arrows(pos4[1],pos4[2],pos8[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.4) #arc1.4

	if(!paired){
	  #arrow from normality to Mann-Whitney
	  col1.4 <- ifelse(min(L$repartition$p.value)<sidak  , 'green', 'grey')
	  arrows(pos1[1], pos1[2],pos10[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.4) #arc1.4
	}else{
	  #arrow from normality to fligner killeen
 	  #Done in arrow 0.3

	  #arrow from fligner to wilcoxon
	  col1.5 <- ifelse(min(L$fligner.test$p.value)>pval  , 'green', 'grey')
	  arrows(pos5[1], pos5[2],pos9[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.5) #arc1.4

	  #arrow from fligner to Mann Whitney
	  col1.6 <- ifelse(min(L$fligner.test$p.value)<=pval  , 'green', 'grey')
	  arrows(pos5[1], pos5[2], pos10[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.6) #arc1.4
	}
   } 

#More than 2 categories#######

   if(C>2){
# Strate 2 : test of variance for 3 categories or more
	#ind is the vector of indexes of the current category
	data.frame(p.value=bartlett.test(x, group)$p.value)-> L$bartlett.test
	data.frame(p.value=lawstat::levene.test(x, group)$p.value) -> L$levene.test
    	data.frame(p.value=fligner.test(x, group)$p.value) -> L$fligner.test

   if(debug.){cat('Strate 2 : test of variance for 3 or more groups ended.\n')}
# Strate 3 : tests of mean for 3 categories or more
	data.frame(p.value=oneway.test(x~group, var.equal=FALSE)$p.value) -> L$oneway.test.vardiff
	data.frame(p.value=oneway.test(x~group,  var.equal=TRUE)$p.value) -> L$anova
	data.frame(p.value=kruskal.test(x, group)$p.value) -> L$kruskal.test
	WRS2::t1way(x~group) -> T1WAY
	data.frame(p.value=T1WAY$p.value) -> L$t1way
	data.frame(p.value=WRS2::med1way(x~group)$p.value) -> L$med1way
   	if(debug.){cat('End of STrata 3 \n'); }

# Strate 4 : Post Hoc 
	as.data.frame(catego(pairwise.t.test(x, group, sd.pool=FALSE) )$groups) -> L$pairwise.t.test.catego
	AOV <- aov(x~group) #store ANOVA result for Tukey and SNK
	SNK <- agricolae::SNK.test(AOV, trt="group", group=TRUE, alpha=pval)$groups 
	data.frame(categories=rownames(SNK), groups=SNK[,2]) -> L$SNK.test
	#TukeyHSD(AOV, ordered=T, conf.level=1-pval) -> L$TukeyHSD
	as.data.frame(catego(pairwise.wilcox.test(x, group))$groups) -> L$pairwise.wilcox.test.catego
	as.data.frame(lincon(x~group)$comp[,c(1:2,6)]) -> L$lincon

# Strate 5 : representation for 3 or more categories
	# Start new device with appropriate size

	graphical.Height <- 70 +2.5*C
	graphical.Width <- 100
	dev.new(width= graphical.Width/10, height= graphical.Height/10, unit='in')
	par(mar=c(0,0,0,0))
	plot(NA, axes=F, xlab='', ylab='',#empty plot
	  xlim=c(0, graphical.Width), ylim=c(graphical.Height,0)
	)

	pos0 <- device.print(0, 0,  title= "----Normality of repartition-------------------------")
	pos1 <- device.print(50, pos0[2],  L$repartition)  
	if(debug.)points(pos1[1], pos1[2], col='red', cex=2)

	pos2 <- device.print(0, pos1[2]+vspace, title="----Homogeneity of variance -------------------------")
	pos3 <- device.print(12, pos2[2], info=L$bartlett.test, title="Bartlett test")
	pos4 <- device.print(50, pos2[2], L$levene.test, "Levene (Brown Forsythe)")
	pos5 <- device.print(88, pos2[2], L$fligner.test, "Fligner Killeen")

	pos6 <- device.print(0, pos4[2]+vspace, title="----Equality of means -----------------------------")
	pos7 <- device.print(5, pos6[2], L$oneway.test.vardiff, "ANOVA var diff")
	pos8 <- device.print(27.5, pos6[2], L$anova, "ANOVA same var")
	pos9 <- device.print(50, pos6[2], L$kruskal.test, "Kruskal")
	pos10 <- device.print(72.5,pos6[2],L$t1way , "AOV trim 20%")
	pos11 <- device.print(93,pos6[2],L$med1way , "AOV median")

	pos12 <- device.print(0, pos11[2]+vspace, title="----Grouping methods -------------------------")
	pos13 <- device.print(pos7[1]+5, pos12[2], L$pairwise.t.test.catego, title="Pair.t.test var diff")
	pos14 <- device.print(pos8[1]+9, pos12[2], info=L$SNK.test, title="SNK")
	#pos15 <- device.print(pos8[1]+15, pos12[2], info=L$TukeyHSD, title="Tukey")
	pos16 <- device.print(pos9[1]+14, pos12[2], L$pairwise.wilcox.test.catego, title="Pairwise.wilcoxon")
	pos17 <- device.print(pos11[1]-5, pos12[2], L$lincon, title='lincon')


	# arrow from normality to bartlett.test()
	col0.1 <- ifelse(min(L$repartition$p.value)>0.05, 'green', 'grey')
	arrows(pos1[1], pos1[2], pos3[1], pos2[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.1) #arc0.1

	# arrow from normality to levene
	col0.2 <- ifelse(min(L$repartition$p.value)>L$repartition$sidak[1] & min(L$repartition$p.value)<pval , 'green', 'grey')
	arrows(pos1[1], pos1[2], pos4[1],pos2[2]- arr.vspace , code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.2) #arc0.2
	#arrow from normality to fligner
	col0.3 <- ifelse(min(L$repartition$p.value)<L$repartition$sidak[1], 'green', 'grey')
	arrows(pos1[1], pos1[2], pos5[1], pos2[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col0.3) #arc0.3
	
	#arrow from bartlett.test to oneway
	col1.1 <- ifelse(L$bartlett.test$p.value>pval  , 'green', 'grey')
	arrows(pos3[1],pos3[2],pos7[1],pos6[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.1) #arc1.1

	#arrow from bartlett.test to anova
	col1.2 <- ifelse(L$bartlett.test$p.value<=pval  , 'green', 'grey')
	arrows(pos3[1],pos3[2],pos8[1],pos6[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.2) #arc1.2

	#arrow from levene-Brown Forsythe to oneway with diff means
	col1.3 <- ifelse(L$levene.test$p.value>pval  , 'green', 'grey')
	arrows(pos4[1],pos4[2],pos7[1],pos6[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.3) #arc1.3


	#arrow from levene-Brown Forsythe to ANOVA with equal means
	col1.4 <- ifelse(L$levene.test$p.value<pval  , 'green', 'grey')
	arrows(pos4[1],pos4[2],pos8[1],pos6[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.4) #arc1.4

	#arrow from fligner to AOV Kruskal
	col1.5 <- ifelse(min(L$fligner.test$p.value)>pval  , 'green', 'grey')
	arrows(pos5[1], pos5[2],pos9[1],pos6[2]- arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.5) #arc1.5

	#arrow from fligner to AOV trim
	col1.6 <- ifelse(min(L$fligner.test$p.value)<=pval  , 'green', 'grey')
	arrows(pos5[1], pos5[2],pos10[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.6) #arc1.4

	#arrow from fligner to AOV Median
	col1.7 <- ifelse(min(L$fligner.test$p.value)<=pval  , 'green', 'grey')
	arrows(pos5[1], pos5[2], pos11[1],pos6[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col1.7) #arc1.4

	#arrow from oneway to pairwise.t.test()
	col2.1 <- ifelse(min(L$oneway.test.vardiff$p.value)<pval  , 'green', 'red')
	arrows(pos7[1], pos7[2], pos13[1], pos12[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col2.1) #arc2.1


	#arrow from oneway to pairwise.t.test()
	col2.2 <- ifelse(min(L$anova$p.value)<pval  , 'green', 'red')
	arrows(pos8[1], pos8[2], pos14[1], pos12[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col2.2) #arc2.2

	#arrow from kruskal to pairwise.wilcoxon.test()
	col2.3 <- ifelse(min(L$kruskal.test$p.value)<pval  , 'green', 'red')
	arrows(pos9[1], pos9[2], pos16[1], pos12[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col2.3) #arc2.3

	#arrow from AOV trim to lincon
	col2.4 <- ifelse(min(L$t1way$p.value)<pval  , 'green', 'red')
	arrows(pos10[1], pos10[2], pos17[1], pos12[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col2.4) #arc2.4

	#arrow from AOV median to lincon
	col2.5 <- ifelse(min(L$med1way$p.value)<pval  , 'green', 'red')
	arrows(pos11[1], pos11[2], pos17[1], pos12[2]-arr.vspace, code=2, length=arr.len, lwd=arr.lwd, angle=arr.angle, col=col2.5) #arc2.5
  } 
return(L)
}

 
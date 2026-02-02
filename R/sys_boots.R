.boots <- function(x,g,ctrl=FALSE,type="mean",var.equal=FALSE,
	conf=0.95,iter=500,alpha=alpha,paired=paired, control=NULL) {
	pvals <- c()
	which(unique(g)==control)-> ind_control
	for (i in 1:iter) {
	  if (type=="mean") {
		  temp <- t.test(sample(x[g==unique(g)[1]],replace=TRUE),sample(x[g==unique(g)[2]],replace=TRUE),var.equal=var.equal,paired=paired)$p.value
	  } else if (type=="median") {
		#temp <- wilcox.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE))$p.value
		  temp <- wilcox.test(sample(x[g==unique(g)[1]],replace=TRUE),sample(x[g==unique(g)[2]],replace=TRUE),exact=FALSE,paired=paired)$p.value
	  } else if (type=="ks") {
		ech1 <- sample(x[g==unique(g)[1]],replace=TRUE)
		ech2 <- sample(x[g==unique(g)[2]],replace=TRUE)
		ech1 <- (ech1 - median(ech1))/sd(ech1)
		ech2 <- (ech2 - median(ech2))/sd(ech2)
		temp <- ks.test(ech1,ech2)$p.value
	  }
	  pvals <- c(pvals,temp)
	}
	pvals <- quantile(pvals,probs=conf,na.rm=T)
	if (ctrl == TRUE) {
	  if (pvals <= alpha) {
		synth <- list()
		starss <- c("","")
		starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",
																  ifelse(pvals <=0.05,"*","")))
		synth$groups <- data.frame(categories=unique(g),group=starss)
		synth$p.value <- pvals
	  } else {
		synth <- list()
		synth$groups <- data.frame(categories=unique(g),group=c("",""))
		synth$p.value <- pvals
	  }
	} else if (ctrl==FALSE) {
	  synth <- list()
	  if (pvals <= alpha) {
		synth$groups <- data.frame(categories=unique(g),groups=c("a","b"))
		synth$p.value <- pvals
	  } else {
		synth$groups <- data.frame(categories=unique(g),groups=c("a","a"))
		synth$p.value <- pvals
	  }
	}
	return(synth)
  }

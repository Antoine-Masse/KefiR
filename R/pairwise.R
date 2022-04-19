
#' Automation function for pairwise calculations.
#'
#' @param x numerical vector
#' @param g category vector
#' @param type 'mean' for pairwise.t.test(p.adjust.method="holm"), 'median' for pairwise.wilcox.test(p.adjust.method="BH"), 'ks' for ks.test(), 'lincon' for lincon() of {WSR2}
#' @param alpha threshold value of p-value to establish the groups.
#' @param control name of the category that will be used as a control to establish differences with '*', '**' and '***'.
#' @param pool.sd switch to allow/disallow the use of a pooled SD.
#' @param silent for displaying or not warnings.
#' @param boot to activate the boostrap on 'mean' and 'median'.
#' @param iter number f iterations (boot==TRUE).
#' @param conf confidence level of bootstrap.
#'
#' @return This function automates the work of the ks.test(), lincon() functions of {WSR2}, pairwise.t.test() and pairwise.wilcox.test() and extracts groups of means or comparisons to a control with the catego() function.
#' @return It pre-sorts the means/medians to ensure that the groups are identified in ascending order.
#' @return It also identifies the robustness of these groups by establishing a bootstrap.
#' @importFrom WRS2 lincon
#' @export
#'
#' @examples
#' data(iris)
#' pairwise(iris[,1],iris[,5],type="mean")# t.test
#' pairwise(iris[,1],iris[,5],type="median",alpha=0.01,boot=TRUE)#wilcox
#' pairwise(iris[,1],iris[,5],type="ks")
#' pairwise(iris[,1],iris[,5],type="lincon")
pairwise <- function(x,g,type="mean",alpha=0.05,control=c(),pool.sd=FALSE,silent=TRUE,boot=FALSE,iter=500,conf=0.95) {
  g <- factor(g)
  init_order <- levels(g)
  which(levels(g)%in%control)-> ind_control
	if (type=="mean") {
		if (length(ind_control)==1) {
			categories <- c(levels(g)[ind_control],levels(g)[-ind_control])
			g <- ordered(g, levels = categories)
			pairwise.t.test(x,g,pool.sd=pool.sd,p.adjust.method="holm")-> result
			groups <- catego(result,control=control)
			if (boot==TRUE) {
				mymat <- groups$p.value
				mydim <- dim(mymat)
				mymat <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
				#barre <- txtProgressBar(min=0,max=iter,width=50)
				#cat("\nBootstrap mean control\n")
				#cat(paste(rep("=",50),collapse=""),"\n")
				for (i in 1: iter) {
				  #setTxtProgressBar(barre,i)
				  x_temp <- x
				  for (j in levels(g)) {

					x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
				  }
				  temp <- pairwise.t.test(x_temp,g,pool.sd=pool.sd,p.adjust.method="holm")$p.value
				  mymat[i,,] <- temp
				}
				#close(barre)
				output <- list()
				apply(mymat,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
				colnames(output$p.value) <- colnames(groups$p.value)
				rownames(output$p.value) <- rownames(groups$p.value)
				groups$bootstrap <- catego(output,control=control)
			  }
		} else {
		  mu <- by(x,g,mean)
		  categories <- levels(g)
		  indices <- order(mu)
		  g <- ordered(g, levels = categories[indices])
		  pairwise.t.test(x,g,pool.sd=pool.sd,p.adjust.method="holm")-> result
		  groups <- catego(result,alpha=alpha)
		  if (boot==TRUE) {
			mymat <- groups$p.value
			mydim <- dim(mymat)
			mymat <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
			#barre <- txtProgressBar(min=0,max=iter,width=50)
			#cat("\nBootstrap mean\n")
			#cat(paste(rep("=",50),collapse=""),"\n")
			for (i in 1: iter) {
			  #setTxtProgressBar(barre,i)
			  x_temp <- x
			  for (j in levels(g)) {
				x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
			  }
			  temp <- pairwise.t.test(x_temp,g,pool.sd=pool.sd,p.adjust.method="holm")$p.value
			  mymat[i,,] <- temp
			}
			#close(barre)
			output <- list()
			apply(mymat,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
			colnames(output$p.value) <- colnames(groups$p.value)
			rownames(output$p.value) <- rownames(groups$p.value)
			groups$bootstrap <- catego(output,alpha=alpha)
		  }
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		return(groups)
	} else if (type == "median" ){
		if (length(ind_control)==1) {
			categories <- c(levels(g)[ind_control],levels(g)[-ind_control])
			g <- ordered(g, levels = categories)
			pairwise.wilcox.test(x,g,p.adjust.method="BH")-> result
			groups <- catego(result,control=control)
			if (boot==TRUE) {
				#print(groups$p.value)
				mymat <- groups$p.value
				mydim <- dim(mymat)
				mymat <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
				#barre <- txtProgressBar(min=0,max=iter,width=50)
				#cat("\nBootstrap median control\n")
				#cat(paste(rep("=",50),collapse=""),"\n")
				for (i in 1: iter) {
				  #setTxtProgressBar(barre,i)
				  x_temp <- x
				  for (j in levels(g)) {

					x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
				  }
				  temp <- pairwise.wilcox.test(x_temp,g,p.adjust.method="BH")$p.value
				  mymat[i,,] <- temp
				}
				#close(barre)
				output <- list()
				apply(mymat,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
				colnames(output$p.value) <- colnames(groups$p.value)
				rownames(output$p.value) <- rownames(groups$p.value)
				groups$bootstrap <- catego(output,control=control)
			}
		} else {
		  mu <- by(x,g,median)
		  #print(mu)
		  categories <- levels(g)
		  indices <- order(mu)
		  g <- ordered(g, levels = categories[indices])
		  pairwise.wilcox.test(x,g,p.adjust.method="BH")-> result
		  groups <- catego(result,alpha=alpha)
		  if (boot==TRUE) {
			#print(groups$p.value)
			mymat <- groups$p.value
			mydim <- dim(mymat)
			mymat <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
			#barre <- txtProgressBar(min=0,max=iter,width=50)
			#cat("\nBootstrap median\n")
			#cat(paste(rep("=",50),collapse=""),"\n")
			for (i in 1: iter) {
			  #setTxtProgressBar(barre,i)
			  x_temp <- x
			  for (j in levels(g)) {
				x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
			  }
			  temp <- pairwise.wilcox.test(x_temp,g,p.adjust.method="BH")$p.value
			  mymat[i,,] <- temp
			}
			#close(barre)
			output <- list()
			apply(mymat,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
			colnames(output$p.value) <- colnames(groups$p.value)
			rownames(output$p.value) <- rownames(groups$p.value)
			groups$bootstrap <- catego(output,alpha=alpha)
		  }
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
		return(groups)
	} else if (type == "ks" ){
		unique_g <- unique(g)
		#print(g)
		for (i in unique_g) {
			x[g==i] <- (x[g==i]-mean(x[g==i]))/sd(x[g==i])
		}
		#return(data.frame(x,g))
		mymat <- matrix(rep(NA,(length(unique_g)-1)^2),nc=(length(unique_g)-1),nr=(length(unique_g)-1))
		rownames(mymat) <- unique_g[2:length(unique_g)] ; colnames(mymat) <- unique_g[1:(length(unique_g)-1)]
		#print(mymat)
		ks_func <- function(x,g,mymat,unique_g) {
			#print(unique_g)
			#print(unique_g[1])
			for (i in 1:(length(unique_g)-1)) {
				for (j in (i+1):length(unique_g)) {
					#cat("i",unique_g[i], " et ")
					#cat("j",unique_g[j]," pour ")
					#print(x[g==unique_g[i]])
					#print(x[g==unique_g[j]])
					pv <- ks.test(x[g==unique_g[i]],x[g==unique_g[j]])$p.value
					#cat("pv",pv,"\n")
					mymat[(j-1),i] <- pv
				}
			}
			return(mymat)
		}
		if (boot == FALSE) {mymat <- suppressWarnings(ks_func(x,g,mymat,unique_g)) ; output <- list() ; output$p.value <- mymat
		} else if (boot==TRUE) {
		  mydim <- dim(mymat)
		  myarray <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
		  for (k in 1: iter) {
			x_temp <- x
			for (j in levels(g)) {
			  x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
			}
			temp <- suppressWarnings(ks_func(x_temp,g,mymat,unique_g))
			myarray[k,,] <- temp
		  }
		  output <- list()
		  apply(myarray,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
		  colnames(output$p.value) <- colnames(mymat)
		  rownames(output$p.value) <- rownames(mymat)
		}
		return(output)
	} else if (type == "lincon" ){
		if (length(ind_control)==1) {
		  categories <- c(levels(g)[ind_control],levels(g)[-ind_control])
		  g <- ordered(g, levels = categories)
		  t1 <- lincon(x~g) #t1way
		  mymat <- matrix(rep(NA,(length(t1$fnames)-1)^2),nc=length(t1$fnames)-1,nr=length(t1$fnames)-1)
		  rownames(mymat) <- t1$fnames[2:length(t1$fnames)] ; colnames(mymat) <- t1$fnames[1:length(t1$fnames)-1]
		  for (i in 1:nrow(t1$comp)) {
			mymat[(t1$comp[i,2]-1),t1$comp[i,1]] <- t1$comp[i,6]
		  }
		  result <- list() ; result$p.value <- mymat
		  groups <- catego(result,control=control)
		} else {
		  mu <- by(x,g,median)
		  #print(mu)
		  categories <- levels(g)
		  indices <- order(mu)
		  g <- ordered(g, levels = categories[indices])
		  t1 <- lincon(x~g) #t1way
		  mymat <- matrix(rep(NA,(length(t1$fnames)-1)^2),nc=length(t1$fnames)-1,nr=length(t1$fnames)-1)
		  rownames(mymat) <- t1$fnames[2:length(t1$fnames)] ; colnames(mymat) <- t1$fnames[1:length(t1$fnames)-1]
		  for (i in 1:nrow(t1$comp)) {
			mymat[(t1$comp[i,2]-1),t1$comp[i,1]] <- t1$comp[i,6]
		  }
		  result <- list() ; result$p.value <- mymat
		  groups <- catego(result,alpha=alpha)
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		return(groups)
	} else if (type == "boot" ){
		if (length(ind_control)==1) { # control = TRUE
			categories <- c(levels(g)[ind_control],levels(g)[-ind_control])
			g <- ordered(g, levels = categories)
			pairwise.boot(x,g,iter=iter)-> result
			groups <- catego(result,control=control)
		} else {
			mu <- by(x,g,meanbp)
			categories <- levels(g)
			indices <- order(mu)
			g <- ordered(g, levels = categories[indices])
			pairwise.boot(x,g,iter=iter)-> result
			groups <- catego(result,alpha=alpha)
		}
		match(init_order,levels(g))-> indices
		groups$groups <- groups$groups[indices,]
		return(groups)
	} else {stop("Error! arg type is 'mean','median', 'ks' or 'lincon'.\n")}
}

#' Calculate an pairwise difference between samples by bootstrap without encountering a central-limit problem
#'
#' @param x numerical vector
#' @param g categorical vector
#' @param iter number of iteration
#' @param mu comparison criterion, by default 'meanbp' (moving average per iteration: an automatic compromise between mean and median), otherwise 'mean', 'median' or 'sd'.
#'
#' @return This function returns an array of p-values like the other peerwise functions.
#' @return It calculates the number of differences between bootstrap samples that do not include 0 (all higher or all lower).
#' @return This number of differences divided by the number of iterations (iter) gives the maximum percentage of times that a convergent difference (all higher or all lower) is found: a confidence.
#' @return This value subtracted from 1 gives an equivalent of the p-value whose precision depends on the number of iterations.
#' @return Note: using the meanbp criterion is more relevant, it allows a compromise between mean and median by avoiding leverage effects.
#' @export
#'
#' @examples
#' # Example 1
#' data(iris)
#' pairwise.boot(iris[,2],iris$Species)
#' # Example 2 by using pairwise(type=="boot")
#' data(mtcars)
#' pairwise(mtcars$mpg[mtcars$carb<=4],mtcars$carb[mtcars$carb<=4],type="boot")
pairwise.boot <- function(x,g,iter=500,mu="meanbp") {
  g <- as.factor(g)
	unique_g <- levels(g)
	mymat <- matrix(rep(NA,(length(unique_g)-1)^2),nc=(length(unique_g)-1),nr=(length(unique_g)-1))
	rownames(mymat) <- unique_g[2:length(unique_g)] ; colnames(mymat) <- unique_g[1:(length(unique_g)-1)]
	for (i in 1:(length(unique_g)-1)) {
		for (j in (i+1):length(unique_g)) {
			ech1 <- x[g==unique_g[i]]
			ech2 <- x[g==unique_g[j]]
			diff <- c()
			for (k in 1:iter) {
				ech1_temp <- sample(ech1,replace=TRUE)
				ech2_temp <- sample(ech2,replace=TRUE)
				if (mu=="meanbp") {diff <- c(diff, (meanbp(ech1_temp)-meanbp(ech2_temp)))}
				else if (mu=="mean") {diff <- c(diff, (mean(ech1_temp)-mean(ech2_temp)))}
				else if (mu=="median") {diff <- c(diff, (median(ech1_temp)-median(ech2_temp)))}
				else if (mu=="sd") {diff <- c(diff, (sd(ech1_temp)-sd(ech2_temp)))}
			}
			length(diff[diff>0]) -> sup
			length(diff[diff<0]) -> inf
			if (sup>inf) {
				pv <- 1-(sup/iter)
			}else {pv <- 1-(inf/iter)}
			mymat[(j-1),i] <- pv*2 # *2 because bilateral
		}
	}
	result <- list() ; result$p.value <- mymat
	return(result)
}




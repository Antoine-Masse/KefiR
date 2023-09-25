#' Determining groups based on the results of a pairwise function
#'
#' @param result output results of the pairwise(), pairwise.t.test() or pairwise.wilcox.test() function
#' @param control name of the category that will be used as a control to establish differences with '*', '**' and '***'.
#' @param alpha threshold of p-value to establish groups of averages or medians of type "a", "ab", "b" ...
#'
#' @importFrom stringr str_sort
#'
#' @return catego() can be used to establish groups with significantly different mean/median values relative to each other or to a reference control. It exploits the results of the functions pairwise.t.test(), pairwise.wilcox.test() and pairwise(type = "median" or "mean").
#' @export
#'
#' @examples
#' data(iris)
#' output <- pairwise.t.test(iris[,1],iris$Species)
#' catego(output)
#' catego(output,control="setosa")
catego <- function(result,control=c(),alpha=0.05) {
  result$p.value -> mymat
  categories <- c(colnames(mymat),rownames(mymat)[nrow(mymat)])
  #print("a") ; print(categories)
  #print("b") ; print(control)
  # ERREUR A CORRIGER : no match.
  which(categories%in%control)-> ind_control
  #which(control%in%levels(categories))-> ind_control
  if (length(ind_control)==1) {
    c(1,mymat[,1]) -> cats
    ifelse(cats<=0.001,"***",ifelse(cats<=0.01,"**",ifelse(cats<=0.05,"*","")))->cats
  } else {
    possibility <- letters[seq( from = 1, to = length(categories))]
    start = 0 ; pos=1 ; cats <- c("a") ; change <- 0
    for (i in 1:ncol(mymat)) {
      #cat("Passage",start,"\n")
      for (j in i:nrow(mymat)) {
        pvals <- mymat[j,i]
        #print(pvals)
        if (start==0) {
          if (pvals > alpha){
            #print("non-signif")
            cats <- c(cats,possibility[pos])
          } else {
            #print("signif")
            cats <- c(cats,possibility[pos+1])
            change <- 1
          }
        } else {
          if (pvals <= alpha){
            #print(length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]])))
            if (length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]]))>0){
              #print("signif avec changement de cat.")
              # a on avait b vs b, l'autre b devient c
              cats[j+1]<- possibility[pos+1]
              change <- 1

            } else {
              #print("signif")
              next
            }
          } else {
            if (length(intersect(strsplit(cats[i],"")[[1]],strsplit(cats[j+1],"")[[1]]))>0) {
              #print("non-signif")
              next
            } else {
              #print("signif mais sensé être de cat différente (situation intermédiaire)")
			  if (str_sort(c(cats[i],cats[j+1]))[1] == cats[i]) {
				cats[i] <- paste0(str_sort(c(cats[i],cats[j+1])),collapse="")
			  } else {
				cats[i+1] <- paste0(str_sort(c(cats[i],cats[j+1])),collapse="")
			  }
            }
          }
        }
      }
      #print(cats)
      #print("########")
      start <- start+1
      if (change == 1) {pos <- pos+1}
    }
  }
  synth <- list()
  groups <- cbind(categories,groups=cats) ; rownames(groups) <- rep("",nrow(groups))
  synth$groups <- groups
  synth$p.value <- result$p.value
  return(synth)
}

#' Automatic average/median/variance comparison function.
#'
#' @param data numerical vector
#' @param cat category vector
#' @param alpha p-value threshold value for all the tests.
#' @param verbose to display the full reasoning of the analysis.
#' @param return allows to return the results of pairwise analysis (p-values and groups).
#' @param paired (under development) to allow the analysis of matched data.
#' @param control name of the category that will eventually be used as a control.
#' @param maxcat maximum number of categories allowed. When this number is high, some tests may return an error message.
#' @param plot to display the distribution of the data.
#' @param silent for displaying or not warnings.
#' @param boot to activate the boostrap on 'mean' and 'median'.
#' @param iter number f iterations (boot==TRUE).
#' @param conf confidence level of bootstrap.
#' @param code allows to display the simplified R source code to be able to do the same R study step by step.
#' @param debug when m.test return error.
#'
#' @return m.test() runs a decision tree to choose the most appropriate test series for sample comparison.
#' @return She chooses the tests, justifies her choices.
#' @return It can output groups of means or a comparison to a control.
#' @return Finally, it will measure the robustness of the results by bootstrap.
#' @importFrom fda.usc fanova.hetero
#' @importFrom agricolae kurtosis
#' @importFrom agricolae skewness
#' @importFrom agricolae SNK.test
#' @importFrom lawstat levene.test
#' @importFrom WRS2 med1way
#' @importFrom WRS2 t1way
#' @importFrom WRS2 lincon
#' @importFrom onewaytests bf.test
#' @importFrom vioplot vioplot
#' @importFrom DescTools DunnettTest
#' @import methods
#' @export
#'
#' @examples
#' data(iris)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[1:100,1],iris[1:100,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,4],iris[,5],verbose=TRUE, plot=FALSE, return=FALSE, boot=FALSE)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE,control="virginica")
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,4],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
m.test <- function (data, cat, alpha=0.05, verbose=TRUE, return=TRUE, paired=FALSE,control=c(),
                    maxcat=50, plot=TRUE,silent=TRUE,boot=TRUE,iter=500,conf=0.95,
                    code=FALSE,debug=FALSE){
  if (code==TRUE) {verbose <- FALSE}
  if (debug==TRUE) {verbose <- FALSE ; code=FALSE}
  if (debug==TRUE) {print("Compilation des fonctions de base")}
  discret.test <- function(vector) {
    return(length(unique(vector))/length(vector))
  }
  boots <- function(data,cat,ctrl=FALSE,type="mean",var.equal=FALSE,conf=0.95,iter=500,alpha=alpha) {
    pvals <- c()
    for (i in 1:iter) {
      if (type=="mean") {
        temp <- t.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE),var.equal=var.equal)$p.value
      } else if (type=="median") {
		#temp <- wilcox.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE))$p.value
		temp <- wilcox.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE),exact=FALSE)$p.value
      } else if (type=="ks") {
        ech1 <- sample(data[cat==unique(cat)[1]],replace=TRUE)
        ech2 <- sample(data[cat==unique(cat)[2]],replace=TRUE)
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
        synth$groups <- data.frame(categories=unique(cat),group=starss)
        synth$p.value <- pvals
      } else {
        synth <- list()
        synth$groups <- data.frame(categories=unique(cat),group=c("",""))
        synth$p.value <- pvals
      }
    } else if (ctrl==FALSE) {
      synth <- list()
      if (pvals <= alpha) {
        synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
        synth$p.value <- pvals
      } else {
        synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
        synth$p.value <- pvals
      }
    }
    return(synth)
  }
  normality <- function(data,cat) {
	subnormality <- function(vector) {
		if (length(vector)<=100) {
			return(shapiro.test(as.numeric(vector))$p.value)
		} else if (length(vector)<=1000) {
			return(jb.norm.test(vector)$p.value)
		} else {
			return(1)
		}
	}
	pvals <- by(data,cat,subnormality)
	return(pvals)
  }
  skew <- function(vector) {return(abs(skewness(vector)))}
  skew2 <- function(vector) {return(skewness.norm.test(vector)$p.value)}
  kurto <- function(vector) {if (is.na(abs(kurtosis(vector)))){return(10)} ; return(abs(kurtosis(vector)))}
  kurto2 <- function(vector) {return(kurtosis.norm.test(vector)$p.value)}
  if (debug==TRUE) {print("Detection of NA, Inf or unvariabilities.")}
  if (paired==TRUE) {stop("Error! The paired analysis is not developped.\n")}
  if (any((is.na(data)))|(any((is.na(cat))))){
    if (verbose==TRUE) {cat("Warning! Missing values.\n")}
    temp <- data.frame(data,cat)
    temp <- na.omit(temp)
    data <- temp[,1]
    cat <- temp[,2]
  }
  if(max(by(data,cat,length),na.rm=T)<3) {
    if (verbose==TRUE) {cat("Error! No enough values in the samples.\n")}
    return(1)
  }
  if(any(!is.finite(data))) {
	stop("Error! Infinite values in data. Analysis impossible.\n")
  }
  if (any(is.na(by(data,cat,length)))) {
	warning("Warning! Some levels do not have corresponding data.")
	cat <- factor(cat)
	if (length(unique(cat))<2) {stop("Not enough levels presenting data.")}
  }
  if(min(by(data,cat,length),na.rm=T)<3) {
    if (verbose==TRUE) {warning("Warning! No enough values for some samples. The categories concerned are ignored.")}
    which(by(data,cat,length)<3)-> ind_temp
    '%notin%' <- Negate('%in%')
    data <- data[cat%notin%names(ind_temp)]
    cat <- cat[cat%notin%names(ind_temp)]
	cat <- as.factor(cat)
	droplevels(cat)->cat
  }
  if(max(by(data,cat,var,na.rm=T),na.rm=T)==0) {
    if (verbose==TRUE) {stop("Error! No variability in samples.")}
    return(1)
  }
  data2 <- data
  cat2 <- cat
  if(min(by(data,cat,var,na.rm=T),na.rm=T)==0) {
    if (verbose==TRUE) {warning("Warning! Some samples do not vary. Non-variable categories are ignored.")}
    which(by(data,cat,var,na.rm=T)==0)-> ind_temp
    '%notin%' <- Negate('%in%')
    data <- data[cat%notin%names(ind_temp)]
    cat <- cat[cat%notin%names(ind_temp)]
	cat <- as.factor(cat)
	droplevels(cat)->cat
  }
  if (length(unique(cat))<=1) {
    if (verbose==TRUE) {stop("Error! Only one category.")}
    return(1)
  }
  if (length(unique(cat))>maxcat) {
    if (verbose==TRUE) {stop("Error! Too much categories.")}
    return(1)
  }
  if (plot==TRUE) {
	if (debug==TRUE) {print("Plot.")}
    boxplot(data~cat,col="cyan")
    vioplot(data~cat,col="#00AA0077",add=TRUE)
    stripchart(data~cat,col="#FF000088",pch=16,vertical=TRUE,add=T,method="jitter",jitter=1/(length(unique(cat))+2))
  }
  discret <- discret.test(data)
  pvals <- normality(data,cat)
  ##########################
  # Correction de Sidak
  if (debug==TRUE) {print("Sidak's correction . Bonferroni would be more secure.")}
  pval <- 1-(1-alpha)^(1/length(unique(cat)))
  ##########################
  if (code==TRUE){
	cat(paste0("alpha1 <- 1-(1-",alpha,")^(1/length(unique(cat))) # Sidak correction of alpha for test repetitions (like Shapiro, Skewness...)\n"))
	if (min(by(data,cat,length))<=100) {
		cat("by(data,cat,shapiro.test)#1) Control of the normality of small samples (<100)\n")
	}
	if (any((by(data,cat,length)<=1000)&(by(data,cat,length)>100))) {
		cat("#1A)\nby(data,cat,jb.norm.test)#1B) Control of the normality of big samples (<1000)\n")
	}
	if (max(by(data,cat,length))>1000) {
		cat("#1) Central limit theory: enough values for some samples to not have to check normality.\n")
	}
  }
  if (min(pvals) > pval) { # NORMAL
    ###################################################
    #			NORMAL
    ###################################################
    if (verbose==TRUE) {cat("1) Control of normality - Shapiro-Wilk test (<=100), Jarque-Bera (<=1000) or nothing (>1000).\n\tAnalyse of values by sample (shapiro.test() & jb.norm.test()) - \n\tThe samples follow the normal law. min(p-value):",min(pvals),"with Sidak correction of alpha to ",pval,"\n")}
    if (length(unique(cat))==2) { # 2 categories
      if (code==TRUE){cat("length(unique(cat)) #2)\n")}
      if (verbose==TRUE) {cat("2) Two categories.\n")}
      formula <- formula(data~cat)
      pvals <- var.test(formula)$p.value
      if (code==TRUE){cat("var.test(data~cat) #3)\n")}
      if (pvals>alpha) {
        ###################################################
        #			NORMAL		2 categories	homogene variance
        ###################################################
        if (verbose==TRUE) {cat("3) Fisher-Snedecor test (var.test()) - Identical sample variances. p-value:",pvals,"\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE)$p.value
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE) #4)\n")}
        if (pvals <= alpha) {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              starss <- c("","")
              starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",ifelse(pvals <=0.05,"*","")))
              synth$groups <- data.frame(categories=unique(cat),group=starss)
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,type="mean",var.equal=TRUE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)}
        } else {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Non-significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),group=c("",""))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)}
        }
      } else {
        ###################################################
        #			NORMAL		2 categories	non-homogene variance
        ###################################################
        if (verbose==TRUE) {cat("3) Fisher-Snedecor test (var.test()) - Non-identical sample variances. p-value:",pvals,"\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],,var.equal=FALSE)$p.value
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE) #4)\n")}
        if (pvals <= alpha) {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              starss <- c("","")
              starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",                                                                       ifelse(pvals <=0.05,"*","")))
              synth$groups <- data.frame(categories=unique(cat),group=starss)
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)}
        } else {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Non-significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),group=c("",""))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)}
        }
      }
    } else { 													# > 2 categories
      ###################################################
      #			NORMAL		>2 categories
      ###################################################
      if (code==TRUE){cat("length(unique(cat))#2)\n")}
      if (verbose==TRUE) {cat("2) More than two categories.\n")}
      pvals <- bartlett.test(data,cat)$p.value
      if (code==TRUE){cat("bartlett.test(data,cat) #3)\n")}
      if (pvals > alpha) {											# Identical variances
        if (verbose==TRUE) {cat("3) Bartlett test (bartlett.test()) - Identical sample variances. p-value:",pvals,"\n")}
        formula <- formula(data~cat)
        mya <- suppressWarnings(aov(data.frame(data,cat), formula=formula))
        if (code==TRUE){cat("mya <- aov(data.frame(data,cat), formula=data~cat) #4)\n")}
        pvals <- summary(mya)[[1]][["Pr(>F)"]][1]
        if (pvals<=alpha) {										# Significant AOV
          if (verbose==TRUE) {cat("4) One-way analysis of variance (aov()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (code==TRUE){cat("library(agricolae)#5a)\nSNK.test(mya,'cat',alpha=",alpha,"))#5b)\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=alpha)
            if (verbose==TRUE) {cat("5) Post-hoc Student and Newman-Keuls tests (pairwise.t.test() & SNK.test() for Newman-Keuls) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE, alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
            # if control == TRUE
			if (length(control)>0) {
				synth$Dunnett <- DunnettTest(data,cat, control = control)
				if (verbose==TRUE) {cat("6) See also post-hoc Dunnett test for the control (DunnettTest() from {DescTools}).\n")}
				if (code==TRUE){cat("library(DescTools)#6a)\nDunnettTest(data,cat,control=",control,")#6b)\n")}
				}
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {warning("Warning! pairwise.t.test() and SNK.test() don't return the same number of groups.")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              warning("Warning! Bootstrap detects weaknesses in the significance of the results.")
            }
            return(synth)
          } else {return(pvals)}
        } else {											#	 Non-significant AOV
          if (verbose==TRUE) {cat("4) One-way analysis of variance (aov()) - Non-significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=alpha)
            if (verbose==TRUE) {cat("5) Post-hoc Student & Newman_Keuls test (pairwise.t.test() & SNK.test() for Newman-Keuls) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
			if (length(control)>0) {
				synth$Dunnett <- DunnettTest(data,cat, control = control)
				if (verbose==TRUE) {cat("6) See also post-hoc Dunnett test for the control (DunnettTest() from {DescTools}).\n")}
				if (code==TRUE){cat("library(DescTools)#6a)\nDunnettTest(data,cat,control=",control,")#6b)\n")}
				}
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(pvals)}# FALSE
        }
      } else {												# Non-identical variances
        if (verbose==TRUE) {cat("3) Bartlett test (bartlett.test()) - Non-identical sample variances. p-value:",pvals,"\n")}
        pvals <- oneway.test(data~cat,var.equal=FALSE)$p.value
        if (code==TRUE){cat("oneway.test(data~cat,var.equal=FALSE) #4)\n")}
        myf <- try(fanova.hetero(data.frame(data,cat = as.factor(cat)),data~cat),silent=silent)
        if (is(myf)=="try-error") {
          #if (verbose==TRUE) {cat("Error on fanova.hetero()\n")}
          pvals2 <- alpha
        } else {pvals2 <- myf$ans[4]}
        if (pvals<=alpha) {
          if (verbose==TRUE) {cat("4) Welch’s heteroscedastic F test (oneway.test(var.equal=FALSE)) Significant differences between samples p-value:",pvals,"\n")}
          if (pvals2 > alpha) {
            if (verbose==TRUE) {cat("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:",pvals2,"\n")}
          }
          if (code==TRUE){cat("result <- pairwise.t.test(data,cat,pool.sd=FALSE)#5a)\nlibrary(KefiR)#5b)\ncatego(result)#5c)\n")}
          if (return==TRUE) {
            synth <- pairwise(data,cat,type="mean",pool.sd=FALSE,alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
            if (verbose==TRUE) {cat("5) Post-hoc Student test (pairwise.t.test(pool.sd=FALSE))\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(pvals)}# TRUE
        } else {
          if (verbose==TRUE) {cat("4) Welch’s heteroscedastic F test (oneway.test(var.equal=FALSE)) Non-significant differences between samples. p-value:",pvals,"\n")}
          if (pvals2 <= alpha) {
            if (verbose==TRUE) {cat("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:",pvals,"\n")}
          }
          if (return==TRUE) {
            synth <- pairwise(data,cat,type="mean",pool.sd=FALSE,alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
            if (verbose==TRUE) {cat("5) Post-hoc Student test (pairwise.t.test(pool.sd=FALSE))\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(pvals)}# FALSE
        }
      }
    }
    ###################################################
    #			NON-NORMAL
    ###################################################
  } else { 												#
    if (verbose==TRUE) {cat("1) Control of normality - Shapiro-Wilk test (<=100), Jarque-Bera (<=1000) or nothing (>1000).\n\tAnalyse of values by sample (shapiro.test() & jb.norm.test()) -\t\nOne or more non-normal samples. min(p-value) : ",min(pvals),"with Sidak correction of alpha to ",pval,"\n")}
	sd_cr <- by(data,cat,sd,na.rm=T) ; median_cr <- by(data,cat,median,na.rm=T)
	data_cr <- data
	for (i in names(median_cr)) {
		data_cr[cat==i] <- (data_cr[cat==i]-median_cr[which(names(median_cr)==i)])/sd_cr[which(names(sd_cr)==i)]
	}
    temp <- pairwise(data_cr,cat,type="ks",silent=silent,boot=FALSE)
    ks_result <- min(unlist(temp$p.value),na.rm=T)
    if (code==TRUE){cat("length(unique(cat))#2)\n")}
    ###################################################
    #			NON-NORMAL		2 categories
    ###################################################
    if (length(unique(cat))==2) { 							# 2 categories
      if (verbose==TRUE) {cat("2) Two categories.\n")}
      sk <- max(by(data,cat,skew))
	  sk2 <- min(by(data,cat,skew2))
      ku <- max(by(data,cat,kurto))
	  ku2 <- min(by(data,cat,kurto2))
	  #jb <- min(by(data,cat,jarquebare))
      tt <- min(by(data,cat,length))
      if (code==TRUE){cat("library(agricolae)\nby(data,cat,skewness)\nby(data,cat,skewness.norm.test)\nby(data,cat,kurtosis)\nby(data,cat,kurtosis.norm.test)\nby(data,cat,length)#3)\n")}
      ###################################################
      #			NON-NORMAL		2 categories		Acceptable for t.test()
      ###################################################
      if ((sk2>pval)&(ku2>pval)&(tt>100)&(discret>0.05)) { #(jb>alpha)
        if (verbose==TRUE) {cat("3) Jarque–Bera test & Skweness & Kurtosis limits and Length of sample (jb.norm.test() & swkeness() & skewness.norm.test() & kurtosis() & kurtosis.norm.test() & length()) - Distribution and length of samples acceptable.\n")
		  #cat("\tJarque-Bera test - normality acceptable (min(p.value)) :",jb,"\n")
          cat("\tSkweness limite (max and absolute):",sk,"\n")
		  cat("\tSkweness bootstrapped (min(p.value)):",sk2,"\n")
          cat("\tKurtosis limite (max and absolute):",ku,"\n")
		  cat("\tKurtosis bootstrapped (min(p.value)):",ku2,"\n")
          cat("\tSample length (minimal):",tt,"\n")
        }
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=FALSE)   #4)\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=FALSE)$p.value
        if (pvals <= alpha) {
          if (verbose==TRUE) {cat("4) Student Test (t.test())- Significant differences between samples.\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              starss <- c("","")
              starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",	ifelse(pvals <=0.05,"*","")))
              synth$groups <- data.frame(categories=unique(cat),group=starss)
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }

              return(synth)
            }
          } else {return(pvals)} # TRUE
        } else {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Non-significant differences between samples.\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),group=c("",""))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)} # FALSE
        }
        ###################################################
        #			NON-NORMAL		2 categories		Non acceptable for t.test()
        ###################################################
      } else {
        ###################################################
        #			NON-NORMAL		2 categories		Non acceptable for t.test()		Discret
        ###################################################
        if (verbose==TRUE) {
          if (discret <= 0.05) {
            cat("3) The number of unique values suggests the presence of discrete data : ",discret*100,"%\n")
            pvals <- mood.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])$p.value
            if (pvals <= alpha) {
              cat("3) Brown and Mood test (mood.test()) - The medianes are different and data need to be centered on the mediane for ansari.test(). p-value : ",pvals,"\n")
              by(data,cat,function(x){return(x-median(x))})->cent_med
              pvals <- ansari.test(unlist(cent_med[1]),unlist(cent_med[2]))$p.value
            } else {
              cat("3) Brown and Mood test (mood.test()) - The medianes are the same. p-value : ",pvals,"\n")
              pvals <- ansari.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])$p.value
            }
            if (pvals < alpha) {
              cat("4) Ansari-Bradley test (ansari.test()) - Data do not have the same variance. p-value: ",pvals,"\n")
            } else {
              cat("4) Ansari-Bradley test (ansari.test()) - Data have the same variance. p-value: ",pvals,"\n")
            }
            ###################################################
            #			NON-NORMAL		2 categories		Non acceptable for t.test()		Wilcox comparison
            ###################################################
          } else {
            if ((sk2<pval)|(ku2<pval)|(tt<100)) { # (jb<alpha)
              cat("3) Skweness & Kurtosis limits and Length of sample (swkeness() & skewness.norm.test() & kurtosis() & kurtosis.norm.test() & length()) - Bad distribution of data (asymmetry, spread) or insufficient length.\n")
			  #cat("\tJarque-Bera test - normality acceptable (min(p.value)) :",jb,"\n")
			  cat("\tSkweness limite (max and absolute):",sk,"\n")
			  cat("\tSkweness bootstrapped (min(p.value)):",sk2,"to compare to Sidak's corrected alpha ",pval,"\n")
			  cat("\tKurtosis limite (max and absolute):",ku,"\n")
			  cat("\tKurtosis bootstrapped (min(p.value)):",ku2,"to compare to Sidak's corrected alpha ",pval,"\n")
			  cat("\tSample length (minimal):",tt,"\n")
            }
            if (ks_result < pval) {
              cat("4) Kolmogorov-Smirnov test (ks.test()) on median-centered and reduced data -\n\tWarning! the data do not have the same distribution. p-value: ",ks_result,"\n\tby comparing with Sidak corrected alpha ",pval,"\n\tThe Mann-Whitney-Wilcoxon test will be less reliable.\n\tWarning! For wilcox.test() : Please, check graphically that the samples have the same distribution.\n")
            } else {
              cat("4) Kolmogorov-Smirnov test (ks.test()) on median-centered and reduced data -\n\tThe samples have the same distribution. p-value: ",ks_result,"\n\tby comparing with Sidak corrected alpha ",pval,"\n\tThe Mann-Whitney-Wilcoxon test will be reliable.\n")
            }
          }
        }
        if (code==TRUE){cat("wilcox.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])  #2)\n")}
        pvals <- suppressWarnings(wilcox.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]]))$p.value
        if (pvals <= alpha) {
          if (verbose==TRUE) {cat("5) Wilcoxon-Mann-Whitney test (wilcox.test()) - Significant differences between samples. p-value: ",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="median",conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              starss <- c("","")
              starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",
                                                                        ifelse(pvals <=0.05,"*","")))
              synth$groups <- data.frame(categories=unique(cat),group=starss)
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)} # TRUE
        } else {
          if (verbose==TRUE) {cat("5) Wilcoxon-Mann-Whitney test (wilcox.test()) - Non-significant differences between samples. p-value : ",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),group=c("",""))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,alpha=alpha)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(pvals)} # FALSE
        }
      }
      ###################################################
      #			NON-NORMAL		>2 categories
      ###################################################
    } else { 											# > 2 categories
      if (verbose==TRUE) {cat("2) More than two categories.\n")}
      sk <- max(by(data,cat,skew))
	  sk2 <- min(by(data,cat,skew2))
      ku <- max(by(data,cat,kurto))
	  ku2 <- min(by(data,cat,kurto2))
	  #jb <- min(by(data,cat,jarquebare))
      if (code==TRUE){
        #cat("#3a)\nby(data,cat,jb.norm.test)#3b)\nlibrary(agricolae)#3c)\nby(data,cat,skewness)#3d)\nby(data,cat,skewness.norm.test)#3e)\nby(data,cat,kurtosis)#3f)\nby(data,cat,kurtosis.norm.test)#3g)\nby(data,cat,length)#3h)\n")
		cat("#library(agricolae)#3a)\nby(data,cat,skewness)#3b)\nby(data,cat,skewness.norm.test)#3c)\nby(data,cat,kurtosis)#3d)\nby(data,cat,kurtosis.norm.test)#3e)\nby(data,cat,length)#3f)\n")
        cat("library(lawstat)#4a)\nlevene.test(data,cat)#4b)\n")
        cat("library(onewaytests)#5a)\nbf.test(data~cat,data=data.frame(data,'cat'=factor(cat)))#5b)\n")
      }
      pvals2 <- suppressWarnings(bf.test(data~cat,data=data.frame(data,"cat"=factor(cat)),verbose=FALSE))$p.value
      if (verbose==TRUE) {
        pvals <- suppressWarnings(levene.test(data,cat))$p.value
        if (((pvals <= alpha)&(pvals2 <= alpha))|((pvals > alpha)&(pvals2 > alpha))) {
          if (verbose==TRUE) {
            cat("3) Consistency in testing of Levene and Brown-Forsyth results.\n")
            cat("\tLevene p-value : ",pvals," - Brown-Forsyth p-value : ",pvals2,"\n")
          }
        } else {
          if (verbose==TRUE) {
            cat("3) Inconsistent testing of Levene and Brown-Forsyth.\n")
            cat("\tLevene p-value : ",pvals," - Brown-Forsyth p-value : ",pvals2,"\n")
			#cat("3') Only Jarque-Bera is taking in account.\n")
			#cat("\tJarque-Bera test - normality acceptable (min(p.value)) :",jb,"\n")
          }
        }
#		if (jb <= alpha) {
#			if (verbose==TRUE) {
#				cat("3') Only Jarque-Bera is taking in account.\n")
#				cat("\tJarque-Bera test - normality non-acceptable (min(p.value)) :",jb,"\n")
#			}
#		} else {
#			if (verbose==TRUE) {
#				cat("3') Only Jarque-Bera is taking in account.\n")
#				cat("\tJarque-Bera test - normality acceptable (min(p.value)) :",jb,"\n")
#			}
#		}
      }
      #print(pvals)
#      if (jb>pval) {					# Acceptable non-normality
#        if (verbose==TRUE) {cat("4) Jarque-Bare & Skweness or Kurtosis limits & Brown-Forsyth test (swkeness() & kurtosis() & bf.test()) - The distribution of values and sample variances are acceptable.\n")
#			  cat("\tSkweness limite (max and absolute):",sk,"\n")
#			  cat("\tSkweness bootstrapped (min(p.value)):",sk2,"\n")
#			  cat("\tKurtosis limite (max and absolute):",ku,"\n")
#			  cat("\tKurtosis bootstrapped (min(p.value)):",ku2,"\n")
#        }
#        if (code==TRUE){cat("mya <- aov(data.frame(data,cat), formula=data~cat) #6)\n")	}
#        formula <- formula(data~cat)
#        mya <- aov(data.frame(data,cat), formula=formula)
#        pvals <- summary(mya)[[1]][["Pr(>F)"]][1]
#        if (pvals<=alpha) {
#          if (verbose==TRUE) {cat("5) Oneway analysis of variance (aov()) - Significant differences between samples. p-value:",pvals,"\n")}
#          if (return==TRUE) {
#            mynk <- SNK.test(mya,"cat",alpha=alpha)
#            if (verbose==TRUE) {cat("6) Post-hoc tests Student and Newman-Keuls (pairwise.t.test() & SNK.test()) \n")}
#            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
#            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
#            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
#            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
#            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
#            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
#            which(unique(cat)==control)-> ind_control
#            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
#            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
#              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
#            }
#            return(synth)
#          } else {return(pvals)} # TRUE
#        } else {
#          if (verbose==TRUE) {cat("5) Oneway test analysis of variance (aov()) - Non-significant differences between samples. p-value:",pvals,"\n")}
#          if (code==TRUE){cat("library(agricolae)#7a)\nprint(SNK.test(mya,'cat',alpha=",alpha,"))#7b)\n")}
#          if (return==TRUE) {
#            mynk <- SNK.test(mya,"cat",alpha=alpha)
#            if (verbose==TRUE) {cat("6) Post-hoc tests Student and Newman-Keuls (pairwise.t.test() & SNK.test()) \n")}
#            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
#            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
#            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
#            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
#            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
#            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
#            which(unique(cat)==control)-> ind_control
#            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1) ) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
#            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
#              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
#            }
#            return(synth)
#          } else {return(pvals)} # FALSE
#        }
#     } else {
        if (code==TRUE){cat("kruskal.test(data,cat) #6)\n")}
        pvals3 <- kruskal.test(data,cat)$p.value
        # Si verbose : blabla
        if (verbose==TRUE) {
          if ((sk2<pval)|(ku2<pval)) {
            cat("4) Skweness & Kurtosis limits (swkeness() & kurtosis()) Bad distribution of data (asymmetry, spread).\n")
			  cat("\tSkweness limite (max and absolute):",sk,"\n")
			  cat("\tSkweness bootstrapped (min(p.value)):",sk2,"\n")
			  cat("\tKurtosis limite (max and absolute):",ku,"\n")
			  cat("\tKurtosis bootstrapped (min(p.value)):",ku2,"\n")
          }
        }
        pvals <- fligner.test(data,cat)$p.value
        if (is.na(pvals)) {
          if (verbose==TRUE) {cat("5) Fligner-Killeen test (fligner.test()) - Error, return NA.\n")}
          return(pvals)# FALSE
        }
        if (verbose==TRUE) {
          if (pvals<=alpha) {
            cat("5) Fligner-Killeen test (fligner.test())Significant differences of variance between samples. p-value:",pvals,"\n")
          } else {
            cat("5) Fligner-Killeen test (fligner.test())Non-significant differences of variance between samples. p-value",pvals,"\n")
          }
          #temp <- pairwise(data,cat,type="ks",silent=silent,boot=boot)$p.value
          #ks_result <- min(unlist(temp),na.rm=TRUE)
          if (ks_result < pval) {
            cat("6) Kolmogorov-Smirnov test (ks.test()) on median-centered and reduced data -\n\tWarning! the samples do not have the same distribution. min(p-value) : ",ks_result,"\n\tby comparing with Sidak corrected alpha ",pval,"\n\tThe Kruskal-Wallis test and Mann-Whitney-Wilcoxon test will be less reliable.\n\tPlease, check graphically the samples distributions.\n")
          } else {
            cat("6) Kolmogorov-Smirnov test (ks.test()) on median-centered and reduced data -\n\tThe samples have the same distribution. min(p-value) : ",ks_result,"\n\tby comparing with Sidak corrected alpha ",pval,"\n\tGood accuracy expected on the tests of Kruskal-Wallis and Mann-Whitney-Wilcoxon\n")
          }
          if (pvals3 <= alpha) {
            cat("7) Kruskal-Wallis test (kruskal.test()) - At least one sample appears to show a difference. p-value:",pvals3,"\n")
          } else {
            cat("7) Kruskal-Wallis test (kruskal.test()) - No different sample a priori. p-value:",pvals3,"\n")
          }
        }
        if ((return==TRUE) | (verbose==TRUE)) {
          if (pvals<=alpha) {
            pvals <- med1way(data~cat)$p.value
            if (is.na(pvals)) {
              if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.\n")}
            } else {
              if (pvals <= alpha) {
                if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Significant differences between the median of samples. p-value:",pvals,"\n")
                  if (pvals3 > alpha) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on medians give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(pvals)}
              } else if (pvals > alpha) {
                if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Non-significant differences between the median of samples. p-value:",pvals,"\n")
                  if (pvals3 <= alpha) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on medians give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(pvals)}
              }
            }
            if (code==TRUE){cat("pairwise.wilcox.test(data,cat,p.adjust.method='BH') #7)\n")}
            if (return==TRUE) {
              if (verbose==TRUE) {cat("9) Wilcoxon-Mann-Whitney test (pairwise.wilcox.test()) \n")}
              synth <- pairwise(data2,cat2,type="median",alpha=alpha,control=control,boot=boot,conf=conf,iter=iter)
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              if (pvals3 <= alpha) {
                return(pvals3)
              } else {return(pvals3)}
            }
          } else {
            pvals <- t1way(data~cat)$p.value
            if (is.na(pvals)) {
              if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.\n")}
            } else {
              if (pvals <= alpha) {
                if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Significant differences between the trimmed samples. p-value:",pvals,"\n")
                  if (pvals3 > alpha) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on trimmed means give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(pvals)}
              } else if (pvals > alpha) {
                if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Non-significant differences between the trimmed samples. p-value:",pvals,"\n")
                  if (pvals3 <= alpha) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on trimmed means give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(pvals)}
              }
            }
            if (code==TRUE){cat("library(WRS2) #7a)\nlincon(data~cat) #7b)\n")}
            if (return==TRUE) {
              if (verbose==TRUE) {cat("9) Correspondance post-hoc on trimmed means (lincon())\n")}
              synth <- pairwise(data,cat,type="lincon",alpha=alpha,control=control)
              return(synth)
            } else {
              if (pvals3 <= alpha) {
                return(pvals3)
              } else {return(pvals3)}
            }
          }
        } else if (return==FALSE) {
          if (pvals3 <= alpha) {return(pvals3)
          }else {return(pvals3)}
        }
#      }
    }
  }
}

#' Determining groups based on the results of a pairwise function
#'
#' @param result output results of the pairwise(), pairwise.t.test() or pairwise.wilcox.test() function
#' @param control name of the category that will be used as a control to establish differences with '*', '**' and '***'.
#' @param pval threshold of p-value to establish groups of averages or medians of type "a", "ab", "b" ...
#'
#' @return catego() can be used to establish groups with significantly different mean/median values relative to each other or to a reference control. It exploits the results of the functions pairwise.t.test(), pairwise.wilcox.test() and pairwise(type = "median" or "mean").
#' @export
#'
#' @examples
#' data(iris)
#' output <- pairwise.t.test(iris[,1],iris$Species)
#' catego(output)
#' catego(output,control="setosa")
catego <- function(result,control=c(),pval=0.05) {
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
          if (pvals > pval){
            #print("non-signif")
            cats <- c(cats,possibility[pos])
          } else {
            #print("signif")
            cats <- c(cats,possibility[pos+1])
            change <- 1
          }
        } else {
          if (pvals <= pval){
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
              cats[i] <- paste0(cats[i],cats[j+1])
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

#' Automation function for pairwise calculations.
#'
#' @param x numerical vector
#' @param g category vector
#' @param type 'mean' for pairwise.t.test(p.adjust.method="holm"), 'median' for pairwise.wilcox.test(p.adjust.method="BH"), 'ks' for ks.test(), 'lincon' for lincon() of {WSR2}
#' @param pval threshold value of p-value to establish the groups.
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
#' pairwise(iris[,1],iris[,5],type="median",pval=0.01,boot=TRUE)#wilcox
#' pairwise(iris[,1],iris[,5],type="ks")
#' pairwise(iris[,1],iris[,5],type="lincon")
pairwise <- function(x,g,type="mean",pval=0.05,control=c(),pool.sd=FALSE,silent=TRUE,boot=FALSE,iter=500,conf=0.95) {
  onoff <- function(x,silent=FALSE,oldw="") {
    if (silent==TRUE) {
      if (x == "on") {
        oldw <- getOption("warn")
        options(warn = -1)
      } else if (x == "off") {
        options(warn = oldw)
      }
    }
    return(oldw)
  }
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
      groups <- catego(result,pval=pval)
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
        groups$bootstrap <- catego(output,pval=pval)
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
      groups <- catego(result,pval=pval)
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
        groups$bootstrap <- catego(output,pval=pval)
      }
    }
    match(init_order,levels(g))-> indices
    groups$groups <- groups$groups[indices,]
    if (boot==TRUE) {groups$bootstrap$groups <- groups$bootstrap$groups[indices,]}
    return(groups)
  } else if (type == "ks" ){
    unique_g <- unique(g)
    for (i in unique_g) {
      x[g==i] <- (x[g==i]-mean(x[g==i]))/sd(x[g==i])
    }
    mymat <- matrix(rep(NA,(length(unique_g)-1)^2),nc=(length(unique_g)-1),nr=(length(unique_g)-1))
    rownames(mymat) <- unique_g[2:length(unique_g)] ; colnames(mymat) <- unique_g[1:(length(unique_g)-1)]
    ks_func <- function(x,g,mymat,unique_g) {
      for (i in 1:(length(unique_g)-1)) {
        for (j in (i+1):length(unique_g)) {
          pv <- ks.test(x[g==unique_g[i]],x[g==unique_g[j]])$p.value
          #print(pv)
          mymat[(j-1),i] <- pv
        }
      }
      return(mymat)
    }
    onoff("on",silent)->oldw
    if (boot == FALSE) {mymat <- ks_func(x,g,mymat,unique_g) ; output <- list() ; output$p.value <- mymat
    } else if (boot==TRUE) {
      mydim <- dim(mymat)
      myarray <- array(rep(NA,iter*mydim[1]*mydim[2]),c(iter,mydim))
      for (k in 1: iter) {
        x_temp <- x
        for (j in levels(g)) {
          x_temp[g==j] <- sample(x_temp[g==j],replace=TRUE)
        }
        temp <- ks_func(x_temp,g,mymat,unique_g)
        myarray[k,,] <- temp
      }
      output <- list()
      apply(myarray,c(2,3),quantile,probs=conf,na.rm=T)->output$p.value
      colnames(output$p.value) <- colnames(mymat)
      rownames(output$p.value) <- rownames(mymat)
    }
    onoff("off",silent,oldw)
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
      groups <- catego(result,pval=pval)
    }
    match(init_order,levels(g))-> indices
    groups$groups <- groups$groups[indices,]
    return(groups)
  } else {stop("Error! arg type is 'mean','median', 'ks' or 'lincon'.\n")}
}

#' Automatic average/median/variance comparison function.
#'
#' @param data numerical vector
#' @param cat category vector
#' @param pval p-value threshold value for all the tests (except ks.test() which is too sensitive).
#' @param verbose to display the full reasoning of the analysis.
#' @param return allows to return the results of pairwise analysis (p-values and groups).
#' @param paired (under development) to allow the analysis of matched data.
#' @param control name of the category that will eventually be used as a control.
#' @param pval_ks minimal p-value of automatic analyse of distribution by ks.test(). The minimum p-value of the ks.test() also fluctuates when boot=TRUE.
#' @param maxcat maximum number of categories allowed. When this number is high, some tests may return an error message.
#' @param plot to display the distribution of the data.
#' @param silent for displaying or not warnings.
#' @param boot to activate the boostrap on 'mean' and 'median'.
#' @param iter number f iterations (boot==TRUE).
#' @param conf confidence level of bootstrap.
#' @param code allows to display the simplified R source code to be able to do the same R study step by step.
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
#' @import methods
#' @export
#'
#' @examples
#' data(iris)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE)
#' m.test(iris[,4],iris[,5],verbose=TRUE, plot=FALSE, return=FALSE, boot=FALSE)
#' m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE,control="virginica")
#' m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
#' m.test(iris[,4],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
m.test <- function (data, cat, pval=0.05, verbose=TRUE, return=TRUE, paired=FALSE,control=c(),
                    pval_ks = 0.01, maxcat=50, plot=TRUE,silent=TRUE,boot=TRUE,iter=500,conf=0.95,
                    code=FALSE){
  if (code==TRUE) {verbose <- FALSE}
  discret.test <- function(vector) {
    return(length(unique(vector))/length(vector))
  }
  onoff <- function(x,silent=FALSE,oldw="") {
    if (silent==TRUE) {
      if (x == "on") {
        oldw <- getOption("warn")
        options(warn = -1)
      } else if (x == "off") {
        options(warn = oldw)
      }
    }
    return(oldw)
  }
  boots <- function(data,cat,ctrl=FALSE,type="mean",var.equal=FALSE,conf=0.95,iter=500,pval=pval) {
    pvals <- c()
    for (i in 1:iter) {
      if (type=="mean") {
        temp <- t.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE),var.equal=var.equal)$p.value
      } else if (type=="median") {
        temp <- wilcox.test(sample(data[cat==unique(cat)[1]],replace=TRUE),sample(data[cat==unique(cat)[2]],replace=TRUE))$p.value
      } else if (type=="ks") {
        ech1 <- sample(data[cat==unique(cat)[1]],replace=TRUE)
        ech2 <- sample(data[cat==unique(cat)[2]],replace=TRUE)
        ech1 <- (ech1 - mean(ech1))/sd(ech1)
        ech2 <- (ech2 - mean(ech2))/sd(ech2)
        temp <- ks.test(ech1,ech2)$p.value
      }
      pvals <- c(pvals,temp)
    }
    pvals <- quantile(pvals,probs=conf,na.rm=T)
    if (ctrl == TRUE) {
      if (pvals <= pval) {
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
      if (pvals <= pval) {
        synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
        synth$p.value <- pvals
      } else {
        synth$groups <- data.frame(categories=unique(cat),groups=c("a","a"))
        synth$p.value <- pvals
      }
    }
    return(synth)
  }
  shapi <- function(vector) {return(shapiro.test(as.numeric(vector))$p.value)}
  skew <- function(vector) {return(abs(skewness(vector)))}
  kurto <- function(vector) {if (is.na(abs(kurtosis(vector)))){return(10)} ; return(abs(kurtosis(vector)))}
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
    return(FALSE)
  }
  if(min(by(data,cat,length),na.rm=T)<3) {
    if (verbose==TRUE) {cat("Warning! No enough values for some samples. The categories concerned are eliminated..\n")}
    which(by(data,cat,length)<3)-> ind_temp
    '%notin%' <- Negate('%in%')
    data <- data[cat%notin%names(ind_temp)]
    cat <- cat[cat%notin%names(ind_temp)]
  }
  if(max(by(data,cat,var,na.rm=T),na.rm=T)==0) {
    if (verbose==TRUE) {cat("Error! No variability in samples.\n")}
    return(FALSE)
  }
  if(min(by(data,cat,var,na.rm=T),na.rm=T)==0) {
    if (verbose==TRUE) {cat("Warning! Some samples do not vary. Non-variable categories are eliminated\n")}
    which(by(data,cat,var,na.rm=T)==0)-> ind_temp
    '%notin%' <- Negate('%in%')
    data <- data[cat%notin%names(ind_temp)]
    cat <- cat[cat%notin%names(ind_temp)]
  }
  if (length(unique(cat))<=1) {
    if (verbose==TRUE) {cat("Error! Only one category.\n")}
    return(FALSE)
  }
  if (length(unique(cat))>maxcat) {
    if (verbose==TRUE) {cat("Error! Too much categories.\n")}
    return(FALSE)
  }
  if (plot==TRUE) {
    boxplot(data~cat,col="cyan")
    vioplot(data~cat,col="#00AA0077",add=TRUE)
    stripchart(data~cat,col="#FF000088",pch=16,vertical=TRUE,add=T,method="jitter",jitter=1/(length(unique(cat))+2))
  }
  discret <- discret.test(data)
  pvals <- by(data,cat,shapi)
  if (code==TRUE){cat("by(data,cat,shapiro.test)#1)\n")}
  if (min(pvals) > pval) { # NORMAL
    ###################################################
    #			NORMAL
    ###################################################
    if (verbose==TRUE) {cat("1) Shapiro-Wilk test (shapiro.test()) - The samples follow the normal law. min(p-value):",min(pvals),"\n")}
    if (length(unique(cat))==2) { # 2 categories
      if (code==TRUE){cat("length(unique(cat)) #2)\n")}
      if (verbose==TRUE) {cat("2) Two categories.\n")}
      formula <- formula(data~cat)
      pvals <- var.test(formula)$p.value
      if (code==TRUE){cat("var.test(data~cat) #3)\n")}
      if (pvals>pval) {
        ###################################################
        #			NORMAL		2 categories	homogene variance
        ###################################################
        if (verbose==TRUE) {cat("3) Fisher-Snedecor test (var.test()) - Identical sample variances. p-value:",pvals,"\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE)$p.value
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE) #4)\n")}
        if (pvals <= pval) {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,pval=pval)
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
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,type="mean",var.equal=TRUE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(TRUE)}
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
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,pval=pval)
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
                                         type="mean",var.equal=TRUE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(FALSE)}
        }
      } else {
        ###################################################
        #			NORMAL		2 categories	non-homogene variance
        ###################################################
        if (verbose==TRUE) {cat("3) Fisher-Snedecor test (var.test()) - Non-identical sample variances. p-value:",pvals,"\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],,var.equal=FALSE)$p.value
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=TRUE) #4)\n")}
        if (pvals <= pval) {
          if (verbose==TRUE) {cat("4) Student test (t.test()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
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
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(TRUE)}
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
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
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
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(FALSE)}
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
      if (pvals > pval) {											# Identical variances
        if (verbose==TRUE) {cat("3) Bartlett test (bartlett.test()) - Identical sample variances. p-value:",pvals,"\n")}
        formula <- formula(data~cat)
        onoff("on",silent)->oldw
        mya <- aov(data.frame(data,cat), formula=formula)
        if (code==TRUE){cat("mya <- aov(data.frame(data,cat), formula=data~cat) #4)\n")}
        onoff("off",silent,oldw)
        pvals <- summary(mya)[[1]][["Pr(>F)"]][1]
        if (pvals<=pval) {										# Significant AOV
          if (verbose==TRUE) {cat("4) One-way analysis of variance (aov()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (code==TRUE){cat("library(agricolae)#5a)\nprint(SNK.test(mya,'cat',alpha=",pval,"))#5b)\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=pval)
            if (verbose==TRUE) {cat("5) Post-hoc Student and Newman-Keuls tests (pairwise.t.test() & SNK.test() for Newman-Keuls) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE, pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
            #print(ind_control)
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(TRUE)}
        } else {											#	 Non-significant AOV
          if (verbose==TRUE) {cat("4) One-way analysis of variance (aov()) - Non-significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=pval)
            if (verbose==TRUE) {cat("5) Post-hoc Student & Newman_Keuls test (pairwise.t.test() & SNK.test() for Newman-Keuls) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(FALSE)}
        }
      } else {												# Non-identical variances
        if (verbose==TRUE) {cat("3) Bartlett test (bartlett.test()) - Non-identical sample variances. p-value:",pvals,"\n")}
        pvals <- oneway.test(data~cat,var.equal=FALSE)$p.value
        if (code==TRUE){cat("oneway.test(data~cat,var.equal=FALSE) #4)\n")}
        myf <- try(fanova.hetero(data.frame(data,cat = as.factor(cat)),data~cat),silent=silent)
        if (is(myf)=="try-error") {
          #if (verbose==TRUE) {cat("Error on fanova.hetero()\n")}
          pvals2 <- pval
        } else {pvals2 <- myf$ans[4]}
        if (pvals<=pval) {
          if (verbose==TRUE) {cat("4) Welch’s heteroscedastic F test (oneway.test(var.equal=FALSE)) Significant differences between samples p-value:",pvals,"\n")}
          if (pvals2 > pval) {
            if (verbose==TRUE) {cat("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:",pvals2,"\n")}
          }
          if (code==TRUE){cat("result <- pairwise.t.test(data,cat,pool.sd=FALSE)#5a)\nlibrary(KefiR)#5b)\ncatego(result)#5c)\n")}
          if (return==TRUE) {
            synth <- pairwise(data,cat,type="mean",pool.sd=FALSE,pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            if (verbose==TRUE) {cat("5) Post-hoc Student test (pairwise.t.test(pool.sd=FALSE))\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(TRUE)}
        } else {
          if (verbose==TRUE) {cat("4) Welch’s heteroscedastic F test (oneway.test(var.equal=FALSE)) Non-significant differences between samples. p-value:",pvals,"\n")}
          if (pvals2 <= pval) {
            if (verbose==TRUE) {cat("Warning! fanova.hetero() does not give the same result as oneway.test. p-value:",pvals,"\n")}
          }
          if (return==TRUE) {
            synth <- pairwise(data,cat,type="mean",pool.sd=FALSE,pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            if (verbose==TRUE) {cat("5) Post-hoc Student test (pairwise.t.test(pool.sd=FALSE))\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(FALSE)}
        }
      }
    }
    ###################################################
    #			NON-NORMAL
    ###################################################
  } else { 												#
    if (verbose==TRUE) {cat("1) Shapiro-Wilk test (shapiro.test()) -  One or more non-normal samples. min(p-value) : ",min(pvals),"\n")}
    temp <- pairwise(data,cat,type="ks",silent=silent,boot=FALSE)
    ks_result <- min(unlist(temp$p.value),na.rm=T)
    if (code==TRUE){cat("length(unique(cat))#2)\n")}
    ###################################################
    #			NON-NORMAL		2 categories
    ###################################################
    if (length(unique(cat))==2) { 							# 2 categories
      if (verbose==TRUE) {cat("2) Two categories.\n")}
      sk <- max(by(data,cat,skew))
      ku <- max(by(data,cat,kurto))
      tt <- min(by(data,cat,length))
      if (code==TRUE){cat("library(agricolae)\nby(data,cat,skewness)\nby(data,cat,kurtosis)\nby(data,cat,length)#3)\n")}
      ###################################################
      #			NON-NORMAL		2 categories		Acceptable for t.test()
      ###################################################
      if ((sk<=2)&(ku<=2)&(tt>100)&(discret>0.05)) {
        if (verbose==TRUE) {cat("3) Skweness & Kurtosis limits and Length of sample (swkeness() & kurtosis() & length()) - Distribution and length of samples acceptable.\n")
          cat("\tSkweness limite (max and absolute):",sk,"\n")
          cat("\tKurtosis limite (max and absolute):",ku,"\n")
          cat("\tSample length (minimal):",tt,"\n")
        }
        if (code==TRUE){cat("t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=FALSE)   #4)\n")}
        pvals <- t.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]],var.equal=FALSE)$p.value
        if (pvals <= pval) {
          if (verbose==TRUE) {cat("4) Student Test (t.test())- Significant differences between samples.\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              synth <- list()
              starss <- c("","")
              starss[-ind_control] <- ifelse(pvals <=0.001,"***",ifelse(pvals <=0.01,"**",																			ifelse(pvals <=0.05,"*","")))
              synth$groups <- data.frame(categories=unique(cat),group=starss)
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=TRUE,
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }

              return(synth)
            }
          } else {return(TRUE)}
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
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
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
                                         type="mean",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(FALSE)}
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
            if (pvals <= pval) {
              cat("3) Brown and Mood test (mood.test()) - The medianes are different and data need to be centered on the mediane for ansari.test(). p-value : ",pvals,"\n")
              by(data,cat,function(x){return(x-median(x))})->cent_med
              pvals <- ansari.test(unlist(cent_med[1]),unlist(cent_med[2]))$p.value
            } else {
              cat("3) Brown and Mood test (mood.test()) - The medianes are the same. p-value : ",pvals,"\n")
              pvals <- ansari.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])$p.value
            }
            if (pvals < pval) {
              cat("4) Ansari-Bradley test (ansari.test()) - Data do not have the same variance. p-value: ",pvals,"\n")
            } else {
              cat("4) Ansari-Bradley test (ansari.test()) - Data have the same variance. p-value: ",pvals,"\n")
            }
            ###################################################
            #			NON-NORMAL		2 categories		Non acceptable for t.test()		Wilcox comparison
            ###################################################
          } else {
            if ((sk<=2)|(ku<=2)|(tt>100)) {
              cat("3) Skweness & Kurtosis limits and length of sample (swkeness() & kurtosis() & length()) Bad distribution of data (asymmetry, spread) or insufficient length.\n")
              cat("\tSkweness limite (max and absolute):",sk,"\n")
              cat("\tKurtosis limite (max and absolute):",ku,"\n")
              cat("\tSample length (minimal):",tt,"\n")
            }
            #onoff("on",silent)->oldw
            #if (boot==FALSE) {pvals <- ks.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])$p.value
            #} else if (boot==TRUE) {
            #		pvals <- boots(data,cat,ctrl=FALSE,type="ks",conf=conf,iter=iter,pval=pval)$p.value
            #}
            #onoff("off",silent,oldw)
            if (ks_result < pval_ks ) {
              cat("4) Kolmogorov-Smirnov test (ks.test()) -  Warning! the data do not have the same distribution. p-value: ",ks_result,"\n\tThe Mann-Whitney-Wilcoxon test will be less reliable.\n\tWarning! For wilcox.test() : Please, check graphically that the samples have the same distribution.\n")
            } else {
              cat("4) Kolmogorov-Smirnov test (ks.test()) - The samples have the same distribution. p-value: ",ks_result,"\n\tThe Mann-Whitney-Wilcoxon test will be reliable.\n")
            }
          }
        }
        onoff("on",silent)->oldw
        if (code==TRUE){cat("wilcox.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])  #2)\n")}
        pvals <- wilcox.test(data[cat==unique(cat)[1]],data[cat==unique(cat)[2]])$p.value
        onoff("off",silent,oldw)
        if (pvals <= pval) {
          if (verbose==TRUE) {cat("5) Wilcoxon-Mann-Whitney test (wilcox.test()) - Significant differences between samples. p-value: ",pvals,"\n")}
          if (return==TRUE) {
            which(unique(cat)==control)-> ind_control
            if (length(ind_control)!=1) {
              synth <- list()
              synth$groups <- data.frame(categories=unique(cat),groups=c("a","b"))
              synth$p.value <- pvals
              if (boot==TRUE) {
                synth$bootstrap <- boots(data,cat,ctrl=FALSE,
                                         type="median",conf=conf,iter=iter,pval=pval)
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
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(TRUE)}
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
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
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
                                         type="median",var.equal=FALSE,conf=conf,iter=iter,pval=pval)
              }
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            }
          } else {return(FALSE)}
        }
      }
      ###################################################
      #			NON-NORMAL		>2 categories
      ###################################################
    } else { 											# > 2 categories
      if (verbose==TRUE) {cat("2) More than two categories.\n")}
      sk <- max(by(data,cat,skew))
      #print(sk)
      ku <- max(by(data,cat,kurto))
      #print(ku)
      if (code==TRUE){
        cat("library(agricolae)#3a)\nby(data,cat,skewness)#3b)\nby(data,cat,kurtosis)#3c)\nby(data,cat,length)#3d)\n")
        cat("library(lawstat)#4a)\nlevene.test(data,cat)#4b)\n")
        cat("library(onewaytests)#5a)\nbf.test(data~cat,data=data.frame(data,'cat'=factor(cat)))#5b)\n")
      }
      onoff("on",silent)->oldw
      pvals2 <- bf.test(data~cat,data=data.frame(data,"cat"=factor(cat)),verbose=FALSE)$p.value
      onoff("off",silent,oldw)
      if (verbose==TRUE) {
        onoff("on",silent)->oldw
        pvals <- levene.test(data,cat)$p.value
        onoff("off",silent,oldw)
        if (((pvals <= pval)&(pvals2 <= pval))|((pvals > pval)&(pvals2 > pval))) {
          if (verbose==TRUE) {
            cat("3) Consistency in testing of Levene and Brown-Forsyth results.\n")
            cat("\tLevene p-value : ",pvals," - Brown-Forsyth p-value : ",pvals2,"\n")
          }
        } else {
          if (verbose==TRUE) {
            cat("3) Inconsistent testing of Levene and Brown-Forsyth .\n\tOnly Brown-Forsyth is taken into account..\n")
            cat("\tLevene p-value : ",pvals," - Brown-Forsyth p-value : ",pvals2,"\n")
          }
        }
      }
      #print(pvals)
      if ((sk<=2)&(ku<=2)&(pvals2>pval)) {					# Acceptable non-normality
        if (verbose==TRUE) {cat("4) Skweness & Kurtosis limits and Brown-Forsyth test (swkeness() & kurtosis() & bf.test()) - The distribution of values and sample variances are acceptable.\n")
          cat("\tSkweness limite (max and absolute):",sk,"\n")
          cat("\tKurtosis limite (max and absolute):",ku,"\n")
        }
        if (code==TRUE){cat("mya <- aov(data.frame(data,cat), formula=data~cat) #6)\n")	}
        formula <- formula(data~cat)
        mya <- aov(data.frame(data,cat), formula=formula)
        pvals <- summary(mya)[[1]][["Pr(>F)"]][1]
        if (pvals<=pval) {
          if (verbose==TRUE) {cat("5) Oneway analysis of variance (aov()) - Significant differences between samples. p-value:",pvals,"\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=pval)
            if (verbose==TRUE) {cat("6) Post-hoc tests Student and Newman-Keuls (pairwise.t.test() & SNK.test()) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1)) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(TRUE)}
        } else {
          if (verbose==TRUE) {cat("5) Oneway test analysis of variance (aov()) - Non-significant differences between samples. p-value:",pvals,"\n")}
          if (code==TRUE){cat("library(agricolae)#7a)\nprint(SNK.test(mya,'cat',alpha=",pval,"))#7b)\n")}
          if (return==TRUE) {
            mynk <- SNK.test(mya,"cat",alpha=pval)
            if (verbose==TRUE) {cat("6) Post-hoc tests Student and Newman-Keuls (pairwise.t.test() & SNK.test()) \n")}
            synth <- pairwise(data,cat,type="mean",pool.sd=TRUE,pval=pval,control=control,boot=boot,conf=conf,iter=iter)
            match(synth$groups[,1],rownames(mynk$groups)) -> ind_temp
            synth$groups <- data.frame(synth$groups , "SNK"=mynk$groups$groups[ind_temp])
            cat1 <- unique(unlist(strsplit(synth$groups$groups,"")))
            cat2 <- unique(unlist(strsplit(mynk$groups$groups,"")))
            rownames(synth$groups) <- rep(c(),nrow(synth$groups))
            which(unique(cat)==control)-> ind_control
            if ((verbose==TRUE) & (length(cat1)!=length(cat2)) & (length(ind_control)!=1) ) {cat("\tWarning! pairwise.t.test() and SNK.test() don't return the same number of groups.\n")}
            if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
              cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
            }
            return(synth)
          } else {return(FALSE)}
        }
      } else {
        if (code==TRUE){cat("kruskal.test(data,cat) #6)\n")}
        pvals3 <- kruskal.test(data,cat)$p.value
        # Si verbose : blabla
        if (verbose==TRUE) {
          if ((sk>2)|(ku>2)) {
            cat("4) Skweness & Kurtosis limits (swkeness() & kurtosis()) Bad distribution of data (asymmetry, spread).\n")
            cat("\tSkweness limite (max and absolute):",sk,"\n")
            cat("\tKurtosis limite (max and absolute):",ku,"\n")
          }
          #onoff("on",silent)->oldw
          #pvals <- levene.test(data,cat)$p.value
          #onoff("off",silent,oldw)
          #if (pvals2 <= pval) {
          #	cat("3) Brown-Forsyth test (bf.test()) Non-identical sample variances.\n")
          #}
        }
        pvals <- fligner.test(data,cat)$p.value
        if (is.na(pvals)) {
          if (verbose==TRUE) {cat("5) Fligner-Killeen test (fligner.test()) - Error, return NA.\n")}
          return(FALSE)
        }
        if (verbose==TRUE) {
          if (pvals<=pval) {
            cat("5) Fligner-Killeen test (fligner.test())Significant differences of variance between samples. p-value:",pvals,"\n")
          } else {
            cat("5) Fligner-Killeen test (fligner.test())Non-significant differences of variance between samples. p-value",pvals,"\n")
          }
          #temp <- pairwise(data,cat,type="ks",silent=silent,boot=boot)$p.value
          #ks_result <- min(unlist(temp),na.rm=TRUE)
          if (ks_result < pval_ks ) {
            cat("6) Kolmogorov-Smirnov test (ks.test()) -  Warning! the samples do not have the same distribution. min(p-value) : ",ks_result,"\n\tThe Kruskal-Wallis test and Mann-Whitney-Wilcoxon test will be less reliable.\n\tPlease, check graphically the samples distributions.\n")
          } else {
            cat("6) Kolmogorov-Smirnov test (ks.test()) -  The samples have the same distribution. min(p-value) : ",ks_result,"\n\tGood accuracy expected on the tests of Kruskal-Wallis and Mann-Whitney-Wilcoxon\n")
          }
          if (pvals3 <= pval) {
            cat("7) Kruskal-Wallis test (kruskal.test()) - At least one sample appears to show a difference. p-value:",pvals3,"\n")
          } else {
            cat("7) Kruskal-Wallis test (kruskal.test()) - No different sample a priori. p-value:",pvals3,"\n")
          }
        }
        if ((return==TRUE) | (verbose==TRUE)) {
          if (pvals<=pval) {
            pvals <- med1way(data~cat)$p.value
            if (is.na(pvals)) {
              if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.\n")}
            } else {
              if (pvals <= pval) {
                if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Significant differences between the median of samples. p-value:",pvals,"\n")
                  if (pvals3 > pval) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on medians give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(TRUE)}
              } else if (pvals > pval) {
                if (verbose==TRUE) {cat("8) Oneway ANOVA of medians (med1way()) - Non-significant differences between the median of samples. p-value:",pvals,"\n")
                  if (pvals3 <= pval) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on medians give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(FALSE)}
              }
            }
            if (code==TRUE){cat("pairwise.wilcox.test(data,cat,p.adjust.method='BH') #7)\n")}
            if (return==TRUE) {
              if (verbose==TRUE) {cat("9) Wilcoxon-Mann-Whitney test (pairwise.wilcox.test()) \n")}
              synth <- pairwise(data,cat,type="median",pval=pval,control=control,boot=boot,conf=conf,iter=iter)
              if ((verbose==TRUE) & (boot==TRUE) & any(synth$bootstrap$groups[,2]!=synth$groups[,2])) {
                cat("\tWarning! Bootstrap detects weaknesses in the significance of the results.\n")
              }
              return(synth)
            } else {
              if (pvals3 <= pval) {
                return(TRUE)
              } else {return(FALSE)}
            }
          } else {
            pvals <- t1way(data~cat)$p.value
            if (is.na(pvals)) {
              if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Failed, return NA. The Kruskal-Wallis test should be used for interpretation.\n")}
            } else {
              if (pvals <= pval) {
                if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Significant differences between the trimmed samples. p-value:",pvals,"\n")
                  if (pvals3 > pval) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on trimmed means give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(TRUE)}
              } else if (pvals > pval) {
                if (verbose==TRUE) {cat("8) One-way ANOVA on trimmed means (t1way()) - Non-significant differences between the trimmed samples. p-value:",pvals,"\n")
                  if (pvals3 <= pval) {
                    cat("\tWarning! The Kruskal-Wallis test and anova on trimmed means give contradictory results.\n")
                  }
                }
                if (return==FALSE) {return(FALSE)}
              }
            }
            if (code==TRUE){cat("library(WRS2) #7a)\nlincon(data~cat) #7b)\n")}
            if (return==TRUE) {
              if (verbose==TRUE) {cat("9) Correspondoncance post-hoc on trimmed means (lincon())\n")}
              synth <- pairwise(data,cat,type="lincon",pval=pval,control=control)
              return(synth)
            } else {
              if (pvals3 <= pval) {
                return(TRUE)
              } else {return(FALSE)}
            }
          }
        } else if (return==FALSE) {
          if (pvals3 <= pval) {return(TRUE)
          }else {return(FALSE)}
        }
      }
    }
  }
}

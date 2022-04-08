#' Validating a bootstrap linear regression model
#'
#' @param reg Linear model
#' @param optional A data.frame is formula is complex with I() function
#' @param plot Enable or disable the display of graphical analysis
#' @param verbose Enable or disable the display of the commented analysis
#' @param conf.level Confidence level for validation of the model
#' @param pval Minimal value accepted for validation of the model and his coefficients
#' @param iter Number of iterations
#'
#' @details bootreg allows to validate a model by bootstrap.
#' @details It will draw several times in the values of the model in order to test its robustness.
#' @details As an output, we can have a validation or not of the model (analysis=F), or we can have a table expressing the variability of the coefficients, their p-values and their maximum fluctuation values.
#' @return A graph of analysis.
#' @export
#'
#' @examples
#' data(mtcars);
#' corrigraph(mtcars);
#' reg<- lm(cyl~disp+hp,data=mtcars);
#' bootreg(reg, verbose=TRUE, plot=TRUE)
bootreg <- function(reg,data=c(),plot=TRUE,verbose=TRUE,conf.level=0.95,pval=0.05,iter=1000) {
  if (length(data)>1) {
	dt <- data
  } else {
	dt <- reg$model
  }
  numind <- nrow(dt)
  indices <- c(1:numind ) # num if individus
  enregistrement <- 0
  erreur <- c() ; predictions <- c() ; verity <- c()
  for (i in c(1:iter)) {
    indices_training <- sample(indices,size=numind ,replace=T)
    indices_test <- setdiff(indices ,indices_training)
    if (length(indices_test) > 0) {
      training  <- dt[indices_training,]	 
      test   <- 	 dt[indices_test,]
      #formula_ <- eval(reg$call[[2]])
	  formula_ <- formula(reg)
      oldw <- getOption("warn")
      options(warn = -1)
      #reg1 <- lm(formula=formula(formula),data=training,subset=indices_training)
	  reg1 <- lm(formula=formula(formula_),data=training)
      #reg1 <- update(reg1,subset=indices_training)
      options(warn = oldw)
      if (enregistrement == 0) {
        pval_mdl <- pf(summary(reg1)$fstatistic[1],summary(reg1)$fstatistic[2],summary(reg1)$fstatistic[3],lower.tail=FALSE)
        names(pval_mdl) <- "p-value of model"
        if (length(c(pval_mdl,summary(reg1)[[4]][,4])) == (length(reg$coefficients)+1)) {
          enregistrement <- enregistrement + 1
          p_values <- c(pval_mdl,summary(reg1)[[4]][,4])  # p_values correspond aux valeurs Pr récupérés pour chaque variable et chaque cycle du bootstrap (ici 1000)
          # Pour extraction des p-values sur glm() : mettre cela : coef(summary(model2))[,4]
          coeff <- reg1$coefficients
          oldw <- getOption("warn")
          options(warn = -1)
          predictions <- c(predictions,predict(reg1,dt)[indices_test])
          options(warn = oldw)
          #verity <- c(verity,test[,1])
		  verity <- c(verity,test[[names(get_all_vars(formula(reg$terms),dt))[1]]])
		  
        }
      } else {
        pval_mdl <- pf(summary(reg1)$fstatistic[1],summary(reg1)$fstatistic[2],summary(reg1)$fstatistic[3],lower.tail=FALSE)
        if (length(c(pval_mdl,summary(reg1)[[4]][,4])) == (length(reg$coefficients)+1)) {
          p_values <- rbind(p_values,c(pval_mdl,summary(reg1)[[4]][,4]))
          coeff <- rbind(coeff,reg1$coefficients) # Coeff correspond aux coefficients récupérés pour chaque variable et chaque cycle du bootstrap (ici 1000)
          oldw <- getOption("warn")
          options(warn = -1)
          predictions <- c(predictions,predict(reg1,dt)[indices_test])
          options(warn = oldw)
          #verity <- c(verity,test[,1])
		  verity <- c(verity,test[[names(get_all_vars(formula(reg$terms),dt))[1]]])
        }
      }
      # Test sur les predictions
      #predictions =c() ; verity = c()

    }}
  coeff <- na.omit(coeff) ; p_values <- na.omit(p_values)
  predverity <- data.frame(predictions ,verity); predverity <- na.omit(predverity)
  predictions  <- predverity[,1] ; verity  <- predverity[,2]
  confiance <- function(x,conf.level=0.99) { # seuil
    temp = sort(x) ; valeur_seuil = round(length(x)*conf.level)
    temp <- temp[valeur_seuil]
    return(temp)
  }
  mode <- function(x) {
    densite <- density(x)
    mode <- densite$x[which(densite$y==max(densite$y))]
    return(mode)
  }
  apply(coeff,2,median) -> coeff_median
  if (plot==T) {
    boxplot_Pr <- function(x,main="") {
      my_min <- min(c(apply(x,2,quantile)[2,]),0.0009)
      if (my_min <=0) {my_min <- 1e-20}
      boxplot(x,log="y",ylim=c(my_min,1),main=main)
      abline(h=0.05,col="red",lwd=2)
      abline(h=0.01,col="orange",lwd=2)
      abline(h=0.001,col="green",lwd=2)
    }
    hist_decrypt <- function(x,breaks=20,main="",sub="",conf.level=0.95) {
      densite <- density(x,na.rm=T) # créer une liste des densités
      hist(x,freq=F,col="#AAFFAA",ylim=c(0,max(densite$y)),breaks=breaks,main=main,cex=0.5,sub=sub) # il faut que freq=F
      lines(densite, col = "red",lwd=3) # Superposer la ligne
      abline(v=confiance(x,conf.level),col="black",lwd=5)
      abline(v=mean(x),col="red",lwd=3)
      abline(v=median(x),col="orange",lwd=2)
      abline(v=mode(x),col="green",lwd=1)
      legend("topright",col=c("black","red","orange","green"),
             c("p-value max in the confidence interval","mean","median","mode"),
             lwd=c(5,3,2,1))
    }
    layout(matrix(1:4,2,2))
    # p-value of the model
    #return(p_values)
    hist_decrypt(p_values[,1],sub="If the repartition is uniform:\ninsignificant model",breaks=20,
                 main="Distribution of the p-values of the model",conf.level=conf.level)
    # p-value of the model & coefficients
    boxplot_Pr(p_values,main="Distribution of the p-values of\nthe model and its coefficients")
    apply(coeff,2,function(x) {(x-median(x))/median(x)*100}) -> percent_coeff
    boxplot(percent_coeff,main="Fluctuation of coefficients (in %)",ylog=NA)
    abs(predictions-verity)-> temp
    by(temp,verity,confiance,conf.level=conf.level) -> CONF
    by(temp,verity,mean) -> MOY
    x <- as.numeric(names(MOY))
    plot(as.numeric(names(MOY)),CONF,type="l",lwd=3,col="red",
         xlim=c(min(x),max(x)),ylim=c(0,max(CONF)),
         xlab="Experimental values",ylab="Predictions",
         main="Average prediction error and\nmaximum error in the confidence interval")
    points(as.numeric(names(MOY)),MOY,type="l",lwd=2,col="black")
  }
  apply(p_values,2,median) -> p.values_median
  apply(p_values,2,confiance,conf.level=conf.level) -> p.values_max
  coeff_model <- c(NA,reg$coefficients)
  coeff_median <- c(NA,coeff_median)
  coeff_IC <- c(NA,apply(coeff,2,int.ech,conf.level=conf.level))
  synth <- rbind(p.values_median,p.values_max,coeff_model,coeff_median,coeff_IC)
  synth <- t(synth)
  synth <- data.frame(synth) ; synth <- cbind(data.frame(rownames(synth)),synth)
  rownames(synth)[1] <- "Model"
  if (verbose==T){return(synth)
  } else {
    if (max(synth[,2])>pval) {return(FALSE)
    } else {return(TRUE)}
  }
}

#' Validating a bootstrap linear regression model
#'
#' @param reg linear model.
#' @param data a data.frame of values if formula is complex with I() function (Optional).
#' @param plot enable or disable the display of graphical analysis.
#' @param verbose enable or disable the display of the commented analysis.
#' @param conf.level confidence level for validation of the model.
#' @param alpha minimal value accepted for validation of the model and his coefficients.
#' @param iter number of iterations.
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
bootreg <- function(reg,data=c(),plot=TRUE,verbose=TRUE,conf.level=0.95,alpha=0.05,iter=1000) {
  if (length(data)>1) {
	dt <- data
  } else {
	dt <- reg$model
  }
  # Contrôler la présence de $
  if (any(unlist(sapply(formula(reg), function(x) grepl("\\$", x))))) {
	options(warn = 0)
    warning("Warning! The character '$' has been detected in the formula. Prefer the syntax lm(Y~x, data=data) over lm(data$Y~data$x), which may cause bootreg() to fail.")
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
      #reg1 <- lm(formula=formula(formula_),data=training,subset=indices_training)
	  #reg1 <- lm(formula = formula(formula_), data = training, subset = (1:nrow(training)) %in% indices_training)
	  reg1 <- lm(formula=formula(formula_),data=training)
      #reg1 <- update(reg1, subset = (1:nrow(training)) %in% indices_training)
      options(warn = oldw)
      if (enregistrement == 0) {
        pval_mdl <- pf(summary(reg1)$fstatistic[1],summary(reg1)$fstatistic[2],summary(reg1)$fstatistic[3],lower.tail=FALSE)
        names(pval_mdl) <- "p-value of model"
        if (length(c(pval_mdl,summary(reg1)[[4]][,4])) == (length(reg$coefficients)+1)) {
          enregistrement <- enregistrement + 1
          p_values <- c(pval_mdl,summary(reg1)[[4]][,4])  # p_values correspond aux valeurs Pr récupérés pour chaque variable et chaque cycle du bootstrap (ici 1000)
          # Pour extraction des p-values sur glm() : mettre cela : coef(summary(model2))[,4]
          coeff <- reg1$coefficients
		  #print(coeff)
          oldw <- getOption("warn")
          options(warn = -1)
          predictions <- c(predictions,predict(reg1,dt)[indices_test])
          options(warn = oldw)
          #verity <- c(verity,test[,1])
		  if (length(all.vars(reg$terms))==length(names(get_all_vars(formula(reg$terms),dt)))){
			verity <- c(verity,test[[names(get_all_vars(formula(reg$terms),dt))[1]]])
		  } else {
			Ynames <- trimws(strsplit(deparse(reg$call[[2]]), "[~+]")[[1]][1])
			verity <- c(verity,test[[Ynames]])
		  }

        }
      } else {
        pval_mdl <- pf(summary(reg1)$fstatistic[1],summary(reg1)$fstatistic[2],summary(reg1)$fstatistic[3],lower.tail=FALSE)
        if (length(c(pval_mdl,summary(reg1)[[4]][,4])) == (length(reg$coefficients)+1)) {
          p_values <- rbind(p_values,c(pval_mdl,summary(reg1)[[4]][,4]))
          coeff <- rbind(coeff,reg1$coefficients) # Coeff correspond aux coefficients récupérés pour chaque variable et chaque cycle du bootstrap (ici 1000)
		  #print(reg1$coefficients)
          oldw <- getOption("warn")
          options(warn = -1)
		  # Prédiction sur les données test (not training)
          predictions <- c(predictions,predict(reg1,dt)[indices_test])
          options(warn = oldw)
          #verity <- c(verity,test[,1])
		  # Données réelles
		  if (length(all.vars(reg$terms))==length(names(get_all_vars(formula(reg$terms),dt)))){
			verity <- c(verity,test[[names(get_all_vars(formula(reg$terms),dt))[1]]])
		  } else {
			Ynames <- trimws(strsplit(deparse(reg$call[[2]]), "[~+]")[[1]][1])
			verity <- c(verity,test[[Ynames]])
		  }
		#cat("npred",nrow(predictions),"nverity",nrow(verity))
		#cat(length(predictions),"\n")
		#	if (length(predict(reg1,dt)[indices_test])!=length(test[[names(get_all_vars(formula(reg$terms),dt))[1]]])) {
		#		cat("prediction\n")
		#		print(predict(reg1,dt)[indices_test])
		#		cat("réalité\n")
		#		print(test[[names(get_all_vars(formula(reg$terms),dt))[1]]])
		#		print(names(get_all_vars(formula(reg$terms),dt)))
		#		break
		#	} 
			#else {
		#		print(names(get_all_vars(formula(reg$terms),dt)))
		#		break
		#	}
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
    layout(matrix(1:3,1,3))
    # p-value of the model
    #return(p_values)
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
         main="Average prediction error (blackà and\nmaximum error in the confidence interval (red)")
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
    if (max(synth[,2])>alpha) {return(FALSE)
    } else {return(TRUE)}
  }
}

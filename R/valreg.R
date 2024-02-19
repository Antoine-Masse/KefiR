#' Function to validate a regression model
#'
#' @param reg A regression model
#' @param verbose To see the detailed balance sheet
#' @param nvar The maximum number of variables allowed.
#' @param boot For checking the model by bootstrap with bootreg() (FALSE or TRUE)
#' @param alpha Maximum value accepted for the p-values of the model and its coefficients.
#' @param conf.level Confidence interval accepted to validate the regression model by bootstrap.
#' @param plot For seeing the graphical analysis of bootreg() (FALSE or TRUE)
#' @param data optional, the data.frame of data if complex model.
#' @param raintest_alpha Minimal value of p-value accepted for Rainbow test
#' @param dwtest_alpha Minimal value of p-value accepted for Durbin-Watson test
#' @param shapiro_alpha Minimal value of p-value accepted for Shapiro-Wilk test
#' @param bptest_alpha Minimal value of p-value accepted for Breush-Pagan test
#'
#' @return This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).
#' @return valreg will therefore validate the regression model, control the p-values and, possibly (boot argument), control the reliability by bootstrap with the bootreg function.
#' @import lmtest
#' @import tseries
#' @importFrom stats cooks.distance
#' @export
#'
#' @examples
#' #Example 1
#' data(iris)
#' reg <- lm(Sepal.Length~.,data=iris[,1:4])
#' valreg(reg,verbose=TRUE)
#' #Example 2
#' data(airquality)
#' reg <- lm(Wind~Temp,data=airquality)
#' valreg(reg,verbose=TRUE)
#' #Example 3
#' data(mtcars)
#' reg <- lm(cyl~gear,data=mtcars)
#' valreg(reg,verbose=TRUE,boot=TRUE)
#' #Example 4
#' data(iris)
#' reg <- lm(Petal.Width~Petal.Length,data=iris[iris$Species=="virginica",])
#' valreg(reg,verbose=TRUE,boot=TRUE,plot=TRUE)
#' plot(reg)
#' #Example 5
#' data(mtcars);
#' reg<- lm(cyl~disp+hp,data=mtcars);
#' corrigraph(mtcars,"cyl");
#' valreg(reg, verbose=TRUE, plot=TRUE)
valreg <- function(reg,verbose=TRUE,nvar=5,boot=TRUE,alpha=0.05,conf.level=0.95,
                   plot=FALSE,data=c(),
                   raintest_alpha=0.05,dwtest_alpha=0.03,shapiro_alpha=0.05,bptest_alpha=0.05) {
  error <- TRUE
  nvar1 <- length(coef(reg))
  if (nvar1 > nvar) {if(verbose==TRUE){cat("More than ",nvar," variables\n")};error = FALSE
  } else if (length(reg$fitted.values)<(nvar+1)) {if(verbose==TRUE){cat("Not enough values in the subset\n")};error = FALSE
  } else if ( (length(summary(reg)[[4]][,4])) != (length(coef(reg))) ) {error = FALSE
  } else {
	# p-values of the model
	if (verbose==TRUE) {cat("01- Analysis of the p-values of the model and its coefficients.\n")}
    pval_mdl <- pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)
    if (pval_mdl > alpha) {if(verbose==TRUE){cat("\tWarning!\n\tBad significance of the model. p-value:",pval_mdl,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tGood significance of the model. p-value:",pval_mdl,"\n")}}
    pval_coeff <- summary(reg)[[4]][,4]
    if (max(pval_coeff) > alpha) {if(verbose==TRUE){cat("\tWarning!\n\tBad significance of the coefficients. max(p.value) :",max(pval_coeff),"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tGood significance of the coefficients. max(pval_coeff) :",max(pval_coeff),"\n")}}
	# Rainbow test - Adequacy
	if (verbose==TRUE) {cat("02- Analysis of the adequacy of model (Equivalence between the global model and the model established on the best points.).\n")}
	# Vérifier si l'intercept est inclus dans le modèle
	if ("(Intercept)" %in% names(coef(reg))) {
		nvar2 <- length(coef(reg)) - 1
	} else {
		nvar2 <- length(coef(reg))
	}
	if (nvar2 > 1) {
		raintest(reg,order.by="mahalanobis")$p.value -> pvalt # test de rainbow Adequacy
		if (pvalt < raintest_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tRainbow test ordered by mahalanobis (raintest()) - Bad adequacy. p.value : ",pvalt,"\n")};error = FALSE
	#} else if (pvalt2 < raintest_pval){if(verbose==TRUE){cat("\tWarning!\n\tRainbow test ordered by Cook's distance (raintest()) - Bad adequacy. p.value : ",pvalt2,"\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tRainbow test (raintest()) - Good adequacy. p.value : ",pvalt,"\n")}}
	} else {
		if (length(data)==0) {
			data <- reg$model
		}
		z <- as.formula(paste0("~",attr (reg$terms, "term.labels")))
		raintest(reg,order.by=z, data=data)$p.value -> pvalt # test de rainbow Adequacy
		if (pvalt < raintest_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tRainbow test ordered by X(raintest()) - Bad adequacy. p.value : ",pvalt,"\n")};error = FALSE
		#} else if (pvalt2 < raintest_pval){if(verbose==TRUE){cat("\tWarning!\n\tRainbow test ordered by Cook's distance (raintest()) - Bad adequacy. p.value : ",pvalt2,"\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tRainbow test (raintest()) - Good adequacy. p.value : ",pvalt,"\n")}}		
	}
	#cooks.distance(reg)->cooksd
	#raintest(reg,order.by=cooksd)$p.value -> pvalt2
    
	# Durbin-Watson test
	if (verbose==TRUE) {cat("03- Analysis of independence of the residuals\n : (Warning : only by sorting the values as a function of Y. Check manually for others criteria.)\n")}
	if (length(data)>0) {ypred <- predict(reg, data)
	}else{ypred <- predict(reg, reg$model)}
    try(dwtest(reg,ypred)$p.value) -> pvalt # Independence of DurbinWatson residuals
    if (pvalt < dwtest_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tDurbin-Watson test (dwtest()) - Bad independence of the residuals. p.value : ",pvalt,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tDurbin-Watson test (dwtest()) - Good independence of the residuals. p.value : ",pvalt,"\n")}}
	# Shapiro-test : distribution of residuals
	if (verbose==TRUE) {cat("04- Analysis of distribution of residuals.\n")}
	if (length(reg$residuals)<=50) {
		shapiro.test(residuals(reg))$p.value->pvalt # Normal distribution of residuals
		if (pvalt < shapiro_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tShapiro-Wilk test (shapiro.test()) - Non-normal distribution of residuals. p.value : ",pvalt,"\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tShapiro-Wilk test (shapiro.test()) - Normal distribution of residuals. p.value : ",pvalt,"\n")}}
	} else if (length(reg$residuals)<=1000) {
		cat("Warning ! no Shapiro-Wilk test because more than 50 values.\n")
		jarque.bera.test(residuals(reg))$p.value->pvalt
		if (pvalt < shapiro_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tJarque-Bera test (jarque.bera.test() of {tseries}) -  - Non-normal distribution of residuals. p.value : ",pvalt,"\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tJarque-Bera test (jarque.bera.test() of {tseries}) - Normal distribution of residuals. p.value : ",pvalt,"\n")}}
	} else {
		cat("Warning ! no Shapiro-Wilk test or Jarque-Bera test because more than 1000 values.\n")
		ks.test(residuals(reg), "pnorm", mean = mean(residuals(reg)), sd = sd(residuals(reg)))$p.value->pvalt
		if (pvalt < shapiro_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tKolmogorov-Smirnov test (ks.test()) -  - Non-normal distribution of residuals. p.value : ",pvalt,"\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tKolmogorov-Smirnov test (ks.test()) - Normal distribution of residuals. p.value : ",pvalt,"\n")}}
	}
	# Variance of residuals
	white <- function(model) {
	  if (any(is.na(model$coefficients))) {stop("NA in coefficients.")}
	  data = model$model
	  # Extraire les termes du modèle initial
	  terms <- attr (model$terms, "term.labels")
	  # Extraire la formule du modèle initial
	  form <- formula (model)
	  # Extraire le nom du jeu de données
	  data_name <- getCall (model)$data
	  # Créer la formule du modèle augmenté
		if (length (terms) > 1) { # Si le modèle initial contient plus d'une variable explicative
		  formula <- paste0 (deparse(form), "+", paste (sapply (terms, function (x) paste0 ("I(", x, "^2)")), collapse = "+"), "+", paste (combn (terms, 2, function (x) paste0 ("I(", x[1], "*", x[2], ")")), collapse = "+"))
		} else { # Si le modèle initial contient une seule variable explicative
		  formula <- paste0 (deparse(form), "+", paste0 ("I(", terms, "^2)"))
		}
	  # Faire le test de White
		model <- lm(as.formula (formula), data = data)
		return(bptest (model)) # Utiliser l'argument data
	}
	if (verbose==TRUE) {cat("05- Analysis of variance of residuals.\n")}
    if (length(coef(reg))>=2 & length(reg$residuals)<=1000) {
      bptest(reg)$p.value -> pvalt # Breush Pagan: constant variance of residuals
      if (pvalt < bptest_alpha) {if(verbose==TRUE){cat("\tWarning!\n\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value : ",pvalt,"\n")};error = FALSE
	  } else {if(verbose==TRUE){cat("\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value : ",pvalt,"\n")}}
    } else if (length(reg$residuals)>1000) {
		cat("\tWarning!\n\tToo many values to justify checking the constancy of the variance. Breush-Pagan test bptest() of {lmtest} may be hypersensitive.\n")
	}
	# Cooks's distance Leverage effect
	cooks.distance(reg)->cooksd
	if (verbose==TRUE) {cat("06- Analysis of leverage effect.\n")}
    if (max(cooksd,na.rm=TRUE) > 1) {if(verbose==TRUE){cat("\tWarning!\n\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance())",max(cooksd,na.rm=TRUE),"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tCook's distance (cooks.distance()) - No leverage effect. max(cooks.distance())",max(cooksd,na.rm=TRUE),"\n")}}
    if (boot == TRUE) {
		if (verbose==TRUE) {cat("07- Analysis of solidity of model by boostrap.\n")}
		if (length(data)>0) {bootreg(reg,verbose=FALSE,alpha=alpha,conf.level=conf.level,plot=plot,data=data) -> bootres
		}else{bootreg(reg,verbose=FALSE,alpha=alpha,conf.level=conf.level,plot=plot) -> bootres}
		if (bootres==FALSE) {if(verbose==TRUE){cat("\tWarning!\n\tBootstrap (bootreg()) - Fragility of the model in boostrap. Please, use bootreg()\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tBootstrap (bootreg()) - Solidity of the model in boostrap.\n")}}
    }
  }
  return(error)
}


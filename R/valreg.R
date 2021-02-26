#' Function to validate a regression model
#'
#' @param reg A regression model
#' @param verbose To see the detailed balance sheet
#' @param nvar The maximum number of variables allowed.
#' @param boot For checking the model by bootstrap with bootreg() (FALSE or TRUE)
#' @param pval Maximum value accepted for the p-values of the model and its coefficients.
#' @param conf.level Confidence interval accepted to validate the regression model by bootstrap.
#' @param plot For seeing the graphical analysis of bootreg() (FALSE or TRUE)
#' @param raintest_pval Minimal value of p-value accepted for Rainbow test
#' @param dwtest_pval Minimal value of p-value accepted for Durbin-Watson test
#' @param shapiro_pval Minimal value of p-value accepted for Shapiro-Wilk test
#' @param bptest_pval Minimal value of p-value accepted for Breush-Pagan test
#'
#' @return This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).
#' @return valreg will therefore validate the regression model, control the p-values and, possibly (boot argument), control the reliability by bootstrap with the bootreg function.
#' @import lmtest
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
valreg <- function(reg,verbose=TRUE,nvar=5,boot=TRUE,pval=0.05,conf.level=0.95,
                   plot=TRUE,
                   raintest_pval=0.05,dwtest_pval=0.03,shapiro_pval=0.05,bptest_pval=0.05) {
  if (plot==T) {boot<-TRUE}
  error <- TRUE
  nvar1 <- length(coef(reg))
  if (nvar1 > nvar) {if(verbose==TRUE){cat("More than ",nvar," variables\n")};error = FALSE
  } else if (length(reg$fitted.values)<(nvar+1)) {if(verbose==TRUE){cat("Not enough values in the subset\n")};error = FALSE
  } else if ( (length(summary(reg)[[4]][,4])) != (length(coef(reg))) ) {error = FALSE
  } else {
	# p-values of the model
	if (verbose==TRUE) {cat("01- Analysis of the p-values of the model and its coefficients.\n")}
    pval_mdl <- pf(summary(reg)$fstatistic[1],summary(reg)$fstatistic[2],summary(reg)$fstatistic[3],lower.tail=FALSE)
    if (pval_mdl > pval) {if(verbose==TRUE){cat("\tWarning!\n\tBad significance of the model. p-value:",pval_mdl,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tGood significance of the model. p-value:",pval_mdl,"\n")}}
    pval_coeff <- summary(reg)[[4]][,4]
    if (max(pval_coeff) > pval) {if(verbose==TRUE){cat("\tWarning!\n\tBad significance of the coefficients. max(pval_coeff) :",max(pval_coeff),"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tGood significance of the coefficients. max(pval_coeff) :",max(pval_coeff),"\n")}}
	# Rainbow test - Adequacy
	if (verbose==TRUE) {cat("02- Analysis of the adaquacy of model (Equivalence between the global model and the model established on the best points.).\n")}
	raintest(reg,order.by="mahalanobis")$p.value -> pvalt # test de rainbow Adequacy
    if (pvalt < raintest_pval) {if(verbose==TRUE){cat("\tWarning!\n\tRainbow test (raintest()) - Bad adequacy. p.value : ",pvalt,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tRainbow test (raintest()) - Good adequacy. p.value : ",pvalt,"\n")}}
	# Durbin-Watson test
	if (verbose==TRUE) {cat("03- Analysis of independence of the residuals.\n")}
    dwtest(reg)$p.value -> pvalt # Independence of DurbinWatson residuals
    if (pvalt < dwtest_pval) {if(verbose==TRUE){cat("\tWarning!\n\tDurbin-Watson test (dwtest()) - Bad independence of the residuals. p.value : ",pvalt,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tDurbin-Watson test (dwtest()) - Good independence of the residuals. p.value : ",pvalt,"\n")}}
	# Shapiro-test : distribution of residuals
	if (verbose==TRUE) {cat("04- Analysis of distribution of residuals.\n")}
    shapiro.test(residuals(reg))$p.value->pvalt # Normal distribution of residuals
    if (pvalt < shapiro_pval) {if(verbose==TRUE){cat("\tWarning!\n\tShapiro-Wilk test (shapiro.test()) - Non-normal distribution of residuals. p.value : ",pvalt,"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tShapiro-Wilk test (shapiro.test()) - Normal distribution of residuals. p.value : ",pvalt,"\n")}}
	# Variance of residuals
	if (verbose==TRUE) {cat("05- Analysis of variance of residuals.\n")}
    if (length(coef(reg))>=2) {
      bptest(reg)$p.value -> pvalt # Breush Pagan: constant variance of residuals
      if (pvalt < bptest_pval) {if(verbose==TRUE){cat("\tWarning!\n\tBreush-Pagan test (bptest()) - Non-constant variance of the residuals. p.value : ",pvalt,"\n")};error = FALSE
	  } else {if(verbose==TRUE){cat("\tBreush-Pagan test (bptest()) - Constant variance of the residuals. p.value : ",pvalt,"\n")}}
    }
	# Cooks's distance Leverage effect
	if (verbose==TRUE) {cat("06- Analysis of leverage effect.\n")}
    cooks.distance(reg)->cooksd
    if (max(cooksd,na.rm=TRUE) > 1) {if(verbose==TRUE){cat("\tWarning!\n\tCook's distance (cooks.distance()) - Leverage effect. max(cooks.distance())",max(cooksd,na.rm=TRUE),"\n")};error = FALSE
	} else {if(verbose==TRUE){cat("\tCook's distance (cooks.distance()) - No leverage effect. max(cooks.distance())",max(cooksd,na.rm=TRUE),"\n")}}
    if (boot == TRUE) {
		if (verbose==TRUE) {cat("07- Analysis of solidity of model by boostrap.\n")}
		bootreg(reg,verbose=FALSE,pval=pval,conf.level=conf.level,plot=plot) -> bootres
		if (bootres==FALSE) {if(verbose==TRUE){cat("\tWarning!\n\tBootstrap (bootreg()) - Fragility of the model in boostrap. Please, use bootreg()\n")};error = FALSE
		} else {if(verbose==TRUE){cat("\tBootstrap (bootreg()) - Solidity of the model in boostrap.\n")}}
    }
  }
  return(error)
}


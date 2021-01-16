#' Function to validate a regression model
#'
#' @param reg : A regression model
#' @param analyse : To see the detailed balance sheet
#' @param nvar : The maximum number of variables allowed.
#'
#' @return This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).
#' @import lmtest
#' @importFrom stats cooks.distance
#' @export
#'
#' @examples
#' #Example 1
#' data(iris)
#' reg <- lm(Sepal.Length~.,data=iris[,1:4])
#' valreg(reg,analyse=TRUE)
#' #Example 2
#' data(airquality)
#' reg <- lm(Wind~Temp,data=airquality)
#' valreg(reg,analyse=TRUE)
#' #Example 3
#' data(mtcars)
#' reg <- lm(cyl~gear,data=mtcars)
#' valreg(reg,analyse=TRUE)
#' #Example 4
#' data(iris)
#' reg <- lm(Petal.Width~Petal.Length,data=iris[iris$Species=="virginica",])
#' valreg(reg,analyse=TRUE)
#' plot(reg)
valreg <- function(reg,analyse=FALSE,nvar=5) {
  error <- "OK"
  nvar1 <- length(coef(reg))
  if (nvar1 > nvar) {if(analyse==TRUE){cat("Plus de 6 variables\n")};error = "error"
  } else if ( (length(summary(reg)[[4]][,4])) != (length(coef(reg))) ) {error = "error"
  } else {
    raintest(reg)$p.value -> pval # test de rainbow Adequacy
    if (pval < 0.05) {if(analyse==TRUE){cat("Bad adequacy.\n")};error = "error"}
    dwtest(reg)$p.value -> pval # Independence of DurbinWatson residues
    if (pval < 0.05) {if(analyse==TRUE){cat("Bad independence of the residues.\n")};error = "error"}
    shapiro.test(residuals(reg))$p.value->pval # Normal distribution of residues
    if (pval < 0.05) {if(analyse==TRUE){cat("Non-normal distribution of residues.\n")};error = "error"}
    if (length(coef(reg))>=2) {
      bptest(reg)$p.value -> pval # Breush Pagan: constant variance of residuals
      if (pval < 0.05) {if(analyse==TRUE){cat("Non-constant variance of the residuals.\n")};error = "error"}
    }
    cooks.distance(reg)->cooksd
    if (max(cooksd,na.rm=TRUE) > 1) {if(analyse==TRUE){cat("Leverage effect.\n")};error = "error"}
  }
  return(error)
}

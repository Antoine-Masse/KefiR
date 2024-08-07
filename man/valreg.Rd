% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/valreg.R
\name{valreg}
\alias{valreg}
\title{Function to validate a regression model}
\usage{
valreg(
  reg,
  verbose = TRUE,
  nvar = 5,
  boot = TRUE,
  alpha = 0.05,
  conf.level = 0.95,
  plot = FALSE,
  data = c(),
  raintest_alpha = 0.05,
  dwtest_alpha = 0.03,
  shapiro_alpha = 0.05,
  bptest_alpha = 0.05
)
}
\arguments{
\item{reg}{A regression model}

\item{verbose}{To see the detailed balance sheet}

\item{nvar}{The maximum number of variables allowed}

\item{boot}{For checking the model by bootstrap with bootreg() (FALSE or TRUE)}

\item{alpha}{Maximum value accepted for the p-values of the model and its coefficients}

\item{conf.level}{Confidence interval accepted to validate the regression model by bootstrap}

\item{plot}{For seeing the graphical analysis of bootreg() (FALSE or TRUE)}

\item{data}{optional, the data.frame of data if complex model.}

\item{raintest_alpha}{Minimal value of p-value accepted for Rainbow test}

\item{dwtest_alpha}{Minimal value of p-value accepted for Durbin-Watson test}

\item{shapiro_alpha}{Minimal value of p-value accepted for Shapiro-Wilk test}

\item{bptest_alpha}{Minimal value of p-value accepted for Breush-Pagan test}
}
\value{
This function allows to run all the tests necessary to validate a regression model (check the normal distribution of the residuals, avoid leverage effects, control the variance of the residuals...).

valreg will therefore validate the regression model, control the p-values and, possibly (boot argument), control the reliability by bootstrap with the bootreg function.
}
\description{
Function to validate a regression model
}
\examples{
# Example 1: Linear model
data(iris)
reg <- lm(Sepal.Length~.,data=iris[,1:4])
valreg(reg,verbose=TRUE)

# Example 2: Mixed model
library(lme4)
data(sleepstudy)
reg_mixed <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
valreg(reg_mixed, verbose=TRUE)

# Example 3: Mixed model
library(lme4)
data(iris)
reg_mixed <- lmer(Sepal.Length ~ Petal.Length + (1|Species), data = iris)
valreg(reg_mixed, verbose=TRUE)
}
\author{
Antoine MASSE (2024)
}

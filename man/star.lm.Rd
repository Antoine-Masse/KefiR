% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/star.R
\name{star.lm}
\alias{star.lm}
\title{Get stars for a linear model}
\usage{
\method{star}{lm}(p.val, ...)
}
\arguments{
\item{p.val}{A linear model object.}

\item{...}{Additional arguments (currently not used).}
}
\value{
A character string representing the significance stars.
}
\description{
Get stars for a linear model
}
\details{
Use the function \code{glance} from the package \code{broom} to get the p-value of
the linear model, and then convert it into stars with the function \code{star}.
}
\examples{
# Example usage
lm_model <- lm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris)
#star.lm(lm_model)
}

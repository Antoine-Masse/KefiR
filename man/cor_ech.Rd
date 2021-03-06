% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_bootstrap.R
\name{cor_ech}
\alias{cor_ech}
\title{Anticipate the sample size that will be needed to detect a correlation (if there is one)}
\usage{
cor_ech(x, y, iter = 100)
}
\arguments{
\item{x}{: The values of my x variable}

\item{y}{: The values of my x variable}

\item{iter}{: Number of iterations of bootstrapping}
}
\value{
: When I do a regression/correlation, I may have an insignificant p-value due to the fact that I don't have enough values. This function exercises a bootstrap that will suggest how many values it would take to have a significant correlation (if any).
}
\description{
Anticipate the sample size that will be needed to detect a correlation (if there is one)
}
\examples{
x <- c(1,5,6,9,10)
y <- c(0,4.5,10,9,30)
plot(x,y,cex=2,pch=16)
cor.test(x,y)$p.value # No significant
cor_ech(x,y)
}

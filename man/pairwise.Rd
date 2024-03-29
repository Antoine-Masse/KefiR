% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise.R
\name{pairwise}
\alias{pairwise}
\title{Automation function for pairwise calculations.}
\usage{
pairwise(
  x,
  g,
  type = "mean",
  alpha = 0.05,
  control = c(),
  pool.sd = FALSE,
  silent = TRUE,
  boot = FALSE,
  iter = 500,
  conf = 0.95
)
}
\arguments{
\item{x}{numerical vector}

\item{g}{category vector}

\item{type}{'mean' for pairwise.t.test(p.adjust.method="holm"), 'median' for pairwise.wilcox.test(p.adjust.method="BH"), 'ks' for ks.test(), 'lincon' for lincon() of {WSR2}}

\item{alpha}{threshold value of p-value to establish the groups.}

\item{control}{name of the category that will be used as a control to establish differences with '\emph{', '\strong{' and '}}'.}

\item{pool.sd}{switch to allow/disallow the use of a pooled SD.}

\item{silent}{for displaying or not warnings.}

\item{boot}{to activate the boostrap on 'mean' and 'median'.}

\item{iter}{number f iterations (boot==TRUE).}

\item{conf}{confidence level of bootstrap.}
}
\value{
This function automates the work of the ks.test(), lincon() functions of {WSR2}, pairwise.t.test() and pairwise.wilcox.test() and extracts groups of means or comparisons to a control with the catego() function.

It pre-sorts the means/medians to ensure that the groups are identified in ascending order.

It also identifies the robustness of these groups by establishing a bootstrap.
}
\description{
Automation function for pairwise calculations.
}
\examples{
data(iris)
pairwise(iris[,1],iris[,5],type="mean")# t.test
pairwise(iris[,1],iris[,5],type="median",alpha=0.01,boot=TRUE)#wilcox
pairwise(iris[,1],iris[,5],type="ks")
pairwise(iris[,1],iris[,5],type="lincon")
}

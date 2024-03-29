% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m.test.R
\name{m.test}
\alias{m.test}
\title{Automatic average/median/variance comparison function.}
\usage{
m.test(
  data,
  cat,
  alpha = 0.05,
  verbose = TRUE,
  return = TRUE,
  paired = FALSE,
  control = c(),
  maxcat = 50,
  plot = TRUE,
  silent = TRUE,
  boot = TRUE,
  iter = 500,
  conf = 0.95,
  code = FALSE,
  debug = FALSE
)
}
\arguments{
\item{data}{numerical vector}

\item{cat}{category vector}

\item{alpha}{p-value threshold value for all the tests.}

\item{verbose}{to display the full reasoning of the analysis.}

\item{return}{allows to return the results of pairwise analysis (p-values and groups).}

\item{paired}{(under development) to allow the analysis of matched data.}

\item{control}{name of the category that will eventually be used as a control.}

\item{maxcat}{maximum number of categories allowed. When this number is high, some tests may return an error message.}

\item{plot}{to display the distribution of the data.}

\item{silent}{for displaying or not warnings.}

\item{boot}{to activate the boostrap on 'mean' and 'median'.}

\item{iter}{number f iterations (boot==TRUE).}

\item{conf}{confidence level of bootstrap.}

\item{code}{allows to display the simplified R source code to be able to do the same R study step by step.}

\item{debug}{when m.test return error.}
}
\value{
m.test() runs a decision tree to choose the most appropriate test series for sample comparison.

She chooses the tests, justifies her choices.

It can output groups of means or a comparison to a control.

Finally, it will measure the robustness of the results by bootstrap.
}
\description{
Automatic average/median/variance comparison function.
}
\examples{
data(iris)
m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE)
m.test(iris[1:100,1],iris[1:100,5],verbose=TRUE, return=TRUE)
m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE)
m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE)
m.test(iris[,4],iris[,5],verbose=TRUE, plot=FALSE, return=FALSE, boot=FALSE)
m.test(iris[,1],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
m.test(iris[,2],iris[,5],verbose=TRUE, return=TRUE,control="virginica")
m.test(iris[,3],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
m.test(iris[,4],iris[,5],verbose=TRUE, return=TRUE,control="setosa")
}

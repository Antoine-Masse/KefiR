% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrigraph.R
\name{corrigraph}
\alias{corrigraph}
\title{igraph of correlated variables global or in relation to y}
\usage{
corrigraph(
  data,
  colY = c(),
  colX = c(),
  type = "x",
  alpha = 0.05,
  exclude = c(0, 0, 0),
  ampli = 4,
  return = FALSE,
  wash = "stn",
  multi = TRUE,
  mu = FALSE,
  prop = FALSE,
  layout = "fr",
  cluster = TRUE,
  verbose = FALSE,
  NAfreq = 1,
  NAcat = FALSE,
  level = 2,
  evolreg = FALSE
)
}
\arguments{
\item{data}{a data.frame}

\item{colY}{a vector of indices or variables to predict. To force the correlogram to display only the variables correlated to a selection of Y.}

\item{colX}{a vector of indices or variables to follow. We will only keep the variables that are connected to them on 1 or more levels (level parameter).}

\item{type}{"x" or "y". To force the display in correlogram mode (colX, type = "x") or in prediction mode (colY, type = "y").}

\item{alpha}{the maximum permissible p-value for the display}

\item{exclude}{the minimum threshold of displayed correlations - or a vector of threshold in this order : c(cor,mu,prop)}

\item{ampli}{coefficient of amplification of vertices}

\item{return}{if return=TRUE, returns the correlation matrix of significant correlation.}

\item{wash}{automatically eliminates variables using differents methods when there are too many variables (method = NA, stn (signal-to-noise ratio), sum, length).}

\item{multi}{to ignore multiple regressions and control only single regressions.}

\item{mu}{to display the effect on median/mean identified by m.test().}

\item{prop}{to display the dependencies between categorical variables identified by GTest().}

\item{layout}{to choose the network organization method - choose "fr", "circle", "kk" or "3d".}

\item{cluster}{to make automatic clustering of variables or not.}

\item{verbose}{to see the comments.}

\item{NAfreq}{from 0 to 1. NA part allowed in the variables. 1 by default (100\% of NA tolerate).}

\item{NAcat}{TRUE or FALSE. Requires recognition of missing data as categories.}

\item{level}{to be used with colY. Number of variable layers allowed (minimum 2, default 5).}

\item{evolreg}{TRUE or FALSE. Not yet available. Allows you to use the evolreg function to improve the predictive ability (R squared) for the variables specified in colY.}
}
\value{
Depending on the parameters:
\describe{
\item{igraph}{A correlation graph network (igraph) of the variables of a data.frame. Non-numeric variables or missing data may be present. Vertices (circles) represent variables, with size indicating connectivity. The color of the edges reflects the nature of the correlation (positive in blue, negative in red). The width of the edges represents the strength of the correlation.}
\item{mu/prop}{If \code{mu} is TRUE or \code{prop} is specified: Connections display mean effects (orange) and dependencies between categorical variables (pink). The edge sizes depend on p-values from kruskal.test() and GTest().}
\item{Y specification}{When \code{colY} is specified: The correlogram identifies X variables correlated to Y, iterating through layers specified by \code{level}. X variables not related to Y are excluded.}
\item{vertex colors}{The color of vertices indicates significant correlations (blue for positive, red for negative, purple for both).}
\item{max predictive capacity}{Values displayed next to Y variables (\code{colY}) indicate the maximum predictive capacity by one or two variables.}
\item{correlation matrix}{If \code{return} is TRUE, the function returns the correlation matrix of significant correlations.}
}
}
\description{
igraph of correlated variables global or in relation to y
}
\examples{
# Example 1
data(swiss)
corrigraph(swiss)
# Example 2
data(airquality)
corrigraph(airquality,layout="3d")
# Example 3
data(airquality)
corrigraph(airquality,c("Ozone","Wind"),type="y")
# Example 4
data(iris)
corrigraph(iris,mu=TRUE)
# Example 5
require(MASS) ; data(Aids2)
corrigraph(Aids2 ,prop=TRUE,mu=TRUE,exclude=c(0.3,0.3,0))
# Example 6
data(airquality)
corrigraph(airquality,c("Ozone","Wind"),type="x")
}

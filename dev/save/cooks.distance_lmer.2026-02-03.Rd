% Generated manually for package check
\name{cooks.distance_lmer}
\alias{cooks.distance_lmer}
\title{Compute Cook's distance for mixed models}
\usage{
cooks.distance_lmer(model, obs = FALSE)
}
\arguments{
\item{model}{A mixed-effects model of class \code{lmerMod} (fitted with lme4).}
\item{obs}{Logical. If TRUE, computes Cook's distance for individual observations as well.}
}
\value{
A named list with Cook's distance values for each random-effect grouping factor,
and optionally for each observation.
}
\description{
Compute Cook's distance for mixed models fitted with \code{lme4}.
}
\examples{
\dontrun{
library(lme4)
data(sleepstudy)
model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
cd <- cooks.distance_lmer(model)
}
}

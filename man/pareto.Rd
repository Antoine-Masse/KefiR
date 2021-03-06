% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pareto.R
\name{pareto}
\alias{pareto}
\title{Pareto chart}
\usage{
pareto(
  x,
  names.arg = c(),
  bar.col = "cyan",
  line.col = "red",
  pch = 16,
  h = 80,
  h.lty = 3,
  main = "",
  xlab = "Defauts",
  ylab = "Frequence (\%)",
  ylab2 = "Cumul",
  mar = c(5, 4, 3, 4)
)
}
\arguments{
\item{x}{A vector}

\item{names.arg}{Vector description}

\item{bar.col}{Color of the bars of the barplot}

\item{line.col}{Color of the cumul line}

\item{pch}{Color of points of the cumul line}

\item{h}{Horizontal value (80 percent by default)}

\item{h.lty}{lty of horizontal line}

\item{main}{Title}

\item{xlab}{X title}

\item{ylab}{Y title}

\item{ylab2}{Y title of cumul}

\item{mar}{Marging}
}
\value{
For displaying a pareto chart
}
\description{
Pareto chart
}
\examples{
valeurs = c(20,10,12,5,2,80)
description = c("Para1","Para2","Para3","Para4","Para5","Para6")
pareto(valeurs,main="Diagramme de Pareto", names.arg=description, bar.col="blue")
}

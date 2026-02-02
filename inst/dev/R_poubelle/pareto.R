#' Pareto chart
#'
#' @param x A vector
#' @param names.arg Vector description
#' @param bar.col Color of the bars of the barplot
#' @param line.col Color of the cumul line
#' @param pch Color of points of the cumul line
#' @param h Horizontal value (80 percent by default)
#' @param h.lty lty of horizontal line
#' @param main Title
#' @param xlab X title
#' @param ylab Y title
#
#' @param ylab2 Y title of cumul
#' @param mar Marging
#'
#' @return For displaying a pareto chart
#' @export
#'
#' @examples
#' valeurs = c(20,10,12,5,2,80)
#' description = c("Para1","Para2","Para3","Para4","Para5","Para6")
#' pareto(valeurs,main="Diagramme de Pareto", names.arg=description, bar.col="blue")
pareto <- function(x, names.arg=c(), bar.col="cyan", line.col="red", pch=16, h=80, h.lty=3, main="", xlab="Defauts", ylab="Frequence (%)", ylab2="Cumul", mar=c(5,4,3,4)) {
  if (length(names.arg) > 0) {
    names.arg <- names.arg[order(x, decreasing = TRUE)]
  }
  x <- sort(x, decreasing = TRUE)
  x <- x * 100 / sum(x)
  cumul <- (cumsum(x) / sum(x)) * 100
  par(mar = mar)
  barplot(x, col = bar.col, axes = FALSE, ylim = c(0, 100), main = main, xlab = xlab, ylab = "", names.arg = names.arg)
  points(seq_along(x), cumul, pch = pch, col = line.col, xlab = "", ylab = "", type = "o")
  abline(h = h, lty = h.lty)
  box()
  axis(2)
  axis(4, c(0, 20, 40, 60, 80, 100), col.axis = line.col, col = line.col)
  mtext(ylab, side = 2, line = 2, cex = 1.2)
  mtext(ylab2, side = 4, col = "red", line = 2, cex = 1.2)

  result <- cbind(x, cumul)
  colnames(result) <- c("frequency", "cumul")
  if (length(names.arg) > 0) {
    rownames(result) <- names.arg
  }
  return(result)
}

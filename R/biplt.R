#' make a biplot that allows you to see the relationship between two components of a PCA, the weight of the respective variables and to color by category
#' 
#' @param pca a list with class "prcomp" from prcomp()
#' @param choices a vector indicating the indices of two desired components that we wish to relate
#' @param col a vector of color or a factor
#' @param threshold the elimination threshold to avoid displaying the variables that have the least influence on the components
#' @param xlab label of x
#' @param ylab label if y
#' @param cex length of points
#' @param pch pitch type
#' @param labels labels for the points
#' @param main title of the plot
#' @param legend logical, if TRUE display a legend
#' @param legend.labels vector of labels for the legend (if NULL, automatically deduced from col)
#' @param legend.pos position of the legend (default "topright")
#' @param legend.pch pch for the legend (if NULL, uses the same as pch)
#' @param legend.title title of the legend
#' @param legend.col vector of colors to use in legend (if NULL, automatically deduced)
#' @return a graphic biplot more esthetic of biplot() function
#'
#' @import utils
#' @importFrom plotrix draw.circle
#'
#' @export
#'
#' @examples
#' # Example 1
#' data(iris)
#' a <- prcomp(iris[,1:4])
#' biplt(a, col=iris$Species, threshold=0.8, choices=c(1,3), legend=TRUE)
#' 
#' # Example 2 with custom colors
#' cols <- rainbow(3)[as.numeric(iris$Species)]
#' biplt(a, col=cols, threshold=0.8, legend=TRUE, legend.labels=levels(iris$Species))
biplt <- function(pca, choices=c(1,2), col="blue", threshold=0.25,
                  xlab="cps", ylab="cps", cex=2, pch=16, 
                  labels=rownames(pca$x),
                  main="",
                  legend=FALSE,
                  legend.labels=NULL,
                  legend.pos="topright",
                  legend.pch=NULL,
                  legend.title="",
                  legend.col=NULL) {
  
  # Récupérer la description donnant la corrélation de chaque variable par composante
  poids_var_x <- pca$rotation[,choices[1]]
  poids_var_y <- pca$rotation[,choices[2]]
  vecnames <- rownames(pca$rotation)
  
  # OPTIONNEL - Calcul du poids global
  poids_var_global <- sqrt(poids_var_x^2 + poids_var_y^2)
  
  # FILTRER : les variables dont le poids est > threshold
  poids_var_x <- poids_var_x[poids_var_global > threshold]
  poids_var_y <- poids_var_y[poids_var_global > threshold]
  vecnames <- vecnames[poids_var_global > threshold]
  
  # Préparer les données pour la légende si nécessaire
  if (legend) {
    # Déterminer les couleurs et labels pour la légende
    if (is.factor(col)) {
      # Si col est un facteur, extraire les niveaux
      col_factor <- col
      if (is.null(legend.labels)) {
        legend.labels <- levels(col_factor)
      }
      if (is.null(legend.col)) {
        # Générer des couleurs pour chaque niveau
        n_levels <- length(levels(col_factor))
        legend.col <- rainbow(n_levels)
      }
      # Convertir le facteur en couleurs
      col <- legend.col[as.numeric(col_factor)]
    } else {
      # Si col est un vecteur de couleurs
      unique_cols <- unique(col)
      
      if (is.null(legend.col)) {
        legend.col <- unique_cols
      }
      
      if (is.null(legend.labels)) {
        # Si aucun label fourni, essayer de déduire ou utiliser les couleurs
        if (length(unique_cols) <= 20) {
          legend.labels <- paste("Groupe", 1:length(unique_cols))
        } else {
          legend.labels <- unique_cols
        }
      }
    }
    
    # Si legend.pch n'est pas fourni, utiliser le même que pch
    if (is.null(legend.pch)) {
      legend.pch <- pch
    }
  }
  
  # Tracer un fond blanc orthonormé
  plot(c(-1,1), c(-1,1), type="n", xlab="", ylab="", asp=1, axes=F, main="")
  grid()
  
  # Tracer le cercle
  draw.circle(0, 0, 1, lwd=3) 
  
  # Tracer les flèches
  arrows(0, 0, poids_var_x, poids_var_y, length=0.1, angle=17, lwd=2)
  
  # Calculer la position du nom de chaque vecteur
  poids_x_bin <- poids_var_x
  poids_x_bin[poids_x_bin < 0] <- -1
  poids_x_bin[poids_x_bin > 0] <- 1
  
  poids_y_bin <- poids_var_y
  poids_y_bin[poids_y_bin < 0] <- -1
  poids_y_bin[poids_y_bin > 0] <- 1
  
  poids_x <- poids_var_x^2
  poids_y <- poids_var_y^2
  poids_z <- poids_x + poids_y
  poids_x <- poids_x / poids_z * poids_x_bin * 0.08
  poids_y <- poids_y / poids_z * poids_y_bin * 0.08
  
  # Tracer les noms des vecteurs
  text(poids_var_x + poids_x, poids_var_y + poids_y, vecnames, cex=0.75, col="red")
  
  par(new = TRUE)
  
  # Graphique des données de PCA
  if (xlab == "cps") {xlab <- paste("Composante", choices[1])}
  if (ylab == "cps") {ylab <- paste("Composante", choices[2])}
  
  axeX <- pca$x[,choices[1]]
  axeY <- pca$x[,choices[2]]
  
  plot(axeX, axeY, pch=pch, cex=cex, col=col, xlab=xlab, ylab=ylab, main=main)
  grid()
  text(axeX, axeY - 0.05*axeY, labels, cex=cex)
  
  # Ajouter la légende si demandé
  if (legend) {
    legend(legend.pos, 
           legend=legend.labels, 
           col=legend.col, 
           pch=legend.pch,
           title=legend.title,
           cex=0.8,
           bg="white")
  }
  
  # Pour superposition des 2 graphiques
  par(new = TRUE)
  
  # Graphique du cercle et flèches par-dessus
  plot(c(-1,1), c(-1,1), type="n", xlab="", ylab="", asp=1, axes=F) 
  arrows(0, 0, poids_var_x, poids_var_y, length=0.1, angle=17, lwd=2)
  draw.circle(0, 0, 1, lwd=3)
  text(poids_var_x + poids_x, poids_var_y + poids_y, vecnames, cex=0.75, col="red")
}
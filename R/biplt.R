#' make a biplot that allows you to see the relationship between two components of a PCA, the weight of the respective variables and to color by category
#' pca,  =c(1,2), ="blue", =0.25,xlab="cps",ylab="cps",cex=2,pch
#' @param pca a list with class "prcomp" from prcomp()
#' @param choices a vector indicating the indices of two desired components that we wish to relate
#' @param col a vector of color
#' @param threshold the elimination threshold to avoid displaying the variables that have the least influence on the components
#' @param xlab label of x
#' @param ylab label if y
#' @param cex length of points
#' @param pch pitch type
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
#' biplt(a,col=iris$Species,threshold=0.8,choices=c(1,3))
biplt <- function(pca,  choices=c(1,2), col="blue", threshold=0.25,xlab="cps",ylab="cps",cex=2,pch=16,labels=rownames(pca$x)) {
	# Récupérer la description donnant la corrélation de chaque variable par composante
	# 1 et 2 pour composantes 1 et 2
	poids_var_x <- pca$rotation[,choices[1]]
	poids_var_y <- pca$rotation[,choices[2]]
	vecnames <- rownames(pca$rotation)
	#OPTIONNEL
	poids_var_global <- sqrt(poids_var_x^2+poids_var_y^2)
	# FILTRER : Ex : les variables dont le poids est > 0.25
	#barplot(poids_var_global,names.arg=vecnames) ; x11()
	threshold= threshold
	poids_var_x <- poids_var_x[poids_var_global>threshold]
	poids_var_y <-poids_var_y[poids_var_global>threshold]
	vecnames <- vecnames[poids_var_global>threshold]
	# Tracer un fond blanc orthonormé
	plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",asp=1,axes=F)  ; grid()
	# Tracer le cercle
	draw.circle(0,0,1,lwd=3) 
	# Tracer les flèches
	arrows(0,0,poids_var_x,poids_var_y,length=0.1,angle=17,lwd=2)
	# Calculer la position du nom de chaque vecteur // aux pointes des flèches
	poids_x_bin <- poids_var_x ; poids_x_bin[poids_x_bin<0]<- -1 ;  poids_x_bin[poids_x_bin>0]<- 1
	poids_y_bin <- poids_var_y ; poids_y_bin[poids_y_bin<0]<- -1 ;  poids_y_bin[poids_y_bin>0]<- 1
	poids_x <- poids_var_x^2 ; poids_y <- poids_var_y^2; poids_z <- poids_x + poids_y
	poids_x <- poids_x / poids_z * poids_x_bin * 0.08 ; poids_y <- poids_y / poids_z * poids_y_bin * 0.08
	# Tracer les noms des vecteurs
	text(poids_var_x+poids_x,poids_var_y+poids_y,vecnames,cex=0.75,col="red")
	par(new = T)
	# Graphique issu de l'étape 4 - données de prcomp()
	if (xlab=="cps") {xlab=paste("Composante",choices[1])}
	if (ylab=="cps") {ylab=paste("Composante",choices[2])}
	axeX <- pca$x[,choices[1]] ; axeY <- pca$x[,choices[2]] 
	plot(axeX,axeY,pch=pch,cex=cex,col=col,xlab=xlab,ylab=ylab);grid()
	text(axeX,axeY-5,labels)
	# Pour superposition des 2 graphiques
	par(new = T)
	# Graphique issu de l'étape 5 
	plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",asp=1,axes=F) 
	arrows(0,0,poids_var_x,poids_var_y,length=0.1,angle=17,lwd=2)
	draw.circle(0,0,1,lwd=3)
	text(poids_var_x+poids_x,poids_var_y+poids_y,vecnames,cex=0.75,col="red")
}










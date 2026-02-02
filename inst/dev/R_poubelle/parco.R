require("plotly")
require("htmlwidgets")
#' Parallel coordinates in a single line of code
#'
#' @param data A data.frame
#' @param Y The Y name of number (for coloration)
#' @param X The X list (optional)
#' @param save For saving the graph (html format)
#' @param file Name of the html page "xxx.html"
#'
#' @return This function allows to explore a large dataset by the parallel coordinates method, using a simplified approach of the plot_ly approach.
#' @import htmlwidgets
#' @importFrom plotly plot_ly
#' @export
#'
#' @examples
#' parco(iris,5)
#' parco(iris,1,2:4)
#' parco(iris[,1:4],"Sepal.Length")
parco <- function(data,Y=c(),X=c(),save=F,file="file.html") {
  colnames(data)<-gsub("-",".",colnames(data))
  # data : une data.frame
  # Y : un numéro de colonne ou un titre
  # X : un vecteur contenant ou des numéros de colonnes ou des titres
  if (is.numeric(Y) == T) {
    mon_titre = paste("Variable :",colnames(data)[Y])
    j <- colnames(data)[Y]
  }else {
    Y<-gsub("-",".",Y)
    mon_titre = paste("Variable :",Y)
    j <- Y
    Y <- which(colnames(data)%in%j)
  }
  if (length(X)==0) {X <- 1:ncol(data)}

  if (is.numeric(X) == T) {#do nothing
  }else {
    which(colnames(data)%in% X) -> X
  }
  if (length(X) !=  length(X[sapply(data[,X],is.numeric)]) ){
    cat("Warning ! Some X are characters.\n")
    X <- which(sapply(data[,X],is.numeric))
  }
  dims <- list() ; k<-0
  for (i in X) {
    k<-k+1
    dims[[k]]<- list(range = range(data[,i],na.rm=T),
                     label = colnames(data)[i],values = formula(paste0("~",colnames(data)[i])))
  }
  # Cette version a été réalisée par Antoine Massé
  # Site  :Aide à l'utilisation de R
  # Merci d'en indiquer la source si vous la recopiée
  # parallele coordinate plot sur variable
  if (is.numeric(data[,Y])) {
    fig <- data %>% plot_ly(type = 'parcoords',
                            line = list(color = ~get(j),
                                        colorscale = 'Jet',
                                        showscale = TRUE,
                                        reversescale = TRUE,
                                        cmin = min( data[[j]]),
                                        cmax = max( data[[j]])),
                            dimensions = dims) %>% layout(title = mon_titre)
  } else {
    # Conversion des catégories en valeurs comprises entre 0 et 1
    categories <- as.character(data[,Y]) ; k <- 0 ; colorscale <- list() ;
    colors <- rainbow(length(categories))
    for (i in unique(categories)) {
      categories[categories==i] <- k
      colorscale[[k+1]] <- c(k,colors[k+1])
      k<-k+1
    }
    categories <- as.numeric(categories)
    fig <- data %>% plot_ly(type = 'parcoords',
                            line = list(color = categories,
                                        colorscale = colorscale),
                            #		showscale = TRUE,
                            #		reversescale = TRUE,
                            #		cmin = min( data[[j]]),
                            #		cmax = max( data[[j]])),
                            dimensions = dims) %>% layout(title = mon_titre)
  }
  if (save==F) {print(fig)}
  if (save==T) {saveWidget(fig, file, selfcontained = F, libdir = "lib")}
}

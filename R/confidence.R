#' Confidence interval of a sample
#'
#' @param vector : A vector
#' @param conf.level : Confidence level (from 0 to 1)
#' @param na.rm : Exclusion of missing values
#'
#' @return Confidence interval on a sample.
#' @export
#'
#' @examples
#' x = c(1,4,3,5,3) ; int.ech(x)
int.ech = function(vector,conf.level=0.95,na.rm=T) { # VERSION 2016
  if (length(vector)==0) { cat("Erreur ! Le vecteur ",substitute(vector),"est vide.\n")}
  else { s = var(vector,na.rm=na.rm)
  n = length(vector)-sum(is.na(vector))
  ddl = n - 1 ; proba = (1-conf.level)*100 ; proba = (100-proba/2)/100
  t_student = qt(proba, ddl) # quantile
  intervalle = t_student * sqrt(s/n)
  moyenne = mean(vector,na.rm=na.rm) ; return(intervalle) }}


#' Confidence interval of a population
#'
#' @param vector : A vector
#' @param conf.level : Confidence level (from 0 to 1)
#' @param sigma : If known : standard deviation of the population
#'
#' @return Confidence interval on a population.
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' int.pop(x,0.95)
int.pop =function(vector,conf.level=0.95 ,sigma=c() ) {
  # vector : un échantillon
  # sigma : écart-type de la population
  if (length(sigma)==0) {sigma=sd(vector)}
  n = length(vector)
  proba =(100-(((1-conf.level)*100)/2))/100
  z = qnorm(proba) ; # Récupération des quantiles de la loi normale ==> Par exemple : si on veut englober 95% des valeurs : estimer la variabilité de la population en prenant le risque de négliger 5% de celle-ci, spot 2,5% en tête de distribution
  intervalle = z*sigma/sqrt(n)
  moyenne = mean(vector)
  return(intervalle) }


#' Confidence interval of a proportion
#'
#' @param proportion : One or more proportions between 0 and 1.
#' @param n : Number of individuals who determined these proportions
#' @param conf.level : Confidence level (from 0 to 1)
#'
#' @return Confidence interval on a proportion.
#' @export
#'
#' @examples
#' int.prop(0.50,100,0.95)
#' #Will give the confidence interval on a proportion of 50% (0.5) calculated on a sample of 100 items.
int.prop = function(proportion,n,conf.level=0.95) {
  # Il faut donner une proportion entre 0 et 1 et une taille d'échantillon
  # proportion peut aussi correspondre à une liste de proportions
  nombre_proportions = length(proportion)
  P = proportion
  proba = (1-((1-conf.level)/2))
  qnorm(proba)
  IC = c(qnorm(proba)*sqrt(P*(1-P)/n))
  return(IC) }


#' Confidence intervals on a matrix of a contingent
#'
#' @param x : a matrix (a two-dimensional contingency table: the entries of x must be non-negative integers).
#' @param margin : a vector giving the margins to split by. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. When x has named dimnames, it can be a character vector selecting dimension names.
#' @param conf.level : Confidence level (from 0 to 1)
#'
#' @return Confidence interval on proportions of table obtained with prop.table().
#' @export
#'
#' @examples
#' contingent = with(airquality, table(cut(Temp, quantile(Temp)), Month))
#' prop.table(contingent) # proportions
#' int.prop.table(contingent) # confidence intervals
#' # by row
#' prop.table(contingent,1) # proportions
#' int.prop.table(contingent,1) # confidence intervals
#' # by col
#' prop.table(contingent,2) # proportions
#' int.prop.table(contingent,2) # confidence intervals

int.prop.table = function(x,margin=NULL,conf.level=0.95) {
    if (is.null(margin)) {
      myIC <- int.prop(prop.table(x),sum(x),conf.level=conf.level)
      myIC <- matrix(myIC,nc=ncol(x),nr=nrow(x))
      myIC <- matrix(myIC,nc=ncol(x),nr=nrow(x))
    } else if (margin==1) {
      contingent <- prop.table(x,1)
      effectifs <- apply(x,1,sum)
      myIC <- c()
      for (i in 1:length(effectifs)) {
        myIC <- c(myIC,int.prop(contingent[i,],effectifs[i],conf.level=conf.level))
      }
      myIC <- matrix(myIC,nc=ncol(x),nr=nrow(x),byrow=TRUE)
    } else if (margin==2) {
      contingent <- prop.table(x,2)
      effectifs <- apply(x,2,sum)
      myIC <- c()
      for (i in 1:length(effectifs)) {
        myIC <- c(myIC,int.prop(contingent[,i],effectifs[i],conf.level=conf.level))
      }
      myIC <- matrix(myIC,nc=ncol(x),nr=nrow(x))
    }
    return(myIC)
  }

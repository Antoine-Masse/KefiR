#' Dynamic expansion on a vector (changing the mean without necessarily changing the range of values)
#'
#' @param liste A vector
#' @param moyenne The new desired average for the vector
#' @param max The new max value desired for the vector
#' @param round Define the rounding level
#'
#' @return Dynamic expansion on a vector (changing the mean without necessarily changing the range of values)
#' @export
#'
#' @examples
#' x <- round(rnorm(50,10,4),0) # evaluation of students
#' hist(x,main="before")
#' exp_dyn(x,12,20)
#' hist(x,main="after")
exp_dyn <- function(liste,moyenne,max,round=0) {
  moy_liste <- (mean(liste,na.rm=T)+median(liste,na.rm=T))/2
  min <- min(liste,na.rm=T)
  liste_temp <- liste
  liste[(liste_temp<=moy_liste)&!is.na(liste)] <- (liste[liste_temp<=moy_liste&!is.na(liste)]-min)/max((liste[liste_temp<=moy_liste&!is.na(liste)]-min),na.rm=T)*(moyenne-min)+min
  #liste[liste_temp<moy_liste&!is.na(liste)] <- (liste[liste_temp<moy_liste&!is.na(liste)]-min)/(moy_list-min)*(moyenne-min)+min
  liste[liste_temp==moy_liste&!is.na(liste)] <- moyenne
  liste[liste_temp>moy_liste&!is.na(liste)] <- (liste[liste_temp>moy_liste&!is.na(liste)]-moy_liste)/max((liste[liste_temp>=moy_liste&!is.na(liste)]-moy_liste),na.rm=T)*(max-moyenne)+moyenne
  liste <- round(liste,round)
  return(liste)
}

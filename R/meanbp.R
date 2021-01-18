#' Moving average (weighted) per iteration
#'
#' @param vector : a vector
#' @param it : numbre of iterations
#' @param na.rm : for ignore missing values
#'
#' @return This function returns an average whose values are weighted according to the distance that separates them from the median initially, and then, at each cycle, with respect to the "average" calculated in the previous cycle.
#' @export
#'
#' @examples
#' x <- rnorm(100,20,5) # simulate a sample of a population with an average of 20%.
#' meanbp(x,100) # 100 because the calculation will only be repeated a maximum of 100 times.
#' #if we compare with mean(x) or median(x): we see that meanbp gives a result closer to the true average which is 20
meanbp = function(vector, it=1000,na.rm=T) {
  #Method Mosteller et Tukey, 1977 :
  # Antoine Masse
  # Version 2021
  if (na.rm==T) {vector <- vector[!is.na(vector)]}
  if (length(vector) >2) {
    EIQ = abs(quantile(vector,probs=0.75)-quantile(vector,probs=0.25))
    if (EIQ ==0) {moyenne = mean(vector)}
    else {
      Z = (vector - median(vector))/(3*EIQ)
      for (i in c(1:length(vector))) {
        if (abs(Z[i])>1) {Z[i]=1}        }
      w = (1-Z^2)^2
      moyenne = sum(vector*w)/sum(w)
      #Methode complement A.M. : Allows you to weight the values according to the distance between them and the previously calculated average (1000 cycles maximum).
      for (j in c(1:it)) {
        moyenne_old = moyenne
        d = abs(vector-moyenne)
        EIQ = (quantile(d,probs=0.75))+(quantile(d,probs=0.25))
        Z = (vector - moyenne)/(3*EIQ)
        for (i in c(1:length(vector))) {
          if (abs(Z[i])>1) {Z[i]=1}
        }
        w = (1-Z^2)^2
        moyenne = sum(vector*w)/sum(w)
        variation = sd(vector)/(moyenne_old-moyenne)
        if (moyenne ==moyenne_old) { break } # si moyenne ne bouge plus
        if (variation > 1000) { break} # si le signal sur bruit est élevé
      }}   }
  else {moyenne = mean(vector)}
  return(moyenne)}
